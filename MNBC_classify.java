import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.GZIPInputStream;

import org.eclipse.collections.api.iterator.IntIterator;
import org.eclipse.collections.api.iterator.MutableIntIterator;
import org.eclipse.collections.api.list.primitive.MutableIntList;
import org.eclipse.collections.api.set.primitive.MutableIntSet;
import org.eclipse.collections.impl.list.mutable.primitive.IntArrayList;
import org.eclipse.collections.impl.set.mutable.primitive.IntHashSet;

public class MNBC_classify { //Previously called MNBC_classify2_onlydelta1000
	private static int k;
	private static int numberOfThreads;
	private static float kmerPenalty = -2000.0F;
	private static float delta = 1500.0F;
	private static String dbDirPath;
	private static String metaFilePath; //Example filename: refSeq_prokaryote_complete_genomes_ok_status_metainfo_300k.txt
	private static String outputFilePath;
	private static boolean readType; //Whether reads are paired-end
	private static String startPath;
	private static String endPath;
	
	private static String[] genomeIds;
	private static float[] logFres;
	private static MutableIntSet[] genomeMinimizers;
	private static HashMap<String, String[]> completeGenomeId2TaxIds;	
	private static HashSet<String> finishedReadIds;
	private static BlockingQueue<String[]> readQueue; //Balance producer and consumers
	private static BlockingQueue<String> resultQueue; //Balance consumers and writer
	private static AtomicInteger erroredConsumerCount = new AtomicInteger();
	
	public static void main(String[] args) {
		for(int i = 0; i < args.length; i++) {
			if(args[i].startsWith("-")) {
				switch(args[i].charAt(1)) {
					case 'k':
						k = Integer.parseInt(args[i + 1]);
						break;
					case 'p':
						kmerPenalty = Float.parseFloat(args[i + 1]);
						break;
					case 'e':
						delta = Float.parseFloat(args[i + 1]);
						break;
					case 'c':
						numberOfThreads = Integer.parseInt(args[i + 1]);
						break;
					case 'd':
						dbDirPath = args[i + 1];
						break;
					case 'm':
						metaFilePath = args[i + 1];
						break;
					case 'o':
						outputFilePath = args[i + 1];
						break;
					case 't':
						readType = args[i + 1].equals("2") ? true : false; //The parameter value itself is "1" or "2"
						startPath = args[i + 2];
						if(readType) {
							endPath = args[i + 3];
						}
						break;
					case 'h':
						printHelpInfo();
						System.exit(0);
				}
			}
		}
		
		long startTime = System.nanoTime();		
		int numberOfCores = Runtime.getRuntime().availableProcessors();
		System.out.println("Number of available cores: " + numberOfCores);
		if(numberOfThreads > numberOfCores) {
			System.out.println("WARNING - Number of available cores " + numberOfCores + " is less than requested number of threads " + numberOfThreads + ", exiting");
			System.exit(1);
		}
		readQueue = new ArrayBlockingQueue<String[]>(numberOfThreads + 10);
		resultQueue = new ArrayBlockingQueue<String>(numberOfThreads + 10);
		
		File outputFile = new File(outputFilePath);
		if(outputFile.exists()) {			
			readBaseOutputFile(outputFile);
			System.out.println(finishedReadIds.size() + " reads have finished previously");
		}
		
		File[] countFiles = new File(dbDirPath).listFiles();
		genomeIds = new String[countFiles.length];
		logFres = new float[countFiles.length];
		genomeMinimizers = new MutableIntSet[countFiles.length];
		
		ExecutorService nested = Executors.newFixedThreadPool(numberOfCores - 1);
		CompletionService<String> pool = new ExecutorCompletionService<String>(nested);
		for(int i = 0; i < countFiles.length; i++) {
			pool.submit(new DBReader(countFiles[i], i));
		}
		
		for(int i = 0; i < countFiles.length; i++) {
			try {
				String outcome = pool.take().get();
				if(outcome.contains("ERROR")) {
					System.out.println(i + "th task failed (" + outcome + "), exiting");
					System.exit(1);
				}
			} catch(Exception e) {
				System.out.println("Exception on " + i + " th returned task, exiting");
				e.printStackTrace();
				System.exit(1);
			}
		}		
		nested.shutdown();
		long endTime = System.nanoTime();
		System.out.println("Read DB in " + + ((endTime - startTime) / 1000000000) + " seconds");
		
		completeGenomeId2TaxIds = readCompleteMeta(metaFilePath);
		
		new Thread(new Producer()).start();
		
		for(int i = 0; i < numberOfThreads; i++) {
			new Thread(new Consumer(i)).start();
		}
		System.out.println("Start classifying");
		
		int completedConsumerCounter = 0;
		try {
			PrintWriter writer = null;
			if(finishedReadIds == null) {
				writer = new PrintWriter(new FileWriter(outputFilePath), true);
				writer.println("Read\tGenome\tSpecies\tGenus\tFamily\tOrder\tClass\tPhylum\tSuperkingdom");
			} else {
				writer = new PrintWriter(new FileWriter(outputFilePath, true), true);
			}			
			
			while((completedConsumerCounter + erroredConsumerCount.get()) < numberOfThreads) {
				String outcome = resultQueue.take();
				if(outcome.endsWith("- finished")) {
					completedConsumerCounter++;
				} else {
					writer.println(outcome);
				}
			}
			writer.close();
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		endTime = System.nanoTime();
		System.out.println("done in " + ((endTime - startTime) / 1000000000) + " seconds");
	}
	
	private static void readBaseOutputFile(File outputFile) {
		finishedReadIds = new HashSet<String>();
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader(outputFile));
			String line = reader.readLine();
			while((line = reader.readLine()) != null) {
				String[] fields = line.split("\t");
				if(fields.length == 9) {
					finishedReadIds.add(fields[0]);
				}
			}
			reader.close();
		} catch(Exception e) {
			System.out.println("ERROR - unable to read previously finished reads");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	private static HashMap<String, String[]> readCompleteMeta(String completeMetaFilePath) {
		HashMap<String, String[]> completeGenomeId2TaxIds = new HashMap<String, String[]>();
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader(completeMetaFilePath));
			String line = reader.readLine();
			while((line = reader.readLine()) != null) {
				String[] fields = line.split("\t");
				completeGenomeId2TaxIds.put(fields[0], new String[] {fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], fields[7]});
			}
			reader.close();
		} catch(Exception e) {
			System.out.println("ERROR: failed reading taxonomy file");
			e.printStackTrace();
			System.exit(1);
		}
		
		return completeGenomeId2TaxIds;
	}
	
	private static class Consumer implements Runnable {
		private int id;
		private String[] read;
		
		public Consumer(int anID) {
			id = anID;
		}

		@Override
		public void run() {
			boolean done = false;			
			
			try {
				if(readType) {
					while(!done) {
						read = readQueue.take();
						if(read.length == 0) {
							readQueue.put(read);						
							resultQueue.put("Consumer " + id + " - finished");
							done = true;
						} else {
							MutableIntSet readMinimizers = new IntHashSet();							
							addKmersOfOneTestFrag(read[1], readMinimizers);
							addKmersOfOneTestFrag(read[2], readMinimizers);
							if(readMinimizers.isEmpty()) {
								resultQueue.put(read[0] + "\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull");
								continue;
							}
							processReadMinimizers(readMinimizers);
						}
					}
				} else {
					while(!done) {
						read = readQueue.take();
						if(read.length == 0) {
							readQueue.put(read);						
							resultQueue.put("Consumer " + id + " - finished");
							done = true;
						} else {
							MutableIntSet readMinimizers = new IntHashSet();							
							addKmersOfOneTestFrag(read[1], readMinimizers);
							if(readMinimizers.isEmpty()) {
								resultQueue.put(read[0] + "\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull");
								continue;
							}							
							processReadMinimizers(readMinimizers);
						}
					}
				}
			} catch(Exception e) {
				System.out.println("Consumer " + id + " - error occurred on read " + read[0] + ", ending");
				e.printStackTrace();
				erroredConsumerCount.incrementAndGet();
			}
		}
		
		private void processReadMinimizers(MutableIntSet readMinimizers) throws Exception {
			TreeMap<Float, MutableIntList> topScores = new TreeMap<>();
			for(int i = 0; i < genomeIds.length; i++) {							
				float score = 0.0F;

				IntIterator it = readMinimizers.intIterator();
				while(it.hasNext()) {
					if(genomeMinimizers[i].contains(it.next())) {
						score += logFres[i];
					} else {
						score += kmerPenalty;
					}
				}

				if(topScores.containsKey(score)) {
					topScores.get(score).add(i);
				} else {
					MutableIntList genomeIdsWithScore = new IntArrayList();
					genomeIdsWithScore.add(i);
					topScores.put(score, genomeIdsWithScore);									
				}
			}
			
			MutableIntList votingGenomes = processTopScores(topScores);							
			String outcome = read[0];															
			if(votingGenomes.size() == 1) {
				//System.out.println("Read " + read[0] + " has 1 voting genome");
				String predictedGenomeId = genomeIds[votingGenomes.getFirst()];
				String[] predictedTaxonIds = completeGenomeId2TaxIds.get(predictedGenomeId);
				outcome += "\t" + predictedGenomeId;
				for(String predictedId : predictedTaxonIds) {
					outcome += "\t" + predictedId;
				}								
			} else {
				//System.out.println("Read " + read[0] + " has " + votingGenomes.size() + " voting genomes");
				HashMap<String, ArrayList<String>> speciesId2GenomeIds = fillSpeciesId2GenomeIds(votingGenomes);								
				String dominantSpecies = null; //Multiple species may have the same greatest number of genomes
				int dominantCount = 0;
				for(Entry<String, ArrayList<String>> species : speciesId2GenomeIds.entrySet()) {
					int count = species.getValue().size();
					if(count >= dominantCount) {
						dominantSpecies = species.getKey();
						dominantCount = count;
					}
				}
				//System.out.println("Dominant species is " + dominantSpecies + " with " + dominantCount + " genomes");
				
				outcome += "\tnull";
				String[] taxonIds = completeGenomeId2TaxIds.get(speciesId2GenomeIds.get(dominantSpecies).get(0));
				for(String taxonId : taxonIds) {
					outcome += "\t" + taxonId;
				}
			}
			resultQueue.put(outcome);
		}
		
		private MutableIntList processTopScores(TreeMap<Float, MutableIntList> topScores) {
			MutableIntList votingGenomes = new IntArrayList();
			Entry<Float, MutableIntList> greatestEntry = topScores.pollLastEntry();
			votingGenomes.addAll(greatestEntry.getValue());
			
			float prevScore = greatestEntry.getKey();
			while(!topScores.isEmpty()) {
				Entry<Float, MutableIntList> curEntry = topScores.pollLastEntry();
				float curScore = curEntry.getKey();
				if(prevScore - curScore > delta) {
					break;
				} else {
					votingGenomes.addAll(curEntry.getValue());
					prevScore = curScore;
				}
			}
			
			return votingGenomes;
		}
		
		private HashMap<String, ArrayList<String>> fillSpeciesId2GenomeIds(MutableIntList votingGenomes) {
			HashMap<String, ArrayList<String>> speciesId2GenomeIds = new HashMap<>();
			MutableIntIterator it = votingGenomes.intIterator();
			while(it.hasNext()) {
				String genomeId = genomeIds[it.next()];
				String speciesId = completeGenomeId2TaxIds.get(genomeId)[0];
				if(speciesId2GenomeIds.containsKey(speciesId)) {
					speciesId2GenomeIds.get(speciesId).add(genomeId);
				} else {
					ArrayList<String> genomeIdsWithSpeciesId = new ArrayList<String>();
					genomeIdsWithSpeciesId.add(genomeId);
					speciesId2GenomeIds.put(speciesId, genomeIdsWithSpeciesId);
				}									
			}
			
			return speciesId2GenomeIds;
		}
		
		private void addKmersOfOneTestFrag(String testFrag, MutableIntSet kmers) {
			int length = testFrag.length();
			if(length < k) {
				return;
			}
			
			//Get minus sequence
			MutableIntSet indicesOfInvalidKmers = new IntHashSet();
			StringBuilder minusSequence = new StringBuilder();			
			for(int i = length - 1; i >= 0; i--) {
				char base = testFrag.charAt(i);
				switch(base) {
					case 'A':
						minusSequence.append('T');
						break;
					case 'C':
						minusSequence.append('G');
						break;
					case 'G':
						minusSequence.append('C');
						break;
					case 'T':
						minusSequence.append('A');
						break;
					default:
						minusSequence.append(base);
						for(int j = i - k + 1; j <= i; j++) {
							indicesOfInvalidKmers.add(j);
						}
				}
			}			
			
			HashSet<String> minimizers = new HashSet<String>();
			int startIndexOfLastKmer = length - k;
			String[] canonicalKmers = new String[startIndexOfLastKmer + 1]; //first get the canonical kmer at each position
			
			if(indicesOfInvalidKmers.isEmpty()) { //all kmers are valid
				for(int i = 0; i <= startIndexOfLastKmer; i++) {
					String plusKmer = testFrag.substring(i, i + k);
					String minusKmer = minusSequence.substring(startIndexOfLastKmer - i, length - i);
					canonicalKmers[i] = (plusKmer.compareTo(minusKmer) < 0) ? plusKmer : minusKmer;
				}
				
				if(length < (2 * k - 1)) { //read length not reaching one window's length
					//left end minimizers, dynamic programming					
					String minimizer = canonicalKmers[0];
					minimizers.add(minimizer);
					for(int v = 1; v < canonicalKmers.length; v++) { //here v is not the u in (u,k) in the paper, it is the index of the last/rightmost kmer in each (u,k) window (i.e. u - 1)
						minimizer = (canonicalKmers[v].compareTo(minimizer) < 0) ? canonicalKmers[v] : minimizer;
						minimizers.add(minimizer);					
					}

					for(String aMinimizer : minimizers) {
						int index = 0;
						
						int lastCharIndex = k - 1;
						for(int i = lastCharIndex; i >= 0; i--) {
							switch(aMinimizer.charAt(i)) {
								case 'C':
									index += (int) Math.pow(4, lastCharIndex - i);
									break;
								case 'G':
									index += 2 * ((int) Math.pow(4, lastCharIndex - i));
									break;
								case 'T':
									index += 3 * ((int) Math.pow(4, lastCharIndex - i));
									break;
							}
						}
						
						kmers.add(index);
					}
					return;
				}
				
				//interior minimizers - first window
				String minimizer = canonicalKmers[0];
				for(int j = 1; j < k; j++) {
					minimizer = (canonicalKmers[j].compareTo(minimizer) < 0) ? canonicalKmers[j] : minimizer;
				}
				minimizers.add(minimizer);

				//interior minimizers - later windows incrementally
				int indexOfLastWindowStartPosition = startIndexOfLastKmer - k + 1;					
				for(int i = 1; i <= indexOfLastWindowStartPosition; i++) {
					if(minimizer.equals(canonicalKmers[i - 1])) {
						minimizer = canonicalKmers[i];
						int indexOfFirstKmerAfterWindow = i + k;
						for(int j = i + 1; j < indexOfFirstKmerAfterWindow; j++) {
							minimizer = (canonicalKmers[j].compareTo(minimizer) < 0) ? canonicalKmers[j] : minimizer;
						}
					} else {
						int indexOfLastKmerInWindow = i + k - 1;
						minimizer = (canonicalKmers[indexOfLastKmerInWindow].compareTo(minimizer) < 0) ? canonicalKmers[indexOfLastKmerInWindow] : minimizer;
					}
					
					minimizers.add(minimizer);					
				}

				//left end minimizers, dynamic programming
				minimizer = canonicalKmers[0];
				minimizers.add(minimizer);
				for(int v = 1; v <= (k - 2); v++) { //here v is not the u in (u,k) in the paper, it is the index of the last/rightmost kmer in each (u,k) window (i.e. u - 1)
					minimizer = (canonicalKmers[v].compareTo(minimizer) < 0) ? canonicalKmers[v] : minimizer;
					minimizers.add(minimizer);					
				}

				//right end minimizers, dynamic programming
				minimizer = canonicalKmers[startIndexOfLastKmer];
				minimizers.add(minimizer);
				for(int v = startIndexOfLastKmer - 1; v >= (startIndexOfLastKmer - k + 2); v--) { //here v is not the u in (u,k) in the paper, it is the index of the first/leftmost kmer in each (u,k) window
					minimizer = (canonicalKmers[v].compareTo(minimizer) < 0) ? canonicalKmers[v] : minimizer;
					minimizers.add(minimizer);					
				}
			} else { //there are invalid kmers
				for(int i = 0; i <= startIndexOfLastKmer; i++) {
					if(!indicesOfInvalidKmers.contains(i)) {
						String plusKmer = testFrag.substring(i, i + k);
						String minusKmer = minusSequence.substring(startIndexOfLastKmer - i, length - i);
						canonicalKmers[i] = (plusKmer.compareTo(minusKmer) < 0) ? plusKmer : minusKmer;
					}
				}
				
				if(length < (2 * k - 1)) { //read length not reaching one window's length
					//left end minimizers, dynamic programming
					String minimizer = "Z";
					for(int v = 0; v < canonicalKmers.length; v++) { //here v is not the u in (u,k) in the paper, it is the index of the last/rightmost kmer in each (u,k) window (i.e. u - 1)
						if(canonicalKmers[v] != null) {
							minimizer = (canonicalKmers[v].compareTo(minimizer) < 0) ? canonicalKmers[v] : minimizer;
							minimizers.add(minimizer);
						}
					}
					
					for(String aMinimizer : minimizers) {
						int index = 0;
						
						int lastCharIndex = k - 1;
						for(int i = lastCharIndex; i >= 0; i--) {
							switch(aMinimizer.charAt(i)) {
								case 'C':
									index += (int) Math.pow(4, lastCharIndex - i);
									break;
								case 'G':
									index += 2 * ((int) Math.pow(4, lastCharIndex - i));
									break;
								case 'T':
									index += 3 * ((int) Math.pow(4, lastCharIndex - i));
									break;
							}
						}
						
						kmers.add(index);
					}
					return;			
				}
				
				//interior minimizers - first window
				String minimizer = "Z";
				for(int j = 0; j < k; j++) {
					if(canonicalKmers[j] != null) {
						minimizer = (canonicalKmers[j].compareTo(minimizer) < 0) ? canonicalKmers[j] : minimizer;
					}
				}
				if(minimizer.length() == k) {
					minimizers.add(minimizer);
				}

				//interior minimizers - later windows incrementally
				int indexOfLastWindowStartPosition = startIndexOfLastKmer - k + 1;					
				for(int i = 1; i <= indexOfLastWindowStartPosition; i++) {
					if(minimizer.equals(canonicalKmers[i - 1])) {
						minimizer = "Z";
						int indexOfFirstKmerAfterWindow = i + k;
						for(int j = i; j < indexOfFirstKmerAfterWindow; j++) {
							if(canonicalKmers[j] != null) {
								minimizer = (canonicalKmers[j].compareTo(minimizer) < 0) ? canonicalKmers[j] : minimizer;
							}
						}
					} else {
						int indexOfLastKmerInWindow = i + k - 1;
						if(canonicalKmers[indexOfLastKmerInWindow] != null) {
							minimizer = (canonicalKmers[indexOfLastKmerInWindow].compareTo(minimizer) < 0) ? canonicalKmers[indexOfLastKmerInWindow] : minimizer;
						}
					}
					
					if(minimizer.length() == k) {
						minimizers.add(minimizer);
					}
				}

				//left end minimizers, dynamic programming
				minimizer = "Z";
				for(int v = 0; v <= (k - 2); v++) { //here v is not the u in (u,k) in the paper, it is the index of the last/rightmost kmer in each (u,k) window (i.e. u - 1)
					if(canonicalKmers[v] != null) {
						minimizer = (canonicalKmers[v].compareTo(minimizer) < 0) ? canonicalKmers[v] : minimizer;
						minimizers.add(minimizer);
					}
				}

				//right end minimizers, dynamic programming
				minimizer = "Z";
				for(int v = startIndexOfLastKmer; v >= (startIndexOfLastKmer - k + 2); v--) { //here v is not the u in (u,k) in the paper, it is the index of the first/leftmost kmer in each (u,k) window
					if(canonicalKmers[v] != null) {
						minimizer = (canonicalKmers[v].compareTo(minimizer) < 0) ? canonicalKmers[v] : minimizer;
						minimizers.add(minimizer);
					}
				}
			}
			
			for(String aMinimizer : minimizers) {
				int index = 0;
				
				int lastCharIndex = k - 1;
				for(int i = lastCharIndex; i >= 0; i--) {
					switch(aMinimizer.charAt(i)) {
						case 'C':
							index += (int) Math.pow(4, lastCharIndex - i);
							break;
						case 'G':
							index += 2 * ((int) Math.pow(4, lastCharIndex - i));
							break;
						case 'T':
							index += 3 * ((int) Math.pow(4, lastCharIndex - i));
							break;
					}
				}
				
				kmers.add(index);
			}
		}
	}
	
	private static class Producer implements Runnable {		
		private int readCounter;
		
		@Override
		public void run() {
			try {
				BufferedReader reader1 = null;
				if(startPath.endsWith(".gz")) {
					reader1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(startPath)), "UTF-8"));
				} else {
					reader1 = new BufferedReader(new FileReader(startPath));
				}
				
				if(readType) {
					BufferedReader reader2 = null;
					if(startPath.endsWith(".gz")) {
						reader2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(endPath)), "UTF-8"));
					} else {
						reader2 = new BufferedReader(new FileReader(endPath));
					}
					
					if(startPath.contains(".fasta")) {
						System.out.println("Start reading paired-end FASTA files " + startPath + " and " + endPath);
						readTestFragFastaFiles(reader1, reader2);
					} else {
						System.out.println("Start reading paired-end FASTQ files " + startPath + " and " + endPath);
						readTestFragFastqFiles(reader1, reader2);
					}					
				} else {				
					if(startPath.contains(".fasta")) {
						System.out.println("Start reading FASTA file " + startPath);
						readTestFragFastaFile(reader1);
					} else {
						System.out.println("Start reading FASTQ file " + startPath);
						readTestFragFastqFile(reader1);
					}					
				}
				
				System.out.println("Finished " + readCounter + " reads");
			} catch(Exception e) {
				System.out.println("Producer - error occurred, exiting");
				e.printStackTrace();
				System.exit(1);
			}			
		}
		
		private void readTestFragFastaFile(BufferedReader reader) throws Exception {
			String line = null;
			String readId = null;
			StringBuilder sequence = new StringBuilder();
			
			if(finishedReadIds == null) {
				while((line = reader.readLine()) != null) {
					line = line.trim();
					if(line.startsWith(">") || line.startsWith("@")) {						
						if(!sequence.isEmpty()) {
							readQueue.put(new String[] {readId, sequence.toString().toUpperCase()});
							readCounter++;
						}
						
						readId = line.substring(1).split("\\s+")[0];
						sequence.setLength(0);
					} else {
						sequence.append(line);
					}
				}				
				if(!sequence.isEmpty()) {
					readQueue.put(new String[] {readId, sequence.toString().toUpperCase()});
					readCounter++;
				}
			} else {
				boolean skip = false;
				while((line = reader.readLine()) != null) {
					line = line.trim();
					if(line.startsWith(">") || line.startsWith("@")) {
						if(!sequence.isEmpty()) {
							readQueue.put(new String[] {readId, sequence.toString().toUpperCase()});
							readCounter++;
						}
						
						readId = line.substring(1).split("\\s+")[0];
						sequence.setLength(0);
						if(finishedReadIds.contains(readId)) {
							skip = true;						
						} else {
							skip = false;
						}
					} else {
						if(!skip) {
							sequence.append(line);
						}						
					}
				}
				if(!sequence.isEmpty()) {
					readQueue.put(new String[] {readId, sequence.toString().toUpperCase()});
					readCounter++;
				}
			}			
			
			reader.close();
			readQueue.put(new String[0]);
		}
		
		private void readTestFragFastqFile(BufferedReader reader) throws Exception {
			String line = null;
			if(finishedReadIds == null) {
				while((line = reader.readLine()) != null) {
					line = line.trim();
					if(line.startsWith(">") || line.startsWith("@")) {
						String readId = line.substring(1).split("\\s+")[0];						
						String readSequence = reader.readLine().trim().toUpperCase();
						reader.readLine();
						reader.readLine();
						readQueue.put(new String[] {readId, readSequence});
						readCounter++;						
					}
				}
			} else {
				while((line = reader.readLine()) != null) {
					line = line.trim();
					if(line.startsWith(">") || line.startsWith("@")) {
						String readId = line.substring(1).split("\\s+")[0];
						if(finishedReadIds.contains(readId)) {
							for(int i = 0; i < 3; i++) {
								reader.readLine();
							}						
						} else {
							String readSequence = reader.readLine().trim().toUpperCase();
							reader.readLine();
							reader.readLine();
							readQueue.put(new String[] {readId, readSequence});
							readCounter++;
						}
					}
				}
			}			
			
			reader.close();
			readQueue.put(new String[0]);
		}
		
		private void readTestFragFastaFiles(BufferedReader reader1, BufferedReader reader2) throws Exception {		
			String line1 = null;
			String line2 = null;
			if(finishedReadIds == null) {
				while((line1 = reader1.readLine()) != null) {
					line1 = line1.trim();
					line2 = reader2.readLine().trim();
					if(line1.startsWith(">") || line1.startsWith("@")) {
						if(!line2.startsWith(">") && !line2.startsWith("@")) {
							System.out.println("Paired-end FASTA format error: " + line1 + " | " + line2);
							System.exit(1);
						}

						String readId = line1.substring(1).split("\\s+")[0];
						if(!readId.equals(line2.substring(1).split("\\s+")[0])) {
							System.out.println("Paired-end FASTA format error: " + line1 + " | " + line2);
							System.exit(1);
						}
						
						readQueue.put(new String[] {readId, reader1.readLine().trim().toUpperCase(), reader2.readLine().trim().toUpperCase()});
						readCounter++;						
					}
				}
			} else {
				while((line1 = reader1.readLine()) != null) {
					line1 = line1.trim();
					line2 = reader2.readLine().trim();
					if(line1.startsWith(">") || line1.startsWith("@")) {
						if(!line2.startsWith(">") && !line2.startsWith("@")) {
							System.out.println("Paired-end FASTA format error: " + line1 + " | " + line2);
							System.exit(1);
						}

						String readId = line1.substring(1).split("\\s+")[0];
						if(!readId.equals(line2.substring(1).split("\\s+")[0])) {
							System.out.println("Paired-end FASTA format error: " + line1 + " | " + line2);
							System.exit(1);
						}

						if(finishedReadIds.contains(readId)) {
							reader1.readLine();
							reader2.readLine();						
						} else {
							readQueue.put(new String[] {readId, reader1.readLine().trim().toUpperCase(), reader2.readLine().trim().toUpperCase()});
							readCounter++;
						}
					}
				}
			}
			if(reader2.readLine() != null) {
				System.out.println("Paired-end FASTA format error at the end: null | " + line2);
				System.exit(1);
			}
			
			reader1.close();
			reader2.close();
			readQueue.put(new String[0]);
		}
		
		private void readTestFragFastqFiles(BufferedReader reader1, BufferedReader reader2) throws Exception {
			String line1 = null;
			String line2 = null;
			if(finishedReadIds == null) {
				while((line1 = reader1.readLine()) != null) {
					line1 = line1.trim();
					line2 = reader2.readLine().trim();
					if(line1.startsWith(">") || line1.startsWith("@")) {
						if(!line2.startsWith(">") && !line2.startsWith("@")) {
							System.out.println("Paired-end FASTQ format error: " + line1 + " | " + line2);
							System.exit(1);
						}

						String readId = line1.substring(1).split("\\s+")[0];
						if(!readId.equals(line2.substring(1).split("\\s+")[0])) {
							System.out.println("Paired-end FASTQ format error: " + line1 + " | " + line2);
							System.exit(1);
						}
						
						String startSeq = reader1.readLine().trim().toUpperCase();
						String endSeq = reader2.readLine().trim().toUpperCase();					
						reader1.readLine();
						reader1.readLine();
						reader2.readLine();
						reader2.readLine();

						readQueue.put(new String[] {readId, startSeq, endSeq});
						readCounter++;						
					}
				}
			} else {
				while((line1 = reader1.readLine()) != null) {
					line1 = line1.trim();
					line2 = reader2.readLine().trim();
					if(line1.startsWith(">") || line1.startsWith("@")) {
						if(!line2.startsWith(">") && !line2.startsWith("@")) {
							System.out.println("Paired-end FASTQ format error: " + line1 + " | " + line2);
							System.exit(1);
						}

						String readId = line1.substring(1).split("\\s+")[0];
						if(!readId.equals(line2.substring(1).split("\\s+")[0])) {
							System.out.println("Paired-end FASTQ format error: " + line1 + " | " + line2);
							System.exit(1);
						}

						if(finishedReadIds.contains(readId)) {
							for(int i = 0; i < 3; i++) {
								reader1.readLine();
								reader2.readLine();
							}
						} else {
							String startSeq = reader1.readLine().trim().toUpperCase();
							String endSeq = reader2.readLine().trim().toUpperCase();					
							reader1.readLine();
							reader1.readLine();
							reader2.readLine();
							reader2.readLine();

							readQueue.put(new String[] {readId, startSeq, endSeq});
							readCounter++;
						}
					}
				}
			}
			if(reader2.readLine() != null) {
				System.out.println("Paired-end FASTQ format error at the end: null | " + line2);
				System.exit(1);
			}
			
			reader1.close();
			reader2.close();
			readQueue.put(new String[0]);
		}
	}
	
	private static class DBReader implements Callable<String> {
		private File countFile;
		private int id;
		
		public DBReader(File aFile, int anID) {
			countFile = aFile;
			id = anID;
		}

		@Override
		public String call() {
			String filename = countFile.getName();
			String[] fields = filename.split("_");
			genomeIds[id] = fields[0] + "_" + fields[1];
			
			genomeMinimizers[id] = new IntHashSet();
			try {
				BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(countFile)), "UTF-8"));
				String line = reader.readLine();
				logFres[id] = (float) Math.log(1.0 / Integer.parseInt(line));				
				
				while((line = reader.readLine()) != null) {
					genomeMinimizers[id].add(Integer.parseInt(line));
				}
				reader.close();
			} catch(Exception e) {
				e.printStackTrace();
				return "ERROR: couldn't read " + filename;
			}
			
			return "Finished reading " + filename;
		}		
	}
	
	private static void printHelpInfo() {
		System.out.println("This MNBC_classify tool (v1.0) classify reads against a reference database.");
		System.out.println("-h:	Show this help menu");
		System.out.println("-k:	K-mer length");
		System.out.println("-c:	Number of threads");
		System.out.println("-d:	Input database directory");
		System.out.println("-m:	Input taxonomy file");
		System.out.println("-o:	Output classification file");
		System.out.println("-t:	Type of reads (paired-end: 2, single-end: 1). Paired-end reads have two following (gzipped) .fasta/.fastq files. Single-end reads have one following (gzipped) .fasta/.fastq file.");
		System.out.println("-p (optional): Penalty for non-existent k-mers (default -2000)");
		System.out.println("-e (optional): Threshold on the difference between adjacent scores (default 1500)");
	}
}
