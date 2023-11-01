import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.GZIPInputStream;

import org.eclipse.collections.api.iterator.IntIterator;
import org.eclipse.collections.api.set.primitive.MutableIntSet;
import org.eclipse.collections.impl.set.mutable.primitive.IntHashSet;

public class MNBCClassify {
	private static int k;
	private static int numberOfThreads;
	private static float kmerPenalty = -2000.0F;
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
	private static BlockingQueue<String[]> readQueue = new ArrayBlockingQueue<String[]>(150); //Balance producer and consumers
	private static BlockingQueue<String> resultQueue = new ArrayBlockingQueue<String>(150); //Balance producer and consumers
	
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
		
		int numberOfCores = Runtime.getRuntime().availableProcessors();
		System.out.println("Number of available cores: " + numberOfCores);
		if(numberOfThreads > numberOfCores) {
			System.out.println("WARNING - Number of available cores " + numberOfCores + " is less than requested number of threads " + numberOfThreads + ", exiting");
			System.exit(1);
		}
		
		File outputFile = new File(outputFilePath);
		if(outputFile.exists()) {			
			readBaseOutputFile(outputFile);
			System.out.println(finishedReadIds.size() + " reads have finished previously");
		}
		
		long startTime = System.nanoTime();		
		File[] countFiles = new File(dbDirPath).listFiles();
		genomeIds = new String[countFiles.length];
		logFres = new float[countFiles.length];
		genomeMinimizers = new MutableIntSet[countFiles.length];
		
		ExecutorService nested = Executors.newFixedThreadPool(numberOfCores - 1);
		CompletionService<String> pool = new ExecutorCompletionService<String>(nested);
		System.out.println("Created a thread pool");
		for(int i = 0; i < countFiles.length; i++) {
			pool.submit(new DBReader(countFiles[i], i));
		}
		System.out.println("Submitted " + countFiles.length + " DBReader tasks");
		
		for(int i = 0; i < countFiles.length; i++) {
			try {
				System.out.println("Waiting to get outcome of " + i + "th returned task...");
				String outcome = pool.take().get();

				if(outcome.contains("ERROR")) {
					System.out.println("This task failed (" + outcome + "), exiting");
					System.exit(1);
				} else {
					System.out.println("This task succeeded (" + outcome + ")");
				}
			} catch(Exception e) {
				System.out.println("Exception on " + i + " th returned task, exiting");
				e.printStackTrace();
				System.exit(1);
			}
		}		
		nested.shutdown();
		
		completeGenomeId2TaxIds = readCompleteMeta(metaFilePath);
		System.out.println("Read meta file, and mapped all " + completeGenomeId2TaxIds.size() + " reference genomes to their taxon IDs from strain to superkingdom");
		
		new Thread(new Producer()).start();
		System.out.println("Started producer thread");
		
		for(int i = 0; i < numberOfThreads; i++) {
			new Thread(new Consumer(i)).start();
		}
		System.out.println("Started " + numberOfThreads + " consumer threads");
		
		int completedConsumerCounter = 0;
		try {
			PrintWriter writer = null;
			if(finishedReadIds == null) {
				writer = new PrintWriter(new FileWriter(outputFilePath), true);
				writer.println("Read\tGenome\tStrain\tSpecies\tGenus\tFamily\tOrder\tClass\tPhylum\tSuperkingdom\tMaxScore\tGenomeCount");
				System.out.println("Created writer to new result file " + outputFilePath);
			} else {
				writer = new PrintWriter(new FileWriter(outputFilePath, true), true);
				System.out.println("Created writer to existing result file " + outputFilePath);
			}			
			
			while(completedConsumerCounter < numberOfThreads) {
				String outcome = resultQueue.take();
				if(outcome.endsWith("- finished")) {
					completedConsumerCounter++;
					System.out.println(outcome);
				} else {
					writer.println(outcome);
				}
			}
			writer.close();
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		long endTime = System.nanoTime();
		System.out.println("done in " + ((endTime - startTime) / 1000000000) + " seconds");
	}
	
	private static void readBaseOutputFile(File outputFile) {
		finishedReadIds = new HashSet<String>();
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader(outputFile));
			String line = reader.readLine();
			while((line = reader.readLine()) != null) {
				String[] fields = line.split("\t");
				if(fields.length == 12) {
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
				String strainName = "null";
				if(fields[9].contains("strain=")) {
					strainName = fields[9].split("=")[1];
				}
				completeGenomeId2TaxIds.put(fields[0], new String[] {strainName, fields[2], fields[3], fields[4], fields[5], fields[6], fields[7], fields[8]});
			}
			reader.close();
		} catch(Exception e) {
			System.out.println(e.getMessage());
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
							String[] fields = read[1].split("-");
							addKmersOfOneTestFrag(fields[0], readMinimizers);
							addKmersOfOneTestFrag(fields[1], readMinimizers);

							float maxScore = Float.NEGATIVE_INFINITY;
							MutableIntSet genomeIdsWithMaxScore = new IntHashSet();
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

								if(score > maxScore) {
									maxScore = score;
									genomeIdsWithMaxScore.clear();
									genomeIdsWithMaxScore.add(i);
								} else if(score == maxScore) {
									genomeIdsWithMaxScore.add(i);
								}							
							}

							String outcome = read[0] + "\t";
							if(genomeIdsWithMaxScore.size() == 1) {
								String predictedGenomeId = genomeIds[genomeIdsWithMaxScore.intIterator().next()];
								String[] predictedTaxonIds = completeGenomeId2TaxIds.get(predictedGenomeId);
								outcome += predictedGenomeId;
								for(String predictedId : predictedTaxonIds) {
									outcome += "\t" + predictedId;
								}
								outcome += "\t" + maxScore + "\t1";
								resultQueue.put(outcome);
							} else {
								int[] genomeIdsWithMaxScoreArray = genomeIdsWithMaxScore.toArray();
								int predictedLevelOfLCA = getLCALevelOfPredictedGenomes(genomeIdsWithMaxScoreArray, completeGenomeId2TaxIds);
								if(predictedLevelOfLCA == -1) {
									outcome += "null\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull";
								} else {
									outcome += "null";
									for(int j = 0; j < predictedLevelOfLCA; j++) {
										outcome += "\tnull";
									}

									String[] taxonIdsOfOnePredictedGenome = completeGenomeId2TaxIds.get(genomeIds[genomeIdsWithMaxScoreArray[0]]);
									for(int j = predictedLevelOfLCA; j < 8; j++) {
										outcome += "\t" + taxonIdsOfOnePredictedGenome[j];
									}									
								}
								outcome += "\t" + maxScore + "\t" + genomeIdsWithMaxScore.size();
								resultQueue.put(outcome);
							}
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
							
							float maxScore = Float.NEGATIVE_INFINITY;
							MutableIntSet genomeIdsWithMaxScore = new IntHashSet();
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

								if(score > maxScore) {
									maxScore = score;
									genomeIdsWithMaxScore.clear();
									genomeIdsWithMaxScore.add(i);
								} else if(score == maxScore) {
									genomeIdsWithMaxScore.add(i);
								}							
							}

							String outcome = read[0] + "\t";
							if(genomeIdsWithMaxScore.size() == 1) {
								String predictedGenomeId = genomeIds[genomeIdsWithMaxScore.intIterator().next()];
								String[] predictedTaxonIds = completeGenomeId2TaxIds.get(predictedGenomeId);
								outcome += predictedGenomeId;
								for(String predictedId : predictedTaxonIds) {
									outcome += "\t" + predictedId;
								}
								outcome += "\t" + maxScore + "\t1";
								resultQueue.put(outcome);
							} else {
								int[] genomeIdsWithMaxScoreArray = genomeIdsWithMaxScore.toArray();
								int predictedLevelOfLCA = getLCALevelOfPredictedGenomes(genomeIdsWithMaxScoreArray, completeGenomeId2TaxIds);
								if(predictedLevelOfLCA == -1) {
									outcome += "null\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull";
								} else {
									outcome += "null";
									for(int j = 0; j < predictedLevelOfLCA; j++) {
										outcome += "\tnull";
									}

									String[] taxonIdsOfOnePredictedGenome = completeGenomeId2TaxIds.get(genomeIds[genomeIdsWithMaxScoreArray[0]]);
									for(int j = predictedLevelOfLCA; j < 8; j++) {
										outcome += "\t" + taxonIdsOfOnePredictedGenome[j];
									}								
								}
								outcome += "\t" + maxScore + "\t" + genomeIdsWithMaxScore.size();
								resultQueue.put(outcome);
							}
						}
					}
				}
			} catch(Exception e) {
				System.out.println("Consumer " + id + " - error occured on read " + read[0] + ", ending");
				e.printStackTrace();
			}
		}
		
		private int getLCALevelOfPredictedGenomes(int[] genomeIdsWithMaxScoreArray, HashMap<String, String[]> completeGenomeId2TaxIds) {
			String[][] taxonIdsOfAllPredictedGenomes = new String[genomeIdsWithMaxScoreArray.length][8];
			for(int i = 0; i < genomeIdsWithMaxScoreArray.length; i++) {			
				String predictedGenomeId = genomeIds[genomeIdsWithMaxScoreArray[i]];
				taxonIdsOfAllPredictedGenomes[i] = completeGenomeId2TaxIds.get(predictedGenomeId);
			}
			
			for(int col = 0; col < 8; col++) {			
				if(taxonIdsOfAllPredictedGenomes[0][col].equals("null")) {
					continue;
				}
				
				boolean isLCA = true; //If LCA, all predicted genomes have the same id at the level
				for(int row = 1; row < genomeIdsWithMaxScoreArray.length; row++) {
					if(!taxonIdsOfAllPredictedGenomes[row][col].equals(taxonIdsOfAllPredictedGenomes[row - 1][col])) {
						isLCA = false;
						break;
					}
				}
				
				if(isLCA) {
					return col;
				}
			}
			
			return -1;
		}
		
		private void addKmersOfOneTestFrag(String testFrag, MutableIntSet kmers) {
			MutableIntSet indicesOfInvalidKmers = checkIfFragOnlyContainsValidCharacters(testFrag);
			HashSet<String> minimizers = new HashSet<String>();
			int startIndexOfLastKmer = testFrag.length() - k;
			String[] canonicalKmers = new String[startIndexOfLastKmer + 1]; //first get the canonical kmer at each position
			if(indicesOfInvalidKmers.isEmpty()) {
				for(int i = 0; i <= startIndexOfLastKmer; i++) {
					String plusKmer = testFrag.substring(i, i + k);
					String minusKmer = getMinusKmer(plusKmer);
					canonicalKmers[i] = (plusKmer.compareTo(minusKmer) < 0) ? plusKmer : minusKmer;
				}
			} else {
				for(int i = 0; i <= startIndexOfLastKmer; i++) {
					if(!indicesOfInvalidKmers.contains(i)) {
						String plusKmer = testFrag.substring(i, i + k);
						String minusKmer = getMinusKmer(plusKmer);
						canonicalKmers[i] = (plusKmer.compareTo(minusKmer) < 0) ? plusKmer : minusKmer;
					}
				}
			}
			
			//interior minimizers, maybe could be optimized
			int indexOfLastWindowStartPosition = startIndexOfLastKmer - k + 1;					
			for(int i = 0; i <= indexOfLastWindowStartPosition; i++) {
				//System.out.println("*************Window index: " + i + "*************");
				//System.out.println("Window sequence: " + chromosome.substring(i, i + k + k - 1));
				String minimizer = "Z";
				if(canonicalKmers[i] != null) {
					minimizer = canonicalKmers[i];
				}
				int indexOfFirstKmerAfterWindow = i + k;						
				for(int j = i + 1; j < indexOfFirstKmerAfterWindow; j++) {
					if(canonicalKmers[j] != null) {
						minimizer = (canonicalKmers[j].compareTo(minimizer) < 0) ? canonicalKmers[j] : minimizer;
					}
				}
				//System.out.println("Minimizer: " + minimizer);
				if(minimizer.length() == k) { //If the invalid letter N occurs in the exact center of a window, all kmers in the window will contain N and be invalid, Z will get through
					minimizers.add(minimizer);
				}
			}
			
			//left end minimizers, dynamic programming
			String minimizer = "Z";
			if(canonicalKmers[0] != null) {
				minimizers.add(canonicalKmers[0]);
				minimizer = canonicalKmers[0];
			}
			for(int v = 1; v <= (k - 2); v++) { //here v is not the u in (u,k) in the paper, it is the index of the last/rightmost kmer in each (u,k) window (i.e. u - 1)
				if(canonicalKmers[v] != null) {
					minimizer = (canonicalKmers[v].compareTo(minimizer) < 0) ? canonicalKmers[v] : minimizer;
					minimizers.add(minimizer);
				}
			}
			
			//right end minimizers, dynamic programming
			minimizer = "Z";
			if(canonicalKmers[startIndexOfLastKmer] != null) {
				minimizers.add(canonicalKmers[startIndexOfLastKmer]);
				minimizer = canonicalKmers[startIndexOfLastKmer];
			}
			for(int v = startIndexOfLastKmer - 1; v >= (startIndexOfLastKmer - k + 2); v--) { //here v is not the u in (u,k) in the paper, it is the index of the first/leftmost kmer in each (u,k) window
				if(canonicalKmers[v] != null) {
					minimizer = (canonicalKmers[v].compareTo(minimizer) < 0) ? canonicalKmers[v] : minimizer;
					minimizers.add(minimizer);
				}
			}
			
			for(String aMinimizer : minimizers) {
				kmers.add(convertKmerToIndex(aMinimizer));
			}
		}
		
		private int convertKmerToIndex(String kmer) {
			int index = 0;
			
			int lastCharIndex = k - 1;
			for(int i = lastCharIndex; i >= 0; i--) {
				switch(kmer.charAt(i)) {
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
			
			return index;
		}
		
		private MutableIntSet checkIfFragOnlyContainsValidCharacters(String frag) {
			MutableIntSet indicesOfInvalidKmers = new IntHashSet();
			
			int length = frag.length();
			for(int i = 0; i < length; i++) {
				char aChar = frag.charAt(i);
				switch(aChar) {
					case 'A':
					case 'C':
					case 'G':
					case 'T':
						continue;
					default:
						for(int j = i - k + 1; j <= i; j++) {
							indicesOfInvalidKmers.add(j);
						}
				}
			}
			
			return indicesOfInvalidKmers;
		}
		
		private String getMinusKmer(String plusKmer) {
			String minusKmer = "";
			
			for(int i = k - 1; i >= 0; i--) {
				switch(plusKmer.charAt(i)) {
					case 'A':
						minusKmer += "T";
						break;
					case 'T':
						minusKmer += "A";
						break;
					case 'G':
						minusKmer += "C";
						break;
					case 'C':
						minusKmer += "G";
						break;
				}
			}
			
			return minusKmer;
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
						System.out.println("Producer - start reading paired-end FASTA files " + startPath + " and " + endPath);
						readTestFragFastaFiles(reader1, reader2);
					} else {
						System.out.println("Producer - start reading paired-end FASTQ files " + startPath + " and " + endPath);
						readTestFragFastqFiles(reader1, reader2);
					}					
				} else {				
					if(startPath.contains(".fasta")) {
						System.out.println("Producer - start reading FASTA file " + startPath);
						readTestFragFastaFile(reader1);
					} else {
						System.out.println("Producer - start reading FASTQ file " + startPath);
						readTestFragFastqFile(reader1);
					}					
				}
				
				System.out.println("Producer - finished reading " + readCounter + " reads");
			} catch(Exception e) {
				System.out.println("Producer - error occurred, exiting");
				e.printStackTrace();
				System.exit(1);
			}			
		}
		
		private void readTestFragFastaFile(BufferedReader reader) throws Exception {
			String line = null;
			if(finishedReadIds == null) {
				while((line = reader.readLine()) != null) {
					line = line.trim();
					if(line.startsWith(">") || line.startsWith("@")) {
						String readId = line.substring(1).split("\\s+")[0];						
						readQueue.put(new String[] {readId, reader.readLine().trim().toUpperCase()});
						readCounter++;						
					}
				}
			} else {
				while((line = reader.readLine()) != null) {
					line = line.trim();
					if(line.startsWith(">") || line.startsWith("@")) {
						String readId = line.substring(1).split("\\s+")[0];
						if(finishedReadIds.contains(readId)) {
							reader.readLine();						
						} else {
							readQueue.put(new String[] {readId, reader.readLine().trim().toUpperCase()});
							readCounter++;
						}
					}
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
						
						readQueue.put(new String[] {readId, reader1.readLine().trim().toUpperCase() + "-" + reader2.readLine().trim().toUpperCase()});
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
							readQueue.put(new String[] {readId, reader1.readLine().trim().toUpperCase() + "-" + reader2.readLine().trim().toUpperCase()});
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

						readQueue.put(new String[] {readId, startSeq + "-" + endSeq});
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

							readQueue.put(new String[] {readId, startSeq + "-" + endSeq});
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
			
			try {
				BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(countFile)), "UTF-8"));
				String line = reader.readLine();
				logFres[id] = (float) Math.log(1.0 / Integer.parseInt(line));
				
				genomeMinimizers[id] = new IntHashSet();
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
		System.out.println("This NBCClassify program (v1.0) classify reads against a reference database.");
		System.out.println("-h:	Show this help menu");
		System.out.println("-k:	K-mer length");
		System.out.println("-c:	Number of threads");
		System.out.println("-d:	Input database directory");
		System.out.println("-m:	Input metainfo file");
		System.out.println("-o:	Final classification file");
		System.out.println("-t:	Type of reads (Paired-end: 2, Single-end: 1). Paired-end reads have two .fasta/.fastq files (can be gzipped) following; single-end reads have one file (can be gzipped).");
		System.out.println("-p (optional): Penalty for non-existent k-mers (default -2000)");
	}
}