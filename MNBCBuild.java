import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.eclipse.collections.api.set.primitive.MutableIntSet;
import org.eclipse.collections.impl.set.mutable.primitive.IntHashSet;

public class MNBCBuild { //Based on NaiveBayesClassifierCount_V3, only use canonical kmers; And only use slurm-jobid.out as the progress file
						//Only use minimizer seeds (w=k, window size=w+k-1), base/kmer ordering can change (here use default alphabetical ACGT order)
	private static int k;
	private static int numberOfThreads;
	private static String referenceGenomeDirPath;	
	private static String outputDirPath;
	private static String previousProgressPath;
	
	public static void main(String[] args) {
		for(int i = 0; i < args.length; i++) {
			if(args[i].startsWith("-")) {
				switch(args[i].charAt(1)) {
					case 'k':
						k = Integer.parseInt(args[i + 1]);
						break;
					case 'c':
						numberOfThreads = Integer.parseInt(args[i + 1]);
						break;
					case 'i':
						referenceGenomeDirPath = args[i + 1];
						break;
					case 'o':
						outputDirPath = args[i + 1];
						break;
					case 'b':
						previousProgressPath = args[i + 1];
						break;
					case 'h':
						printHelpInfo();
						System.exit(0);
				}
			}
		}
		
		long startTime = System.nanoTime();
		HashSet<String> previouslyCompletedGenomes = null;
		if(previousProgressPath != null) {
			System.out.println("***************************************************************");
			System.out.println("In the previous killed run, the following reference genomes have been completed:");			
			previouslyCompletedGenomes = readPreviousProgressFile();
			System.out.println("***************************************************************");
		}
		
		int numberOfCores = Runtime.getRuntime().availableProcessors();
		System.out.println("Number of available cores: " + numberOfCores);
		if(numberOfThreads > numberOfCores) {
			System.out.println("WARNING - Number of available cores " + numberOfCores + " is less than requested number of threads " + numberOfThreads);
			numberOfThreads = numberOfCores - 1;
			System.out.println("Will use " + numberOfThreads + " threads instead");
		}
		
		ExecutorService nested = Executors.newFixedThreadPool(numberOfThreads);
		CompletionService<String> pool = new ExecutorCompletionService<String>(nested);
		System.out.println("Created a thread pool of size " + numberOfThreads);
		
		int taskCounter = 0;
		File[] trainingGenomes = new File(referenceGenomeDirPath).listFiles();
		if(previousProgressPath == null) {
			for(int i = 0; i < trainingGenomes.length; i++) {
				pool.submit(new ReferenceGenomeProcessor(trainingGenomes[i], i));
			}
			taskCounter = trainingGenomes.length;			
		} else {
			for(int i = 0; i < trainingGenomes.length; i++) {
				if(!previouslyCompletedGenomes.contains(trainingGenomes[i].getName())) {
					pool.submit(new ReferenceGenomeProcessor(trainingGenomes[i], i));
					taskCounter++;
				}				
			}			
		}		
		System.out.println("Submitted all " + taskCounter + " tasks to thread pool");
		
		for(int i = 0; i < taskCounter; i++) {
			try {
				System.out.println("Trying to get outcome of " + i + "th returned task...");
				String outcome = pool.take().get();
				
				if(outcome.contains("Finished")) {
					System.out.println("Congratulations! This task is successful: " + outcome);
				} else {
					System.out.println("Unfortunately this task isn't successful: " + outcome);
				}
			} catch(ExecutionException e) {
				System.out.println("ExecutionException - thrown on " + i + " th returned task");
				e.printStackTrace();
			} catch(InterruptedException e) {
				System.out.println("InterruptedException - thrown on " + i + " th returned task");
				e.printStackTrace();
			}
		}
		
		nested.shutdown();
		long endTime = System.nanoTime();
		System.out.println("done in " + ((endTime - startTime) / 1000000000) + " seconds");		
	}
	
	private static void printHelpInfo() {
		System.out.println("This NBCBuild program builds a reference genome database from a set of genome sequence files by k-mer counting.");
		System.out.println("-h:	Show this help menu");
		System.out.println("-k:	K-mer length");
		System.out.println("-c:	Number of threads");
		System.out.println("-i:	Input directory containing the (filtered) sequence files of reference genomes downloaded from NCBI RefSeq");
		System.out.println("-o:	Output database directory");
		System.out.println("-b (optional): Log file of the previous abnormally killed run (.out file in Slurm)");
	}
	
	private static HashSet<String> readPreviousProgressFile() {
		HashSet<String> previouslyCompletedGenomes = new HashSet<String>();
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader(previousProgressPath));
			String line = null;
			while((line = reader.readLine()) != null) {
				if(line.contains("Finished")) {
					String[] fields = line.split("\\s+");
					previouslyCompletedGenomes.add(fields[13]);
					System.out.println(line);
				}
			}
			reader.close();
		} catch(Exception e) {
			System.out.println("Exception - can't read the previous Progress file " + previousProgressPath);
			e.printStackTrace();
			System.exit(1);
		}
		
		return previouslyCompletedGenomes;
	}
	
	private static class ReferenceGenomeProcessor implements Callable<String> {
		private File referenceGenome;
		private int id;
		
		public ReferenceGenomeProcessor(File aReferenceGenome, int anID) {
			referenceGenome = aReferenceGenome;
			id = anID;
		}

		@Override
		public String call() {			
			//long startTime = System.nanoTime();
			String filename = referenceGenome.getName();
			System.out.println("Task " + id + " - start processing genome " + filename + "...");
			HashSet<String> minimizers = new HashSet<String>();
			int kmerTotalCount = 0; //Total number of valid kmers in both strands
			
			try {
				ArrayList<String> chromosomes = readGenomeFile(referenceGenome);
				System.out.println("Task " + id + " - read " + chromosomes.size() + " chromosomes");
				
				for(int m = 0; m < chromosomes.size(); m++) {
					System.out.println("Task " + id + " - start processing " + m + "th chromosome...");
					String chromosome = chromosomes.get(m);
					
					MutableIntSet indicesOfInvalidKmers = checkIfFragOnlyContainsValidCharacters(chromosome);
					int indexOfLastKmerStartPosition = chromosome.length() - k;
					kmerTotalCount += indexOfLastKmerStartPosition + 1 - indicesOfInvalidKmers.size();
					
					String[] canonicalKmers = new String[indexOfLastKmerStartPosition + 1]; //first get the canonical kmer at each position
					if(indicesOfInvalidKmers.isEmpty()) {
						for(int i = 0; i <= indexOfLastKmerStartPosition; i++) {
							String plusKmer = chromosome.substring(i, i + k);
							String minusKmer = getMinusKmer(plusKmer);
							canonicalKmers[i] = (plusKmer.compareTo(minusKmer) < 0) ? plusKmer : minusKmer;
						}
					} else {
						for(int i = 0; i <= indexOfLastKmerStartPosition; i++) {
							if(!indicesOfInvalidKmers.contains(i)) {
								String plusKmer = chromosome.substring(i, i + k);
								String minusKmer = getMinusKmer(plusKmer);
								canonicalKmers[i] = (plusKmer.compareTo(minusKmer) < 0) ? plusKmer : minusKmer;
							}
						}
					}
					
					//interior minimizers, maybe could be optimized
					int indexOfLastWindowStartPosition = indexOfLastKmerStartPosition - k + 1;					
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
					if(canonicalKmers[indexOfLastKmerStartPosition] != null) {
						minimizers.add(canonicalKmers[indexOfLastKmerStartPosition]);
						minimizer = canonicalKmers[indexOfLastKmerStartPosition];
					}
					for(int v = indexOfLastKmerStartPosition - 1; v >= (indexOfLastKmerStartPosition - k + 2); v--) { //here v is not the u in (u,k) in the paper, it is the index of the first/leftmost kmer in each (u,k) window
						if(canonicalKmers[v] != null) {
							minimizer = (canonicalKmers[v].compareTo(minimizer) < 0) ? canonicalKmers[v] : minimizer;
							minimizers.add(minimizer);
						}
					}
				}
				
				kmerTotalCount = kmerTotalCount * 2;
				System.out.println("Task " + id + " - valid k-mer count: " + kmerTotalCount + ", minimizer count: " + minimizers.size());				
			} catch(Exception e) {
				e.printStackTrace();
				return "Task " + id + " - Exception - can't find file " + filename;
			}
			
			try {
				PrintWriter writer = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDirPath + "/" + filename.substring(0, filename.length() - (filename.endsWith(".gz") ? 7 : 4)) + "_count.gz")), "UTF-8")), true);
				writer.println(kmerTotalCount);
				for(String minimizer : minimizers) {
					//System.out.println("Minimizer " + minimizer);
					writer.println(convertKmerToIndex(minimizer));					
				}
				writer.close();
			} catch(IOException e) {
				e.printStackTrace();
				return "ExceptionIO - can't open or write count file for reference genome: " + filename;
			}
			//long endTime = System.nanoTime();
			//long runningTime = (endTime - startTime) / 1000000000;
			
			return "Task " + id + " - Finished the genome count file " + filename/* + " in " + runningTime + " seconds"*/;
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
		
		private ArrayList<String> readGenomeFile(File genomeFile) throws FileNotFoundException, IOException {
			ArrayList<String> chromosomes = new ArrayList<String>();		
			String chromosome = "";
			
			BufferedReader reader = null;
			if(genomeFile.getName().endsWith(".gz")) {
				reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(genomeFile)), "UTF-8"));
			} else {
				reader = new BufferedReader(new FileReader(genomeFile));
			}
								
			String line = null;
			while((line = reader.readLine()) != null) {
				if(line.startsWith(">")) {
					if(!chromosome.isEmpty()) {						
						chromosomes.add(chromosome.toUpperCase());						
						chromosome = "";
					}
				} else {					
					chromosome = chromosome.concat(line);
				}
			}			
			if(!chromosome.isEmpty()) {
				chromosomes.add(chromosome.toUpperCase());
			}				
			reader.close();			
			
			return chromosomes;
		}
	}
}
