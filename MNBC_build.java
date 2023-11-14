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

public class MNBC_build { //Based on NaiveBayesClassifierCount_V3, only use canonical kmers; And only use slurm-jobid.out as the progress file
						//Only use minimizer seeds (w=k, window size=w+k-1), base/kmer ordering can change (here use default alphabetical ACGT order)
	private static int k;
	private static int numberOfThreads;
	private static int lengthThreshold = 0;
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
					case 'f':
						lengthThreshold = Integer.parseInt(args[i + 1]);
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
		int numberOfCores = Runtime.getRuntime().availableProcessors();
		System.out.println("Number of available cores: " + numberOfCores);
		if(numberOfThreads > numberOfCores) {
			System.out.println("WARNING - Number of available cores " + numberOfCores + " is less than requested number of threads " + numberOfThreads + ", exiting");
			System.exit(1);
		}
		
		HashSet<String> previouslyCompletedGenomes = null;
		if(previousProgressPath != null) {
			System.out.println("***************************************************************");
			System.out.println("In the previous killed run, the following reference genomes have been completed:");			
			previouslyCompletedGenomes = readPreviousProgressFile();
			System.out.println("***************************************************************");
		}
		
		ExecutorService nested = Executors.newFixedThreadPool(numberOfThreads);
		CompletionService<String> pool = new ExecutorCompletionService<String>(nested);
		System.out.println("Created a thread pool");
		
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
		System.out.println("Submitted all " + taskCounter + " tasks");
		
		for(int i = 0; i < taskCounter; i++) {
			try {
				System.out.println("Waiting to get outcome of " + i + "th returned task...");
				String outcome = pool.take().get();
				
				if(outcome.contains("Finished")) {
					System.out.println("Congratulations! This task is successful: " + outcome);
				} else {
					System.out.println("Unfortunately this task failed: " + outcome);
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
		System.out.println("This MNBC_build tool (v1.0) builds a reference database from a set of sequence files.");
		System.out.println("-h:	Show this help menu");
		System.out.println("-k:	K-mer length");
		System.out.println("-c:	Number of threads");		
		System.out.println("-i:	Input directory containing the (gzipped) files of reference sequences (e.g. *.fna.gz downloaded from RefSeq");
		System.out.println("-o:	Output database directory");
		System.out.println("-f (optional): The minimum filtering threshold on the sequence length (an integer not smaller than 0, the default value is 0). Sequences with lengths below this threshold will be ignored, and all plasmids will be ignored.");
		System.out.println("-b (optional): Log file of the previous prematurely killed run (i.e. Slurm .out file)");
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
					
					MutableIntSet indicesOfInvalidKmers = new IntHashSet();
					String minusSequence = "";
					int length = chromosome.length();
					for(int i = length - 1; i >= 0; i--) {
						switch(chromosome.charAt(i)) {
							case 'A':
								minusSequence += "T";
								break;
							case 'C':
								minusSequence += "G";
								break;
							case 'G':
								minusSequence += "C";
								break;
							case 'T':
								minusSequence += "A";
								break;
							default:
								for(int j = i - k + 1; j <= i; j++) {
									indicesOfInvalidKmers.add(j);
								}
						}
					}					
					
					int startIndexOfLastKmer = length - k;
					String[] canonicalKmers = new String[startIndexOfLastKmer + 1]; //first get the canonical kmer at each position
					kmerTotalCount += startIndexOfLastKmer + 1 - indicesOfInvalidKmers.size();
					
					if(indicesOfInvalidKmers.isEmpty()) {						
						for(int i = 0; i <= startIndexOfLastKmer; i++) {
							String plusKmer = chromosome.substring(i, i + k);
							String minusKmer = minusSequence.substring(startIndexOfLastKmer - i, length - i);
							canonicalKmers[i] = (plusKmer.compareTo(minusKmer) < 0) ? plusKmer : minusKmer;
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
					} else {
						for(int i = 0; i <= startIndexOfLastKmer; i++) {
							if(!indicesOfInvalidKmers.contains(i)) {
								String plusKmer = chromosome.substring(i, i + k);
								String minusKmer = minusSequence.substring(startIndexOfLastKmer - i, length - i);
								canonicalKmers[i] = (plusKmer.compareTo(minusKmer) < 0) ? plusKmer : minusKmer;
							}
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
				return "ExceptionIO - can't write count file for reference genome: " + filename;
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
		
		private ArrayList<String> readGenomeFile(File genomeFile) throws FileNotFoundException, IOException {
			ArrayList<String> chromosomes = new ArrayList<String>();		
			String chromosome = "";
			int chromosomeLength = 0;
			boolean retain = true;
			
			BufferedReader reader = null;
			if(genomeFile.getName().endsWith(".gz")) {
				reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(genomeFile)), "UTF-8"));
			} else {
				reader = new BufferedReader(new FileReader(genomeFile));
			}
								
			String line = null;
			while((line = reader.readLine()) != null) {
				if(line.startsWith(">")) {
					line = line.toLowerCase();
					
					if(chromosomeLength != 0) {						
						if(chromosomeLength >= lengthThreshold) {
							chromosomes.add(chromosome.toUpperCase());							
						}
						chromosomeLength = 0;
						chromosome = "";						
					}
					
					if(line.contains("plasmid")) {
						retain = false;
						continue;
					} else {
						retain = true;
						continue;
					}
				} else {
					if(retain) {
						chromosomeLength += line.length();
						chromosome = chromosome.concat(line);
						continue;
					}
				}
			}			
			if(chromosomeLength != 0) {
				if(chromosomeLength >= lengthThreshold) {
					chromosomes.add(chromosome.toUpperCase());
				}
			}			
			reader.close();			
			
			return chromosomes;
		}
	}
}
