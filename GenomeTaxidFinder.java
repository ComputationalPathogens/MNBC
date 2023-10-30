import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class GenomeTaxidFinder {
	public static void main(String[] args) {
		String genomeDirPath = args[0]; //Directory containing downloaded RefSeq reference genomes of interest, these genomes will be used as the reference database
		String refseqAssemblySummaryPath = args[1]; //Example filename: assembly_summary_refseq.txt
		String nodeDmpPath = args[2]; //Example filename: nodes.dmp
		String outputPath = args[3]; //Example filename: refSeq_prokaryote_complete_genomes_ok_status_metainfo.txt
		
		ArrayList<String> refseqAssemblyIDs = readGenomesSummaryTable(genomeDirPath);
		HashMap<String, String[]> refseqAssemblyID2Taxid = readRefseqAssemblySummary(refseqAssemblySummaryPath);
		HashMap<String, String> taxid2TaxLevel = new HashMap<String, String>();
		HashMap<String, String> taxid2ParentTaxid = readNodeDmp(nodeDmpPath, taxid2TaxLevel);
		
		refseqAssemblyID2Taxid.put("GCF_009387985.1", new String[] {"1", "2601652", "Bacillus phage 000TH010"}); //Patch for viruses
		refseqAssemblyID2Taxid.put("GCF_009388005.1", new String[] {"1", "2601660", "Bacillus phage 049ML001"}); //Patch for viruses
		refseqAssemblyID2Taxid.put("GCF_003613715.1", new String[] {"1", "2315627", "Bacillus phage Ray17"}); //Patch for viruses
		refseqAssemblyID2Taxid.put("GCF_008221585.1", new String[] {"1", "2591379", "Bacillus phage vB_BspS_SplendidRed"}); //Patch for viruses
		
		try {
			PrintWriter writer = new PrintWriter(outputPath);
			writer.print("RefSeq assembly ID\ttaxid.sub-species\ttaxid.species\ttaxid.genus\ttaxid.family\ttaxid.order\ttaxid.class\ttaxid.phylum\ttaxid.superkingdom\tOrganism name\n");
			generateGenomeTaxonomyTableRows(writer, refseqAssemblyIDs, refseqAssemblyID2Taxid, taxid2TaxLevel, taxid2ParentTaxid);
			writer.close();
		} catch(Exception e) {
			System.out.println("Error in writing results to " + outputPath);
			e.printStackTrace();
			System.exit(3);
		}
		
		System.out.println("done");
	}
	
	private static void generateGenomeTaxonomyTableRows(PrintWriter writer, ArrayList<String> refseqAssemblyIDs, HashMap<String, String[]> refseqAssemblyID2Taxid, HashMap<String, String> taxid2TaxLevel, HashMap<String, String> taxid2ParentTaxid) {
		for(String assemblyID : refseqAssemblyIDs) {
			String row = assemblyID;
			String[] ranks = new String[8];
			//System.out.println("assemblyID " + assemblyID);
			String[] taxidAndName = refseqAssemblyID2Taxid.get(assemblyID);
			
			if(taxidAndName[0].equals("1")) {
				ranks[0] = "";
				ranks[1] = taxidAndName[1];
				String supposedSpeciesLevel = taxid2TaxLevel.get(ranks[1]);
				if(!supposedSpeciesLevel.equals("species")) {
					System.out.println("Error: " + ranks[1] + " is supposed to be a species according to refseqAssemblySummary, but node.dmp reports " + supposedSpeciesLevel);
					System.exit(2);
				}
				
				fillRanksArray(ranks[1], ranks, taxid2TaxLevel, taxid2ParentTaxid);				
			} else {
				ranks[0] = taxidAndName[1];				
				/*String supposedtrainOrSubstrainLevel = taxid2TaxLevel.get(taxidAndName[1]);				
				if(supposedtrainOrSubstrainLevel.equals("strain")) {
					ranks[0] = taxidAndName[1];					
				} else {
					String parentSupposedStrainTaxid = taxid2ParentTaxid.get(taxidAndName[1]);
					String supposedStrainLevel = taxid2TaxLevel.get(parentSupposedStrainTaxid);
					if(supposedStrainLevel.equals("strain")) {
						ranks[0] = parentSupposedStrainTaxid;
					} else {
						ranks[0] = taxidAndName[1];						
						warning.println("Warning: " + taxidAndName[1] + " subspecies id is used as a strain id for " + taxidAndName[2] + "; parent species id " + parentSupposedStrainTaxid + " at level " + supposedStrainLevel);						
					}
				}*/
				
				fillRanksArray(ranks[0], ranks, taxid2TaxLevel, taxid2ParentTaxid);
			}
			
			for(String rank : ranks) {
				row += "\t" + rank;
			}
			row += "\t" + taxidAndName[2] + "\n";
			
			//System.out.println("One row result: " + row);
			writer.print(row);
		}
	}
	
	private static void fillRanksArray(String initial, String[] ranks, HashMap<String, String> taxid2TaxLevel, HashMap<String, String> taxid2ParentTaxid) {
		String currentTaxid = initial;
		while(!currentTaxid.isEmpty()) {
			String parentTaxid = taxid2ParentTaxid.get(currentTaxid);
			String level = taxid2TaxLevel.get(parentTaxid);
			
			if(level.equals("species")) {
				ranks[1] = parentTaxid;
				currentTaxid = parentTaxid;
				continue;
			} else if(level.equals("genus")) {
				ranks[2] = parentTaxid;
				currentTaxid = parentTaxid;
				continue;
			} else if(level.equals("family")) {
				ranks[3] = parentTaxid;
				currentTaxid = parentTaxid;
				continue;
			} else if(level.equals("order")) {
				ranks[4] = parentTaxid;
				currentTaxid = parentTaxid;
				continue;
			} else if(level.equals("class")) {
				ranks[5] = parentTaxid;
				currentTaxid = parentTaxid;
				continue;
			} else if(level.equals("phylum")) {
				ranks[6] = parentTaxid;
				currentTaxid = parentTaxid;
				continue;
			} else if(level.equals("superkingdom")) {
				ranks[7] = parentTaxid;
				currentTaxid = "";
			} else {
				currentTaxid = parentTaxid;
				continue;
			}										
		}
	}
	
	private static HashMap<String, String> readNodeDmp(String nodeDmpPath, HashMap<String, String> taxid2TaxLevel) {
		HashMap<String, String> taxid2ParentTaxid = new HashMap<String, String>();
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader(nodeDmpPath));
			String line = null;
			while((line = reader.readLine()) != null) {
				String[] fields = line.split("\\|");
				String taxid = fields[0].trim();
				taxid2ParentTaxid.put(taxid, fields[1].trim());
				taxid2TaxLevel.put(taxid, fields[2].trim());
			}
			reader.close();
		} catch(Exception e) {
			System.out.println("Error in reading node.dmp file: " + nodeDmpPath);
			e.printStackTrace();
			System.exit(1);
		}
		
		System.out.println("Finished reading node.dmp file : " + nodeDmpPath);
		return taxid2ParentTaxid;
	}
	
	private static HashMap<String, String[]> readRefseqAssemblySummary(String refseqAssemblySummaryPath) {
		HashMap<String, String[]> refseqAssemblyID2Taxid = new HashMap<String, String[]>(); //String[0]: indicating if taxid is at species level (1) or strain level (0), String[1]: species or strain taxid, String[2]: organism name
		
		try {			
			BufferedReader reader = new BufferedReader(new FileReader(refseqAssemblySummaryPath));
			String line = reader.readLine();
			line = reader.readLine();
			while((line = reader.readLine()) != null) {
				String[] fields = line.split("\t");
				if(fields[5].equals(fields[6])) {
					refseqAssemblyID2Taxid.put(fields[0], new String[] {"1", fields[5], fields[7] + " " + fields[8]});
				} else {
					refseqAssemblyID2Taxid.put(fields[0], new String[] {"0", fields[5], fields[7] + " " + fields[8]});
				}				
			}
			reader.close();
		} catch(Exception e) {
			System.out.println("Error in reading Refseq assembly summary file: " + refseqAssemblySummaryPath);
			e.printStackTrace();
			System.exit(1);
		}
		
		System.out.println("Finished reading Refseq assembly summary file : " + refseqAssemblySummaryPath);
		return refseqAssemblyID2Taxid;
	}
	
	private static ArrayList<String> readGenomesSummaryTable(String genomeDirPath) {
		ArrayList<String> refseqAssemblyIDs = new ArrayList<String>();
		
		File[] genomes = new File(genomeDirPath).listFiles();
		for(File genome : genomes) {
			String[] fields = genome.getName().split("_");
			refseqAssemblyIDs.add(fields[0] + "_" + fields[1]);
		}
		
		System.out.println("Finished reading genome IDs from " + genomeDirPath);
		return refseqAssemblyIDs;
	}
}
