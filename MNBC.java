/**
 * 
 * @author Ruipeng Lu (ruipeng.lu@inspection.gc.ca)
 *
 */

public class MNBC {
	public static void main(String[] args) {
		if(args == null || args.length == 0) {
			help();
			System.exit(0);
		}
		
		if(args[0].equals("taxonomy")) {
			MNBC_taxonomy.execute(args);
		} else if(args[0].equals("build")) {
			MNBC_build.execute(args);
		} else if(args[0].equals("classify")) {
			MNBC_classify.execute(args);
		} else {
			help();
		}
	}
	
	private static void help() {
		System.out.println("MNBC (v1.1): a taxonomic read classifier");
		System.out.println("Step 1: generate taxonomy file -- Run 'MNBC taxonomy -h' for help");
		System.out.println("Step 2: build reference database -- Run 'MNBC build -h' for help");
		System.out.println("Step 3: classify reads -- Run 'MNBC classify -h' for help");
	}
}
