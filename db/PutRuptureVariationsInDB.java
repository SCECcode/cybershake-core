


public class PutRuptureVariationsInDB {
	private static String filename;
	
	public static void main(String[] args) {
		if (args.length<1) {
			System.out.println("Usage: PutRuptureVariationsInDB < ruptVarFile | tree >");
			System.exit(0);
		}
		
		filename = args[0];
		
	}
}
