package util;

import java.io.File;

public class BSAFileUtil {
	
	private static String siteName;
	public static int getSourceIDFromFile(File file) {
		return getIDFromTokens(file, 0);
	
	}
	
	public static int getRuptureIDFromFile(File file) {
		return getIDFromTokens(file,1);
	}
	
	private static int getIDFromTokens(File file, int indexToToken) {

		String filename = file.getName();
		int lastIndex = filename.lastIndexOf(siteName);

		String endName = filename.substring(lastIndex+siteName.length()+1);
		// System.out.println("endName: " + endName);

		String[] tokens = endName.split("_|.bsa");
		/*for (int i=0; i<tokens.length; i++) {
			System.out.println("token " + i + ": " + tokens[i]);
		}*/

		return Integer.parseInt(tokens[indexToToken]);
		
	}

	public static int getSourceIDFromRuptureVariationFile(File file) {
		return getIDFromTokens(file, 0);
	}
	
	public static int getRuptureIDFromRuptureVariationFile(File file) {
		return getIDFromTokens(file, 1);
	}

	public static int getRupVarIDFromRuptureVariationFile(File file) {
		return getIDFromTokens(file, 2);
	}
	
	public static void setSiteName(String newSiteName) {
		siteName = newSiteName;
		
	}






}
