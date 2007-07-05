package util;

import java.io.File;
import java.util.StringTokenizer;

public class BSAFileUtil {
	
	private static String prefix = "PeakVals_USC";
	private static String rupVarPrefix = "PeakVals_USC";
	
	public static int getSourceIDFromFile(File file) {
		return getIDFromTokens(file, 0);
	
	}
	
	public static int getRuptureIDFromFile(File file) {
		return getIDFromTokens(file,1);
	}
	
	private static int getIDFromTokens(File file, int indexToToken) {
		if (file.getName().startsWith(prefix)) {
			String filename = file.getName();
			String endName = filename.substring(prefix.length()+1);
			
			String[] tokens = endName.split("_|.bsa");
			for (int i=0; i<tokens.length; i++) {
				//System.out.println("token " + i + ": " + tokens[i]);
			}
			
			
			//System.out.println("endName: " + endName);
			
			return Integer.parseInt(tokens[indexToToken]);
		}
		else {
			return -1;	
		}
		
	}

	public static int getSourceIDFromRuptureVariationFile(File file) {
		return getIDFromRuptureVariationsFileTokens(file, 0);
	}
	
	public static int getRuptureIDFromRuptureVariationFile(File file) {
		return getIDFromRuptureVariationsFileTokens(file, 1);
	}

	public static int getRupVarIDFromRuptureVariationFile(File file) {
		return getIDFromRuptureVariationsFileTokens(file, 2);
	}
	
	private static int getIDFromRuptureVariationsFileTokens(File file, int indexToToken) {
		if (file.getName().startsWith(rupVarPrefix)) {
			String filename = file.getName();
			String endName = filename.substring(rupVarPrefix.length()+1);
			
			String[] tokens = endName.split("_|.bsa");
			for (int i=0; i<tokens.length; i++) {
				//System.out.println("token " + i + ": " + tokens[i]);
			}
			
			
			//System.out.println("endName: " + endName);
			
			return Integer.parseInt(tokens[indexToToken]);
		}
		else {
			return -1;	
		}
	}






}
