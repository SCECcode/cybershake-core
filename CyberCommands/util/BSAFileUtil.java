package util;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

public class BSAFileUtil {
	
	public static String siteName;
	public static String pathLog;
	public static ArrayList<String> totalFilenameList;
	//public static File[] totalFileArray;
	public static ArrayList<File> totalFileList;
	
	public static ArrayList<File> createTotalFileList(File saFile) {
		createTotalFileListHelper(saFile);
		return totalFileList;
		
	}
	
	public static void createTotalFileListHelper(File saFile) {
		File[] safilesList = saFile.listFiles(new BSAFilenameFilter());
		File[] sadirsList = saFile.listFiles(new DirFileFilter());
		
		System.out.println("Path: " + saFile.getPath() + saFile.getName());
		if (pathLog == null) {
			pathLog = saFile.getPath() + saFile.getName() + " ";
		}
		else {
			pathLog += saFile.getPath() + saFile.getName() + " ";
		}
		
		totalFileList.addAll(Arrays.asList(safilesList));
		for (int filesIndex=0; filesIndex < safilesList.length; filesIndex++) {
			System.out.println("\tFilename: " + safilesList[filesIndex].getName());
			totalFilenameList.add(safilesList[filesIndex].getName());
/*			File[] fileArrayDest = new File[safilesList.length + totalFileArray.length];
			System.arraycopy(safilesList, 0, fileArrayDest, 0, safilesList.length);
			System.arraycopy(totalFileArray, 0, fileArrayDest, safilesList.length, totalFileArray.length);
			totalFileArray = fileArrayDest;*/
		}
		
		if (sadirsList != null) {
			for (int dirsIndex=0; dirsIndex<sadirsList.length; dirsIndex++) {
				createTotalFileListHelper(sadirsList[dirsIndex]);
			}
		}
		else {
			return;
		}
		
	}
	
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
