package util;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

public class BSAFileUtil {
	
	public static String pathLog;
	public static ArrayList<String> totalFilenameList;
	//public static File[] totalFileArray;
	public static ArrayList<File> totalFileList;
	
	public static ArrayList<File> createTotalFileList(File saFile) {
		totalFilenameList = new ArrayList<String>();
		totalFileList = new ArrayList<File>();
		createTotalFileListHelper(saFile);
		return totalFileList;
		
	}
	
	private static void createTotalFileListHelper(File saFile) {
		File[] safilesList = saFile.listFiles(new BSAFilenameFilter());
		File[] sadirsList = saFile.listFiles(new NonCVSDirFileFilter());
		if (!saFile.getName().equals("CVS")) {
			//System.out.println("Path: " + saFile.getPath());
		}
		
		if (pathLog == null && !saFile.getName().equals("CVS")) {
			pathLog = saFile.getPath() + " ";
		}
		else if (!saFile.getName().equals("CVS")) {
			pathLog += saFile.getPath() + " ";
		}
/*		if (totalFileList == null) {
			totalFileList = new ArrayList<File>();
		}
		else {
			System.out.println("BSAFileUtil.totalFileList is NOT null");
		}*/
		totalFileList.addAll(Arrays.asList(safilesList));
		for (int filesIndex=0; filesIndex < safilesList.length; filesIndex++) {
			/*System.out.println("\tFilename: " + safilesList[filesIndex].getName());*/
			totalFilenameList.add(safilesList[filesIndex].getName());
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
	
	public static int getSourceIDFromFile(File file, String siteName) {
		return getIDFromTokens(file, 0, siteName);
	
	}
	
	public static int getRuptureIDFromFile(File file, String siteName) {
		return getIDFromTokens(file, 1, siteName);
	}
	
	private static int getIDFromTokens(File file, int indexToToken, String siteName) {

		String filename = file.getName();
		//System.out.println("BSAFileUtil::getIDFromTokens: filename: " + filename);
		int lastIndex = filename.lastIndexOf(siteName);
		//System.out.println("BSAFileUtil::getIDFromTokens: lastIndex: " + lastIndex);

		String endName = filename.substring(lastIndex+siteName.length()+1);
		// System.out.println("endName: " + endName);

		String[] tokens = endName.split("_|.bsa");
		/*for (int i=0; i<tokens.length; i++) {
			System.out.println("token " + i + ": " + tokens[i]);
		}*/

		return Integer.parseInt(tokens[indexToToken]);
		
	}

	public static int getSourceIDFromRuptureVariationFile(File file, String siteName) {
		return getIDFromTokens(file, 0, siteName);
	}
	
	public static int getRuptureIDFromRuptureVariationFile(File file, String siteName) {
		return getIDFromTokens(file, 1, siteName);
	}

	public static int getRupVarIDFromRuptureVariationFile(File file, String siteName) {
		return getIDFromTokens(file, 2, siteName);
	}






}
