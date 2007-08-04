package util;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

public class BSAFileUtil {
	
	public static String pathLog;
	public static ArrayList<String> totalFilenameList;
	//public static File[] totalFileArray;
	public static ArrayList<File> totalFileList;
	
	public static boolean isInDebugMode = false;
	
	public static ArrayList<File> createTotalFileList(File saFile) {
		totalFilenameList = new ArrayList<String>();
		totalFileList = new ArrayList<File>();
		createTotalFileListHelper(saFile);
		return totalFileList;
		
	}
	
	private static void createTotalFileListHelper(File saFile) {
		File[] safilesList = saFile.listFiles(new BSAFilenameFilter());
		File[] sadirsList = saFile.listFiles(new NonCVSDirFileFilter());
		
		if (isInDebugMode) {
			System.out.println("Printing all files in " + saFile.getName());
			File[] totalSAFilesList = saFile.listFiles();
			for (int i=0; i<totalSAFilesList.length; i++) {
				System.out.println("totalSAFilesList[" + i + "]: " + totalSAFilesList[i].getName()); 
			}
		}
		
		if (isInDebugMode) {
			for (int i=0; i < safilesList.length; i++) {
				System.out.println("safilesList[" + i + "]: " + safilesList[i].getName());
			}			
		}
		
		if (isInDebugMode) {
			for (int i=0; i < sadirsList.length; i++) {
				System.out.println("sadirsList[" + i + "]: " + sadirsList[i].getName());
			}			
		}

		
		if (!saFile.getName().equals("CVS") && isInDebugMode) {
			System.out.println("Path: " + saFile.getPath());
		}
		
		if (pathLog == null && !saFile.getName().equals("CVS")) {
			pathLog = saFile.getPath() + " ";
		}
		else if (pathLog != null && !saFile.getName().equals("CVS")) {
			pathLog += saFile.getPath() + " ";
		}
		
		totalFileList.addAll(Arrays.asList(safilesList));
		for (int filesIndex=0; filesIndex < safilesList.length; filesIndex++) {
			if (isInDebugMode) {
				System.out.println("\tFilename: " + safilesList[filesIndex].getName());
			}
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
		if (isInDebugMode) {
			System.out.println("BSAFileUtil::getIDFromTokens: filename: " + filename);
		}
		int lastIndex = filename.lastIndexOf(siteName);
		
		if (isInDebugMode) {
			System.out.println("BSAFileUtil::getIDFromTokens: lastIndex: " + lastIndex);
		}

		String endName = filename.substring(lastIndex+siteName.length()+1);
		
		if (isInDebugMode) {
			System.out.println("endName: " + endName);
		}
		

		String[] tokens = endName.split("_|.bsa");

		if (isInDebugMode) {
			for (int i=0; i<tokens.length; i++) {
				System.out.println("token " + i + ": " + tokens[i]);
			}			
		}


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

	public static boolean isInDebugMode() {
		return isInDebugMode;
	}

	public static void setInDebugMode(boolean isInDebugMode) {
		BSAFileUtil.isInDebugMode = isInDebugMode;
	}






}
