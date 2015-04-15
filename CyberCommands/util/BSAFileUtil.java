package util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import commands.CyberLoadamps.Mode;

public class BSAFileUtil {
	
	private static final int RUP_VAR_ID_POS = 2;
	private static final int RUPTURE_ID_POS = 1;
	private static final int SOURCE_ID_POS = 0;
	public static String pathLog;
	public static ArrayList<String> totalFilenameList;
	public static ArrayList<File> totalFileList;
	
	public static boolean isInDebugMode = false;
	private static boolean hasSGTVariationCharacter = false;
	
	public static ArrayList<File> createTotalFileList(File saFile, Mode fileMode) {
		totalFilenameList = new ArrayList<String>();
		totalFileList = new ArrayList<File>();
		if (fileMode==Mode.ZIP) {
			createTotalFileListZipHelper(saFile);
		} else if (fileMode==Mode.BSA || fileMode==Mode.HEAD){
			//BSA or HEAD mode
			createTotalFileListHelper(saFile);
		} else if (fileMode==Mode.ROTD) {
			createTotalFileListRotd(saFile);
		}
		return totalFileList;
		
	}
	
	private static void createTotalFileListRotd(File saFile) {
		File[] rotdList;
		if (saFile.isDirectory()) {
			rotdList = saFile.listFiles(new RotDFilenameFilter());
		} else {
			//just inserting 1 file		totalFileList.addAll(Arrays.asList(zipfilesList));
			rotdList = new File[1];
			rotdList[0] = saFile;
		}
		
		if (isInDebugMode) {
			System.out.println("Printing all files in " + saFile.getName());
			for (int i=0; i<rotdList.length; i++) {
				System.out.println(rotdList[i]);
			}
		}
		
		totalFileList.addAll(Arrays.asList(rotdList));
		for (int filesIndex=0; filesIndex < rotdList.length; filesIndex++) {
			if (isInDebugMode) {
				System.out.println("\tFilename: " + rotdList[filesIndex].getName());
			}
			totalFilenameList.add(rotdList[filesIndex].getName());
		}
	}
	
	private static void createTotalFileListZipHelper(File saFile) {
		File[] zipfilesList;
		File[] sadirsList;
		if (saFile.isDirectory()) {
			zipfilesList = saFile.listFiles(new ZipFilenameFilter());
			sadirsList = saFile.listFiles(new NonCVSDirFileFilter());
		} else {
			//Just inserting 1 zip file
			zipfilesList = new File[1];
			zipfilesList[0] = saFile;
			sadirsList = null;
		}
		
		if (isInDebugMode) {
			System.out.println("Printing all files in " + saFile.getName());
			File[] totalSAFilesList = saFile.listFiles();
			for (int i=0; i<totalSAFilesList.length; i++) {
				System.out.println("totalSAFilesList[" + i + "]: " + totalSAFilesList[i].getName()); 
			}
		}
		
		if (isInDebugMode) {
			for (int i=0; i < zipfilesList.length; i++) {
				System.out.println("zipfilesList[" + i + "]: " + zipfilesList[i].getName());
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
		
		totalFileList.addAll(Arrays.asList(zipfilesList));
		for (int filesIndex=0; filesIndex < zipfilesList.length; filesIndex++) {
			if (isInDebugMode) {
				System.out.println("\tFilename: " + zipfilesList[filesIndex].getName());
			}
			totalFilenameList.add(zipfilesList[filesIndex].getName());
		}
		
		if (sadirsList != null) {
			for (int dirsIndex=0; dirsIndex<sadirsList.length; dirsIndex++) {
				createTotalFileListZipHelper(sadirsList[dirsIndex]);
			}
		}
		else {
			return;
		}
		
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
		return getIDFromTokens(file.getName(), SOURCE_ID_POS, siteName);
	
	}
	
	public static int getRuptureIDFromFile(File file, String siteName) {
		return getIDFromTokens(file.getName(), RUPTURE_ID_POS, siteName);
	}
	
	private static int getIDFromTokens(String filename, int indexToToken, String siteName) {
		if (isInDebugMode) {
			System.out.println("BSAFileUtil::getIDFromTokens: filename: " + filename);
		}
		int startingIndexOfSiteName = filename.lastIndexOf(siteName);
		int firstUnderscoreAfterLastIndex = filename.indexOf("_",startingIndexOfSiteName);
		int secondUnderscoreAfterLastIndex = filename.indexOf("_", firstUnderscoreAfterLastIndex+1);
		
		if (isInDebugMode) {
			System.out.println("BSAFileUtil::getIDFromTokens: firstUnderscoreAfterLastIndex: " + firstUnderscoreAfterLastIndex);
			System.out.println("BSAFileUtil::getIDFromTokens: secondUnderscoreAfterLastIndex: " + secondUnderscoreAfterLastIndex);
		}
		
		int sgtVariationCharLength = 0;
		
		if (firstUnderscoreAfterLastIndex != -1) {
			if (isInDebugMode) {
				System.out.println("Checking for non-digit characters representing non-one SGT Variation");
			}
			
			for (int i=firstUnderscoreAfterLastIndex+1; i<secondUnderscoreAfterLastIndex; i++) {
				if (isInDebugMode) {
					System.out.println(i + ": " + filename.charAt(i));
				}
				
				sgtVariationCharLength++;
				
				if (isInDebugMode) {
					System.out.println("Checking if " + filename.charAt(i) + " is a digit");
				}
				if (!Character.isDigit(filename.charAt(i))) {
					hasSGTVariationCharacter  = true;
					break;
				}
			}
		}
		
		if (isInDebugMode) {
			System.out.println("sgtVariationCharLength: " + sgtVariationCharLength);	
		}
		
		int lastIndexOfSGTVariationChar = firstUnderscoreAfterLastIndex+sgtVariationCharLength;
		
		
		if (isInDebugMode) {
			System.out.println("BSAFileUtil::getIDFromTokens: startingIndexOfSiteName: " + startingIndexOfSiteName);
			System.out.println("BSAFileUtil::getIDFromTokens: lastIndexOfSGTVariationChar: " + lastIndexOfSGTVariationChar);
		}
		
		String endName = new String();

		if (hasSGTVariationCharacter) {
			if (isInDebugMode) {
				System.out.println("BSAFileUtil::getIDFromTokens: " + filename + " has sgt variation characters");
			}
			endName = filename.substring(lastIndexOfSGTVariationChar+2);
		}
		else {
			endName = filename.substring(startingIndexOfSiteName+siteName.length()+1);
		}
		
		
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
		return getIDFromTokens(file.getName(), SOURCE_ID_POS, siteName);
	}
	
	public static int getRuptureIDFromRuptureVariationFile(File file, String siteName) {
		return getIDFromTokens(file.getName(), RUPTURE_ID_POS, siteName);
	}

	public static int getRupVarIDFromRuptureVariationFile(File file, String siteName) {
		return getIDFromTokens(file.getName(), RUP_VAR_ID_POS, siteName);
	}

	public static int getSourceIDFromRuptureVariationZipEntry(ZipEntry entry, String siteName) {
		return getIDFromTokens(entry.getName(), SOURCE_ID_POS, siteName);
	}
	
	public static int getRuptureIDFromRuptureVariationZipEntry(ZipEntry entry, String siteName) {
		return getIDFromTokens(entry.getName(), RUPTURE_ID_POS, siteName);
	}

	public static int getRupVarIDFromRuptureVariationZipEntry(ZipEntry entry, String siteName) {
		return getIDFromTokens(entry.getName(), RUP_VAR_ID_POS, siteName);
	}

	public static boolean isInDebugMode() {
		return isInDebugMode;
	}

	public static void setInDebugMode(boolean isInDebugMode) {
		BSAFileUtil.isInDebugMode = isInDebugMode;
	}






}
