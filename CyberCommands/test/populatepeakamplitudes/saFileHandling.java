package test.populatepeakamplitudes;

import java.io.File;

import util.BSAFileUtil;
import util.SwapBytes;

public class saFileHandling {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		File bsaFile = new File("safiles/PeakVals_allUSC_127_6.bsa");
		System.out.println("bsaFile path and name: " + bsaFile.getPath() + bsaFile.getName());
		
		//bsaFileUtilTest();
		//printAllFilesinSAFiles();
		bsaFileUtilForRuptureVariationsFileTest();
		

	}

	private static void bsaFileUtilForRuptureVariationsFileTest() {
		File bsaFile = new File("safiles/rupturevariations/PeakVals_USC_127_6_1.bsa");
		int source_ID = BSAFileUtil.getSourceIDFromRuptureVariationFile(bsaFile);
		int rupture_ID = BSAFileUtil.getRuptureIDFromRuptureVariationFile(bsaFile);
		int rupVar_ID = BSAFileUtil.getRupVarIDFromRuptureVariationFile(bsaFile);
		System.out.println("Source_ID: " + source_ID);
		System.out.println("Rupture_ID: " + rupture_ID);
		System.out.println("RupVar_ID: " + rupVar_ID);
		
	}

	private static void bsaFileUtilTest() {
		File bsaFile = new File("safiles/PeakVals_allUSC_127_6.bsa");
		int source_ID = BSAFileUtil.getSourceIDFromFile(bsaFile);
		int rupture_ID = BSAFileUtil.getRuptureIDFromFile(bsaFile);
		System.out.println("Source_ID: " + source_ID);
		System.out.println("Rupture_ID: " + rupture_ID);
	}

	private static void printAllFilesinSAFiles() {
		File safiles = new File("safiles");
		File[] safilesList =  safiles.listFiles();
		for (int i=0; i<safilesList.length; i++) {
			System.out.println("Filename: " + safilesList[i].getName());
		}
	}

}
