package data;

import java.io.File;
import java.util.ArrayList;
import java.util.zip.ZipEntry;

import util.BSAFileUtil;
import util.NumberHelper;
import util.SABinary2Float;
import util.SwapBytes;


public class SARuptureFromRuptureVariationFile extends SARuptureFromFile {

	public SARuptureVariation rupVar = new SARuptureVariation();

	public SARuptureFromRuptureVariationFile(File bsaFile, String siteName) {
		byte[] byteArrayFromFile = SwapBytes.swapBinaryFileToArrayForFloats(bsaFile.getPath());
		if (byteArrayFromFile == null) {
			throw new NullPointerException("byteArrayFromFile is null; bsaFile.getPath(): " + bsaFile.getPath());
		}
		ArrayList<Float> floats = SABinary2Float.convertByteArrayToArrayOfFloats(byteArrayFromFile);
		super.createRupVars(floats);
		sourceID = BSAFileUtil.getSourceIDFromRuptureVariationFile(bsaFile, siteName);
		ruptureID = BSAFileUtil.getRuptureIDFromRuptureVariationFile(bsaFile, siteName);
		rupVar.variationNumber = BSAFileUtil.getRupVarIDFromRuptureVariationFile(bsaFile, siteName);
		createPeriods(floats);
		try {
			rupVar.computeGeomAvg();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public SARuptureFromRuptureVariationFile(byte[] data, String siteName, ZipEntry saZipEntry) {
		ArrayList<Float> floats = SABinary2Float.convertByteArrayToArrayOfFloats(data);
		ArrayList<Float> swappedFloats = SwapBytes.swapArrayOfFloatsToArrayOfFloats(floats);
		super.createRupVars(swappedFloats);
		sourceID = BSAFileUtil.getSourceIDFromRuptureVariationZipEntry(saZipEntry, siteName);
		ruptureID = BSAFileUtil.getRuptureIDFromRuptureVariationZipEntry(saZipEntry, siteName);
		rupVar.variationNumber = BSAFileUtil.getRupVarIDFromRuptureVariationZipEntry(saZipEntry, siteName);
		createPeriods(swappedFloats);
		try {
			rupVar.computeGeomAvg();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public SARuptureFromRuptureVariationFile() {
	}

	public void createPeriods(ArrayList<Float> floats) {
		int maxPeriod = 27;
		int periodCount = 1;
		String componentDirection = "East";
		
		for (int i=0; i<floats.size(); i++) {
			if (periodCount > maxPeriod) {
				periodCount = 1;
				if (componentDirection.equals("East")) {
					// Swap component direction
					componentDirection = "North";
				}
				else if (componentDirection.equals("North")){
					// Swap component direction
					componentDirection = "East";
				}
			}
			
			//System.out.println("SA for Rupture Variation " + rupVar.variationNumber + ", " + componentDirection + " Component, Period " + periodCount + " : " + floats.get(i) );

			if (componentDirection.equals("North")) {
				rupVar.northComp.periods[periodCount-1] = floats.get(i);
			}
			else if (componentDirection.equals("East")) {
				rupVar.eastComp.periods[periodCount-1] = floats.get(i);
			}
			periodCount++;
		}
	}
}

