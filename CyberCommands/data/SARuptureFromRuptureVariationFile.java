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
		super.createRupVars(floats, false);
		sourceID = BSAFileUtil.getSourceIDFromRuptureVariationFile(bsaFile, siteName);
		ruptureID = BSAFileUtil.getRuptureIDFromRuptureVariationFile(bsaFile, siteName);
		rupVar.variationNumber = BSAFileUtil.getRupVarIDFromRuptureVariationFile(bsaFile, siteName);
		createPeriods(floats, false);
		try {
			rupVar.computeGeomAvg();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public SARuptureFromRuptureVariationFile(byte[] data, String siteName, ZipEntry saZipEntry) {
		byte[] swappedBytes = SwapBytes.swapByteArrayToByteArrayForFloats(data);
		ArrayList<Float> floats = SABinary2Float.convertByteArrayToArrayOfFloats(swappedBytes);
		super.createRupVars(floats, false);
		sourceID = BSAFileUtil.getSourceIDFromRuptureVariationZipEntry(saZipEntry, siteName);
		ruptureID = BSAFileUtil.getRuptureIDFromRuptureVariationZipEntry(saZipEntry, siteName);
		rupVar.variationNumber = BSAFileUtil.getRupVarIDFromRuptureVariationZipEntry(saZipEntry, siteName);
		createPeriods(floats, false);
		try {
			rupVar.computeGeomAvg();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public SARuptureFromRuptureVariationFile(byte[] data, String siteName, BSAHeader head) {
		byte[] swappedBytes = SwapBytes.swapByteArrayToByteArrayForFloats(data);
		ArrayList<Float> floats = SABinary2Float.convertByteArrayToArrayOfFloats(swappedBytes);
		super.createRupVars(floats, true);
		sourceID = head.source_id;
		ruptureID = head.rupture_id;
		rupVar.variationNumber = head.rup_var_id;
		createPeriods(floats, true);
		try {
			rupVar.computeGeomAvg();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public SARuptureFromRuptureVariationFile() {
	}

	public void createPeriods(ArrayList<Float> floats, boolean headers) {
		int maxPeriod = SAPeriods.num_periods;
		if (headers) {
			maxPeriod = SAPeriods.num_head_periods;
		}
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

