package data;

import java.io.File;
import java.util.ArrayList;

import util.BSAFileUtil;
import util.SABinary2Float;
import util.SwapBytes;

public class SARuptureFromFile {

	public int numRuptureVariations = 0;
	public String site = "";
	protected int sourceID = -1;
	protected int ruptureID = -1;
	public ArrayList<SARuptureVariation> rupVars = new ArrayList<SARuptureVariation>();

	public SARuptureFromFile() {
		
	}
	
	public SARuptureFromFile(ArrayList<Float> floats) {
		createRupVars(floats, false);
	}

	public SARuptureFromFile(File bsaFile, String siteName) {
		byte[] byteArrayFromFile = SwapBytes.swapBinaryFileToArrayForFloats(bsaFile.getPath());
		ArrayList<Float> floats = SABinary2Float.convertByteArrayToArrayOfFloats(byteArrayFromFile);
		createRupVars(floats, false);
		sourceID = BSAFileUtil.getSourceIDFromFile(bsaFile, siteName);
		ruptureID = BSAFileUtil.getRuptureIDFromFile(bsaFile, siteName);
	}

	protected void createRupVars(ArrayList<Float> floats, boolean header) {
		int maxPeriod = SAPeriods.num_periods;
		if (header) {
			maxPeriod = SAPeriods.num_head_periods;
		}
		int ruptureVarCount = 0;
		int periodCount = 1;
		SARuptureVariation rupVar = new SARuptureVariation();
		rupVar.variationNumber = 0;
		String componentDirection = "East";
		for (int i=0; i<floats.size(); i++) {
			if (periodCount > maxPeriod) {
				periodCount = 1;
				if (componentDirection.equals("East")) {
					componentDirection = "North";
				}
				else if (componentDirection.equals("North")){
					// Switch component direction
					componentDirection = "East";
					
					// Add the rupture variation to the rupture variation list
					/*System.out.println("Adding rupVar " + ruptureVarCount + " to rupVars: ");*/
					rupVars.add(rupVar);
					
					// Create a new rupture variation
					rupVar = new SARuptureVariation();
					ruptureVarCount++;
					rupVar.variationNumber = ruptureVarCount;
				}

			}
			/*System.out.println("SA for Rupture Variation " + ruptureVarCount + ", " + componentDirection + " Component, Period " + periodCount + " : " + floats.get(i) );*/
			if (componentDirection.equals("North")) {
				rupVar.northComp.periods[periodCount-1] = floats.get(i);
			}
			else if (componentDirection.equals("East")) {
				rupVar.eastComp.periods[periodCount-1] = floats.get(i);
			}
			periodCount++;
		}
		// Add the final rupture variation to the rupture variation list 
		/*System.out.println("Adding rupVar " + ruptureVarCount + " to rupVars: ");*/
		rupVars.add(rupVar);
		
		numRuptureVariations = ruptureVarCount+1;
		
	}

	public void createSARuptureFromArray(ArrayList<Float> floats) {
		int maxPeriod = SAPeriods.num_periods;
		int ruptureVarCount = 1;
		int periodCount = 1;
		String componentDirection = "East";
		for (int i=0; i<floats.size(); i++) {
			if (periodCount > maxPeriod) {
				periodCount = 1;
				if (componentDirection.equals("East")) {
					componentDirection = "North";
				}
				else if (componentDirection.equals("North")){
					componentDirection = "East";
					ruptureVarCount++;
				}

			}
			/*System.out.println("SA for Rupture Variation " + ruptureVarCount + ", " + componentDirection + " Component, Period " + periodCount + " : " + floats.get(i) );*/
			periodCount++;
		}
		numRuptureVariations = ruptureVarCount;
	}

	public void computeAllGeomAvgComponents() {
		for (int i=0;i<rupVars.size();i++){
			try {
				rupVars.get(i).computeGeomAvg();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	public int getRuptureID() {
		return ruptureID;
	}

	public void setRuptureID(int ruptureID) {
		this.ruptureID = ruptureID;
	}

	public int getSourceID() {
		return sourceID;
	}

	public void setSourceID(int sourceID) {
		this.sourceID = sourceID;
	}
}
