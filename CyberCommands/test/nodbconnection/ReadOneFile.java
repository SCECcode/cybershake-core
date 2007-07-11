package test.nodbconnection;

import java.io.File;
import java.util.ArrayList;

import junit.framework.JUnit4TestAdapter;

import data.SARuptureVariation;
import data.SARuptureFromFile;

import util.SABinary2Float;
import util.SwapBytes;

import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class ReadOneFile {
	private static SARuptureFromFile saRupture;

	/**
	 * @param args
	 */
	@BeforeClass public static void runBeforeAllTests() {
		byte[] byteArrayFromFile = SwapBytes.swapBinaryFileToArrayForFloats("safiles/norupvarids/PeakVals_allUSC_100_0.bsa");
		ArrayList<Float> floats = SABinary2Float.convertByteArrayToArrayOfFloats(byteArrayFromFile);
		//saRupture = new SARuptureFromFile(floats);
		saRupture = new SARuptureFromFile(new File("safiles/norupvarids/PeakVals_allUSC_100_0.bsa"), "USC");
		//printEastandNorthComponents();
		saRupture.computeAllGeomAvgComponents();
		//printGeomAvgComponents();
		
	}
	
	@Test public void testIDsFromFilename() {
		assertEquals(100,saRupture.getSourceID());
		assertEquals(0,saRupture.getRuptureID());
	}

	@Test public void geomAvgComponents() {
		//System.out.println("number of rupture variations: " + saRupture.rupVars.size());
		/*assertEquals(8,saRupture.rupVars.size());*/
		for (int rupVarIter=0;rupVarIter<saRupture.rupVars.size();rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);
			SARuptureVariation currRupVar = saRupture.rupVars.get(rupVarIter);
			for (int periodIter=0;periodIter<currRupVar.geomAvgComp.periods.length;periodIter++) {
				/*System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
						+ ", Geometrically Averaged Component, Period " + (periodIter+1) + " : " 
						+ currRupVar.geomAvgComp.periods[periodIter] );*/
				assertTrue(currRupVar.geomAvgComp.periods[periodIter] > 0);
			}
			
		}
	}

	@Test public void eastandNorthComponents() {
		/*System.out.println("number of rupture variations: " + saRupture.rupVars.size());*/
		assertEquals(8,saRupture.rupVars.size());
		for (int rupVarIter=0;rupVarIter<saRupture.rupVars.size();rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);
			SARuptureVariation currRupVar = saRupture.rupVars.get(rupVarIter);
			for (int periodIter=0;periodIter<currRupVar.eastComp.periods.length;periodIter++) {
				/*System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
						+ ", East Component, Period " + (periodIter+1) + " : " 
						+ currRupVar.eastComp.periods[periodIter] );*/
				assertTrue(currRupVar.eastComp.periods[periodIter] > 0);
			}
			for (int periodIter=0;periodIter<currRupVar.northComp.periods.length;periodIter++) {
				/*System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
						+ ", North Component, Period " + (periodIter+1) + " : " 
						+ currRupVar.northComp.periods[periodIter] );*/
				assertTrue(currRupVar.northComp.periods[periodIter] > 0);
			}
			
		}
	}
	
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (ReadOneFile.class);
	}
}
