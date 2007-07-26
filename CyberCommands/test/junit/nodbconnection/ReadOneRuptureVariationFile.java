package test.junit.nodbconnection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;

import junit.framework.JUnit4TestAdapter;

import org.junit.BeforeClass;
import org.junit.Test;

import data.SARuptureFromRuptureVariationFile;
import data.SARuptureVariation;

public class ReadOneRuptureVariationFile {
	private static SARuptureFromRuptureVariationFile saRupture;
	private static boolean showOutputForPeriodValues=false;

	/**
	 * @param args
	 */
	@BeforeClass public static void runBeforeAllTests() {
		//saRupture = new SARuptureFromFile(floats);
		saRupture = new SARuptureFromRuptureVariationFile(new File("safiles/loadampstest/PeakVals_USC_104_12_16.bsa"), "USC");
		//printEastandNorthComponents();
		saRupture.computeAllGeomAvgComponents();
		//printGeomAvgComponents();

	}

	@Test public void testIDsFromFilename() {
		assertEquals(104,saRupture.getSourceID());
		assertEquals(12,saRupture.getRuptureID());
	}

	@Test public void geomAvgComponents() {
		//System.out.println("number of rupture variations: " + saRupture.rupVars.size());
		assertEquals(1,saRupture.rupVars.size());
		for (int rupVarIter=0;rupVarIter<saRupture.rupVars.size();rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);
			SARuptureVariation currRupVar = saRupture.rupVars.get(rupVarIter);
			for (int periodIter=0;periodIter<currRupVar.geomAvgComp.periods.length;periodIter++) {
				if (showOutputForPeriodValues) {
					System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
							+ ", Geometrically Averaged Component, Period " + (periodIter+1) + " : " 
							+ currRupVar.geomAvgComp.periods[periodIter] );
					assertTrue(currRupVar.geomAvgComp.periods[periodIter] > 0);
					assertEquals(16, currRupVar.variationNumber);
				}
			}

		}
	}

	@Test public void eastAndNorthComponents() {
		/*System.out.println("number of rupture variations: " + saRupture.rupVars.size());*/
		assertEquals(1,saRupture.rupVars.size());
		for (int rupVarIter=0;rupVarIter<saRupture.rupVars.size();rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);
			SARuptureVariation currRupVar = saRupture.rupVars.get(rupVarIter);
			for (int periodIter=0;periodIter<currRupVar.eastComp.periods.length;periodIter++) {
				if (showOutputForPeriodValues) {
					System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
							+ ", East Component, Period " + (periodIter+1) + " : " 
							+ currRupVar.eastComp.periods[periodIter] );
					assertTrue(currRupVar.eastComp.periods[periodIter] > 0);
				}
			}
			for (int periodIter=0;periodIter<currRupVar.northComp.periods.length;periodIter++) {
				if (showOutputForPeriodValues) {
					System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
							+ ", North Component, Period " + (periodIter+1) + " : " 
							+ currRupVar.northComp.periods[periodIter] );
					assertTrue(currRupVar.northComp.periods[periodIter] > 0);
					assertEquals(16, currRupVar.variationNumber);
				}
			}

		}
	}

	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (ReadOneRuptureVariationFile.class);
	}
}
