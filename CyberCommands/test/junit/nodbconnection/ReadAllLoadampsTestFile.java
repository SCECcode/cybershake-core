package test.junit.nodbconnection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.ArrayList;

import org.junit.BeforeClass;
import org.junit.Test;

import data.SARuptureFromRuptureVariationFile;
import data.SARuptureVariation;

import junit.framework.JUnit4TestAdapter;

import util.BSAFileUtil;

public class ReadAllLoadampsTestFile {

	private static ArrayList<File> totalFileList;
	private static boolean showOutputForPeriodValues = true;

	@BeforeClass public static void runBeforeAllTests() {
		totalFileList = BSAFileUtil.createTotalFileList(new File("safiles/loadampstest"));
	}

	@Test public void loadAllSARuptureVariations() {
		for (int i=0; i<totalFileList.size(); i++) {
			SARuptureFromRuptureVariationFile saRupture = new SARuptureFromRuptureVariationFile(totalFileList.get(i),"USC");
			System.out.println("Source_ID: " + saRupture.getSourceID() + ", Rupture_ID: " + saRupture.getRuptureID());
			eastAndNorthComponents(saRupture);
			saRupture.computeAllGeomAvgComponents();
			geomAvgComponents(saRupture);
		}
	}

	public void eastAndNorthComponents(SARuptureFromRuptureVariationFile saRupture) {
		/*System.out.println("number of rupture variations: " + saRupture.rupVars.size());*/
		assertEquals(1,saRupture.rupVars.size());
		for (int rupVarIter=0;rupVarIter<saRupture.rupVars.size();rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);
			SARuptureVariation currRupVar = saRupture.rupVars.get(rupVarIter);
			for (int periodIter=0;periodIter<currRupVar.eastComp.periods.length;periodIter++) {
				if (showOutputForPeriodValues) {
					System.out.println("Rupture Variation " + currRupVar.variationNumber 
							+ ", East Component, Period " + (periodIter+1) + " : " 
							+ currRupVar.eastComp.periods[periodIter] );
					assertTrue(currRupVar.eastComp.periods[periodIter] > 0);
				}
			}
			for (int periodIter=0;periodIter<currRupVar.northComp.periods.length;periodIter++) {
				if (showOutputForPeriodValues) {
					System.out.println("Rupture Variation " + currRupVar.variationNumber 
							+ ", North Component, Period " + (periodIter+1) + " : " 
							+ currRupVar.northComp.periods[periodIter] );
					assertTrue(currRupVar.northComp.periods[periodIter] > 0);
				}
			}

		}
	}

	public void geomAvgComponents(SARuptureFromRuptureVariationFile saRupture) {
		//System.out.println("number of rupture variations: " + saRupture.rupVars.size());
		assertEquals(1,saRupture.rupVars.size());
		for (int rupVarIter=0;rupVarIter<saRupture.rupVars.size();rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);
			SARuptureVariation currRupVar = saRupture.rupVars.get(rupVarIter);
			for (int periodIter=0;periodIter<currRupVar.geomAvgComp.periods.length;periodIter++) {
				if (showOutputForPeriodValues) {
					System.out.println("Rupture Variation " + currRupVar.variationNumber 
							+ ", Geometrically Averaged Component, Period " + (periodIter+1) + " : " 
							+ currRupVar.geomAvgComp.periods[periodIter] );
					assertTrue(currRupVar.geomAvgComp.periods[periodIter] > 0);
				}
			}

		}
	}




	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (ReadAllLoadampsTestFile.class);
	}

}
