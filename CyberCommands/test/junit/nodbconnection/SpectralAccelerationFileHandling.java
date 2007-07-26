package test.junit.nodbconnection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;

import junit.framework.JUnit4TestAdapter;

import org.junit.BeforeClass;
import org.junit.Test;
import util.BSAFilenameFilter;
import util.BSAFileUtil;


public class SpectralAccelerationFileHandling {

	private static File safiles;
	private static File[] safilesList;
	private static String siteName;

	@BeforeClass public static void runOnceBeforeAllTests() {
		siteName = "USC";
		safiles = new File("safiles/loadampstest");
		safilesList = safiles.listFiles(new BSAFilenameFilter());
	}
	
	@Test public void fileListTest() {
		assertNotNull(safiles.listFiles());
	}
	
	@Test public void bsaFileFilterNotNull() {
		assertNotNull(safilesList);
	}

	@Test public void bsaFileUtilForRuptureVariationsFileTest() {
		File bsaFile = new File("safiles/rupturevariations/PeakVals_USC_413_13_8.bsa");
		int source_ID = BSAFileUtil.getSourceIDFromRuptureVariationFile(bsaFile, siteName);
		int rupture_ID = BSAFileUtil.getRuptureIDFromRuptureVariationFile(bsaFile, siteName);
		int rupVar_ID = BSAFileUtil.getRupVarIDFromRuptureVariationFile(bsaFile, siteName);
		assertEquals(413, source_ID);
		assertEquals(13,rupture_ID);
		assertEquals(8,rupVar_ID);
	}

	@Test public void bsaFileUtilTest() {
		File bsaFile = new File("safiles/PeakVals_allUSC_127_6.bsa");
		int source_ID = BSAFileUtil.getSourceIDFromFile(bsaFile, siteName);
		int rupture_ID = BSAFileUtil.getRuptureIDFromFile(bsaFile, siteName);
		assertEquals(127, source_ID);
		assertEquals(6,rupture_ID);
	}
	
	@Test public void printAllFilesinFileArray() {
		
		for (int i=0; i<safilesList.length; i++) {
			assertTrue(safilesList[i].getName().endsWith(".bsa"));
			//System.out.println("Filename: " + safilesList[i].getName());
		}
	}

	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (SpectralAccelerationFileHandling.class);
	}

}
