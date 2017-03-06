package test.junit.nodbconnection.fileprocessing;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.BeforeClass;
import org.junit.Test;

import util.BSAFileUtil;

import data.SARuptureFromRuptureVariationFile;

import junit.framework.JUnit4TestAdapter;


public class TestSGTVariationBSAFilename {
	private static final String WNGC_BSA_FILENAME = "safiles/extraletterformat/PeakVals_WNGC_B_464_10_1.bsa";
	private static final String PAS_BSA_FILENAME = "safiles/extraletterformat/PeakVals_PAS_B_21_0_103.bsa";
	private static SARuptureFromRuptureVariationFile saRupturePAS;
	private static SARuptureFromRuptureVariationFile saRuptureWNGC;
	private static File bsaFilePAS;
	private static File bsaFileWNGC;
	
	@BeforeClass public static void setUp() {
		bsaFilePAS = new File(PAS_BSA_FILENAME); 
		bsaFileWNGC = new File(WNGC_BSA_FILENAME);
		saRupturePAS = new SARuptureFromRuptureVariationFile(new File(PAS_BSA_FILENAME), "PAS");
		saRuptureWNGC = new SARuptureFromRuptureVariationFile(new File(WNGC_BSA_FILENAME), "WNGC");
	}
	
	@Test public void testCharacterIsDigit() {
		assertTrue(Character.isDigit("1".charAt(0)));
	}
	
	@Test public void testIDsFromFilenameUsingBSAUtil() {
		BSAFileUtil.setInDebugMode(true);
		assertEquals(21,BSAFileUtil.getSourceIDFromFile(bsaFilePAS, "PAS"));
		assertEquals(0,BSAFileUtil.getRuptureIDFromFile(bsaFilePAS, "PAS"));
		assertEquals(103,BSAFileUtil.getRupVarIDFromRuptureVariationFile(bsaFilePAS, "PAS"));

		assertEquals(464,BSAFileUtil.getSourceIDFromFile(bsaFileWNGC, "WNGC"));
		assertEquals(10,BSAFileUtil.getRuptureIDFromFile(bsaFileWNGC, "WNGC"));
		assertEquals(1,BSAFileUtil.getRupVarIDFromRuptureVariationFile(bsaFileWNGC, "WNGC"));
		
	}
	
	@Test public void testIDsFromFilenameUsingSARupture() {
		assertEquals(21,saRupturePAS.getSourceID());
		assertEquals(0,saRupturePAS.getRuptureID());
		assertEquals(103,saRupturePAS.rupVar.variationNumber);
		
		assertEquals(464,saRuptureWNGC.getSourceID());
		assertEquals(10,saRuptureWNGC.getRuptureID());
		assertEquals(1,saRuptureWNGC.rupVar.variationNumber);
	}
	
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter(TestSGTVariationBSAFilename.class);
	}
}
