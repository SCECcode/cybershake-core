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
	private static SARuptureFromRuptureVariationFile saRupture;
	private static File bsaFile;
	
	@BeforeClass public static void setUp() {
		bsaFile = new File("safiles/extraletterformat/PeakVals_PAS_B_21_0_103.bsa"); 
		saRupture = new SARuptureFromRuptureVariationFile(new File("safiles/extraletterformat/PeakVals_PAS_B_21_0_103.bsa"), "PAS");		
	}
	
	@Test public void testCharacterIsDigit() {
		assertTrue(Character.isDigit("1".charAt(0)));
	}
	
	@Test public void testIDsFromFilenameUsingBSAUtil() {
		BSAFileUtil.setInDebugMode(true);
		assertEquals(21,BSAFileUtil.getSourceIDFromFile(bsaFile, "PAS"));
		assertEquals(0,BSAFileUtil.getRuptureIDFromFile(bsaFile, "PAS"));
		assertEquals(103,BSAFileUtil.getRupVarIDFromRuptureVariationFile(bsaFile, "PAS"));
	}
	
	@Test public void testIDsFromFilenameUsingSARupture() {
		assertEquals(21,saRupture.getSourceID());
		assertEquals(0,saRupture.getRuptureID());
		assertEquals(103,saRupture.rupVar.variationNumber);
	}
	
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter(TestSGTVariationBSAFilename.class);
	}
}
