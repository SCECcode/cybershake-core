package test.junit.processing;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.junit.Test;

import processing.RuptureVariationFileInserter;

import junit.framework.JUnit4TestAdapter;

public class TestRuptureVariationFileInserter {
	
	private static final String INTENSITY = "intensity";
	private static final String SGT_VARIATION_ID_1 = "1";
	private static final String SGT_VARIATION_ID_523 = "523";
	private static final String SITENAME = "USC";
	private static final String PATHNAME = "safiles";

	@Test public void testConstructor() {
		try {
			RuptureVariationFileInserter rvfi = new RuptureVariationFileInserter(PATHNAME, SITENAME, SGT_VARIATION_ID_1, INTENSITY);
			assertEquals(Integer.parseInt(SGT_VARIATION_ID_1),rvfi.getCurrentSGT_Variation_ID());
			assertEquals(PATHNAME, rvfi.getPathName());
			assertEquals(SITENAME, rvfi.getSiteName());
			
			rvfi = new RuptureVariationFileInserter(PATHNAME, SITENAME, SGT_VARIATION_ID_523, INTENSITY);
			assertEquals(Integer.parseInt(SGT_VARIATION_ID_523),rvfi.getCurrentSGT_Variation_ID());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (TestRuptureVariationFileInserter.class);
	}
	
}
