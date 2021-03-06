package test.junit.nodbconnection.fileprocessing;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;

import junit.framework.JUnit4TestAdapter;

import org.junit.BeforeClass;
import org.junit.Test;

import util.BSAFileUtil;

public class TraverseDirectories {
	
	@BeforeClass public static void runBeforeAllTests() {
		File saFile = new File("safiles/traversetest");
		BSAFileUtil.setInDebugMode(false);
		BSAFileUtil.createTotalFileList(saFile);
		
		
	}
	
	@Test public void pathNames() {
		//System.out.println("pathLog: " + BSAFileUtil.pathLog);
		assertEquals("safiles\\traversetest safiles\\traversetest\\deeperpath ", BSAFileUtil.pathLog);
	}
	
	@Test public void checkFileExtensions() {
		for (int i=0; i<BSAFileUtil.totalFilenameList.size(); i++) {
			assertTrue(BSAFileUtil.totalFilenameList.get(i).endsWith(".bsa"));
		}
		
	}
	
	@Test public void totalFileArray() {
		assertEquals(BSAFileUtil.totalFilenameList.size(),BSAFileUtil.totalFileList.size());
	}

	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter(TraverseDirectories.class);
	}
}
