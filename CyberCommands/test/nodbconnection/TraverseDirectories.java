package test.nodbconnection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.ArrayList;
import junit.framework.JUnit4TestAdapter;

import org.junit.BeforeClass;
import org.junit.Test;

import util.BSAFileUtil;
import util.BSAFilenameFilter;
import util.DirFileFilter;

public class TraverseDirectories {
	
	@BeforeClass public static void runBeforeAllTests() {
		BSAFileUtil.totalFilenameList = new ArrayList<String>();
		BSAFileUtil.totalFileList = new ArrayList<File>();
		File saFile = new File("safiles/traversetest");
		BSAFileUtil.createTotalFileList(saFile);
		
	}
	
	@Test public void pathNames() {
		System.out.println("pathLog: " + BSAFileUtil.pathLog);
		assertEquals("safiles\\traversetesttraversetest safiles\\traversetest\\deeperpathdeeperpath ", BSAFileUtil.pathLog);
	}
	
	@Test public void checkFileExtensions() {
		for (int i=0; i<BSAFileUtil.totalFilenameList.size(); i++) {
			assertTrue(BSAFileUtil.totalFilenameList.get(i).endsWith(".bsa"));
		}
		
	}
	
	@Test public void totalFileArray() {
		assertEquals(BSAFileUtil.totalFilenameList.size(),BSAFileUtil.totalFileList.size());
	}
	
	private static void getSAFiles(File saFile) {
		File[] safilesList = saFile.listFiles(new BSAFilenameFilter());
		File[] sadirsList = saFile.listFiles(new DirFileFilter());
		
		System.out.println("Path: " + saFile.getPath() + saFile.getName());
		if (BSAFileUtil.pathLog == null) {
			BSAFileUtil.pathLog = saFile.getPath() + saFile.getName() + " ";
		}
		else {
			BSAFileUtil.pathLog += saFile.getPath() + saFile.getName() + " ";
		}
		
		
		for (int filesIndex=0; filesIndex < safilesList.length; filesIndex++) {
			System.out.println("\tFilename: " + safilesList[filesIndex].getName());
			if (BSAFileUtil.totalFilenameList == null) {
				System.out.println("null");
			}
			BSAFileUtil.totalFilenameList.add(safilesList[filesIndex].getName());
		}
		
		if (sadirsList != null) {
			for (int dirsIndex=0; dirsIndex<sadirsList.length; dirsIndex++) {
				getSAFiles(sadirsList[dirsIndex]);
			}
		}
		else {
			return;
		}
		
	}

	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter(TraverseDirectories.class);
	}
}
