package test.commands;

import static org.junit.Assert.*;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import junit.framework.JUnit4TestAdapter;
import commands.CyberLoadamps;

public class TestCyberLoadamps {
	
	private static final String[] NO_OPTIONS = {""};
	private static final String[] P_OPTION_ONLY = {"-p","safiles"};
	private static final String[] P_SITE_OPTIONS_ONLY  = {"-p","safiles","-site","USC"};
	private static final String[] P_SITE_SGT_OPTIONS_ONLY  = {"-p","safiles","-site","USC","-sgt","1"};
	private static final String[] P_OPTION_MISSING = {"-p"};
	private static final String[] HELP_OPTION = {"-help"}; 
	
	private ByteArrayOutputStream bytes;
	private PrintStream console;
	
	@Before public void setUp() {
		bytes = new ByteArrayOutputStream();
		console = System.out;
		System.setOut(new PrintStream(bytes));
	}
	
	@After public void tearDown() {
		System.setOut(console);
	}
	
	@Test public void testNoArgsMessage() {
		CyberLoadamps.main(NO_OPTIONS);
		assertEquals(CyberLoadamps.getNO_P_OPTION_MESSAGE(),bytes.toString().trim());		
	}
	
	@Test public void testPArgOnly() {
		CyberLoadamps.main(P_OPTION_ONLY);
		assertEquals(CyberLoadamps.getNO_SITE_OPTION_MESSAGE(),bytes.toString().trim());
	}
	
	@Test public void testPandSiteArgsOnly() {
		CyberLoadamps.main(P_SITE_OPTIONS_ONLY);
		assertEquals(CyberLoadamps.getNO_SGT_OPTION_MESSAGE(),bytes.toString().trim());
	}
	
	@Test public void testPandSiteandSGTArgsOnly() {
		CyberLoadamps.main(P_SITE_SGT_OPTIONS_ONLY);
		assertEquals(CyberLoadamps.getNO_SERVER_OPTION_MESSAGE(),bytes.toString().trim());
	}
	
	@Test public void testMissingArgumentException() {
		CyberLoadamps.main(P_OPTION_MISSING);
		assertTrue(bytes.toString().contains("no argument for:"));
	}
	
	@Test public void testHelpOption() {
		CyberLoadamps.main(HELP_OPTION);
		assertTrue(bytes.toString().contains("usage: CyberLoadamps"));
	}
	
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (TestCyberLoadamps.class);
	}
}
