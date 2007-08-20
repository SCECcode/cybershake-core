package test.commands;


import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import junit.framework.JUnit4TestAdapter;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import commands.CyberRegion;

public class TestCyberSites {
	
	private static final String[] HELP_OPTION = {"-help"};
	private static final String[] NO_OPTIONS = {""};
	private static final String[] SERVER_OPTION_ONLY = {"-server","intensity"};

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
	
	@Test public void testHelpOption() {
		CyberRegion.main(HELP_OPTION);
		assertTrue(bytes.toString().contains("usage: CyberRegion"));
	}
	
	@Test public void testNoServerOption() {
		CyberRegion.main(NO_OPTIONS);
		assertEquals(CyberRegion.getNO_SERVER_OPTION_MESSAGE(),bytes.toString().trim());
	}
	
	@Test public void testServerOptionOnly() {
		CyberRegion.main(SERVER_OPTION_ONLY);
		assertEquals(CyberRegion.getNO_SITES_SPECIFIED_MESSAGE(),bytes.toString().trim());
	}

	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (TestCyberSites.class);
	}
	
}
