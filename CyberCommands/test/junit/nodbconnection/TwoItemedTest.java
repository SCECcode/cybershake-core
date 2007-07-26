package test.junit.nodbconnection;

import static org.junit.Assert.assertEquals;

import org.junit.BeforeClass;
import org.junit.Test;

import test.fitnesse.TwoItemed;

import junit.framework.JUnit4TestAdapter;

public class TwoItemedTest {
	private static String parseString = "2 3";
	private String multiplyString = "6 6";
	private static TwoItemed twoItemed;
	
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter(TwoItemedTest.class);
	}
	
	@BeforeClass public static void runBeforeAllTests() {
		twoItemed = TwoItemed.parse(parseString);
	}
	
	@Test public void testParseAndToString() {
		System.out.println(twoItemed.toString());
		System.out.println(parseString);
		assertEquals(parseString,twoItemed.toString());
		
		twoItemed.multiply();
		assertEquals(multiplyString, twoItemed.toString());
		
	}
	
	@Test public void testEquals() {
		TwoItemed expectedTwoItemed = new TwoItemed(6,6);
		assertEquals(expectedTwoItemed,twoItemed);
	}
}
