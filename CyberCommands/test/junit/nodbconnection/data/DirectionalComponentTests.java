package test.junit.nodbconnection.data;

import static org.junit.Assert.assertEquals;
import junit.framework.JUnit4TestAdapter;

import org.junit.Test;

import data.DirectionalComponent;

public class DirectionalComponentTests {
	
	private String parseString = "3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3 3.3";

	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (DirectionalComponentTests.class);
	}
	
	@Test public void testParseAndToString() {
		DirectionalComponent dc = DirectionalComponent.parse(parseString);
		assertEquals(dc.toString(), parseString);
	}
	
	@Test public void testEqualsAndNotEquals() {
		DirectionalComponent expectedDC = new DirectionalComponent();
		DirectionalComponent actualDC = new DirectionalComponent();
		
		for (int i=0;i<DirectionalComponent.periodsLength;i++) {
			expectedDC.periods[i] = new Float(5);	
			actualDC.periods[i] = new Float(5);
		}
		
		assertEquals(expectedDC,actualDC);
	}
	
}

