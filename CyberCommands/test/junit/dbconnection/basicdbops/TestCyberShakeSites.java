package test.junit.dbconnection.basicdbops;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;
import org.hibernate.classic.Session;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import mapping.CyberShakeSites;
import junit.framework.JUnit4TestAdapter;


public class TestCyberShakeSites {
	
	private static Session sess;
	private static SessionFactory sf;

	@BeforeClass public static void runBeforeAllTests() {
		sf = new Configuration().configure("intensity.cfg.xml").buildSessionFactory();
		sess = sf.openSession();
		assertTrue(sess.isConnected());
		assertTrue(sess.isOpen());
	}
	
	@Test public void testComparison() {
		
		CyberShakeSites siteUSC = (CyberShakeSites) sess.createQuery(
				"from CyberShakeSites as sites where sites.csShortName = 'USC'")
				.uniqueResult();
		
		CyberShakeSites siteFromComparison = (CyberShakeSites) sess.createQuery(
				"from CyberShakeSites as sites where sites = ?")
				.setEntity(0, siteUSC)
				.uniqueResult();
		
		assertEquals(siteUSC, siteFromComparison);
				
	}
	
	@AfterClass public static void tearDown() {
		sess.close();
		assertFalse(sess.isConnected());
		assertFalse(sess.isOpen());
	}

	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (TestCyberShakeSites.class);
	}

}
