package test.junit.dbconnection.basicdbops;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Iterator;
import java.util.List;

import junit.framework.JUnit4TestAdapter;
import mapping.CyberShakeSiteRegions;
import mapping.CyberShakeSites;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;
import org.hibernate.criterion.Restrictions;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

public class TestQueryForRegion {
	
	
	
	private static SessionFactory sessFactory;
	private static Session sess;


	@BeforeClass public static void setUp() {
		sessFactory = new Configuration().configure("intensity.cfg.xml").buildSessionFactory();
		sess = sessFactory.openSession();
	}
	
	@Test public void checkConnection() {
		assertTrue(sess.isConnected());
	    assertTrue(sess.isOpen());
	}
	
	@Test public void queryUsingSiteID() {
		List regions = sess.createCriteria(CyberShakeSiteRegions.class).
			add(Restrictions.eq("id.csSiteId", new Integer(18))).
			list();
		
		for (int i=0; i<regions.size(); i++) {
			CyberShakeSiteRegions region = (CyberShakeSiteRegions)regions.get(i);
			assertEquals("USC",region.getCyberShakeSites().getCsShortName().trim());
		}
	}
	
	@Test public void queryUsingHQL() {
		
		List sites = sess.createQuery("from CyberShakeSites").list();
		Iterator result = sess.createQuery("select count(*) from CyberShakeSites").list().iterator();
		Long count = null;
		
		while (result.hasNext()) {
			 count = (Long) result.next();
		}
		
		assertNotNull(count);
		assertEquals(count.intValue(),sites.size());
		
		for (int i=0; i<sites.size(); i ++) {
			CyberShakeSites site = (CyberShakeSites)sites.get(i);
			assertNotNull(site);
		}
	}
	
	@AfterClass public static void tearDown() {
		sess.close();
		assertFalse(sess.isConnected());
	    assertFalse(sess.isOpen());
	}
	
	
	
	
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (TestQueryForRegion.class);
	}
}
