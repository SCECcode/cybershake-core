package test.junit.dbconnection.basicdbops;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.hibernate.Hibernate;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;
import org.junit.Test;

import junit.framework.JUnit4TestAdapter;

public class SQLQueryForSiteID {
	
	@Test public void sqlQueryForSiteID() {
		SessionFactory sessFactory = new Configuration().configure("intensity.cfg.xml").buildSessionFactory();
		Session retrieveSiteIDSess = sessFactory.openSession();
		
		String siteIDQuery = "SELECT CS_Site_ID FROM CyberShake_Sites WHERE CS_Site_Name = 'USC'";
		//System.out.println(query);
		List siteIDList = retrieveSiteIDSess.createSQLQuery(siteIDQuery).addScalar("CS_Site_ID", Hibernate.INTEGER).list();
		Object siteIDObject = siteIDList.get(0); 
		int siteID = (Integer)siteIDObject;
		
		retrieveSiteIDSess.close();
		
		assertEquals(18,siteID);
	}
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (SQLQueryForSiteID.class);
	}
}
