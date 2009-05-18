package test.junit.dbconnection.basicdbops;

import org.junit.Test;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import junit.framework.JUnit4TestAdapter;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;


public class IntensityDBConnection {
	@Test public void intensityDatabaseConnection() {
		// 2. Fire up Hibernate        
	    SessionFactory sf = new Configuration().configure("intensity.cfg.xml").buildSessionFactory();
	    
	    // 3. Save PeakAmplitude instance and close connection
	    Session sess = sf.openSession();
	    
	    assertTrue(sess.isConnected());
	    assertTrue(sess.isOpen());
	    
	    sess.close();		
	    
	    assertFalse(sess.isConnected());
	    assertFalse(sess.isOpen());
	}
	
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (SurfaceDBConnection.class);
	}
}
