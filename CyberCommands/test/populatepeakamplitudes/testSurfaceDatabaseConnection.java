package test.populatepeakamplitudes;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.Transaction;
import org.hibernate.cfg.Configuration;

public class testSurfaceDatabaseConnection {
	public static void main(String[] args) {
		// 2. Fire up Hibernate        
	    SessionFactory sf = new Configuration().configure("surface.cfg.xml").buildSessionFactory();
	    
	    // 3. Save PeakAmplitude instance and close connection
	    Session sess = sf.openSession();
	    sess.close();		
	}

}
