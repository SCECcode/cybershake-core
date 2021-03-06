package test.junit.dbconnection.basicdbops;

import junit.framework.JUnit4TestAdapter;
import mapping.PeakAmplitude;
import mapping.PeakAmplitudePK;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;
import org.junit.Test;

public class IntensityDBAndInsertOnePA {

	/**
	 * @param args
	 */
	@Test public void testIntensityDBAndInsertOnePA() {
		// 1. Create the PeakAmplitudes Class
		PeakAmplitude pa = new PeakAmplitude();
		PeakAmplitudePK paPK = new PeakAmplitudePK();
		paPK.setERF_ID(1);
		paPK.setSite_ID(1);
		paPK.setSource_ID(1);
		paPK.setRupture_ID(1);
		paPK.setRup_Var_ID(1);
		paPK.setRup_Var_Scenario_ID(1);
		paPK.setIM_Type("a");
		pa.setPaPK(paPK);
		pa.setIM_Value(1);
		pa.setUnits("feet");
		
		// 2. Fire up Hibernate        
        SessionFactory sf = new Configuration().configure("intensity.cfg.xml").buildSessionFactory();
        
        // 3. Save PeakAmplitude instance and close connection
        Session sess = sf.openSession();
        //Transaction t = sess.beginTransaction();
        sess.save(pa);
        //t.commit();
        sess.close();

	}
	
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (IntensityDBAndInsertOnePA.class);
	}

}
