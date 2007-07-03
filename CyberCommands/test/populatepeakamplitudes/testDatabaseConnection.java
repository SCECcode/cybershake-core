package test.populatepeakamplitudes;

import mapping.PeakAmplitudes;
import mapping.PeakAmplitudesPK;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.Transaction;
import org.hibernate.cfg.Configuration;

public class testDatabaseConnection {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// 1. Create the PeakAmplitudes Class
		PeakAmplitudes pa = new PeakAmplitudes();
		PeakAmplitudesPK paPK = new PeakAmplitudesPK();
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
        SessionFactory sf = new Configuration().configure().buildSessionFactory();
        
        // 3. Save PeakAmplitude instance and close connection
        Session sess = sf.openSession();
        Transaction t = sess.beginTransaction();
        sess.save(pa);
        t.commit();
        sess.close();

	}

}
