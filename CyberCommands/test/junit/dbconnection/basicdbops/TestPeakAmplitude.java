	package test.junit.dbconnection.basicdbops;
	
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;

import junit.framework.JUnit4TestAdapter;
import mapping.PeakAmplitude;
import mapping.PeakAmplitudePK;

import org.hibernate.Criteria;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;
import org.hibernate.criterion.Restrictions;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import data.SAPeriods;
import data.SARuptureFromRuptureVariationFile;

	public class TestPeakAmplitude {
		
		private static SessionFactory sf;
		private static Session sess;
		private static PeakAmplitude expectedPA;
		private static SARuptureFromRuptureVariationFile saRupture;
		private static SAPeriods saperiods;
		private boolean isInDebugMode = false;

		@BeforeClass public static void runBeforeAllTests() {
			// Open connection
			sf = new Configuration().configure("intensity.cfg.xml").buildSessionFactory();
			sess = sf.openSession();
			assertTrue(sess.isConnected());
			assertTrue(sess.isOpen());
			
			// Initialize data objects
			saperiods = new SAPeriods();
			saRupture = new SARuptureFromRuptureVariationFile(new File("safiles/loadampstest/PeakVals_USC_104_12_3.bsa"),"USC");
			expectedPA = new PeakAmplitude();
			PeakAmplitudePK paPK = new PeakAmplitudePK();
			
			// Set values for the PeakAmplitudes Class
			paPK.setERF_ID(29);
			paPK.setSite_ID(18);
			assertEquals(104,saRupture.getSourceID());
			assertEquals(12,saRupture.getRuptureID());
			assertEquals(3,saRupture.rupVar.variationNumber);
			paPK.setSource_ID(saRupture.getSourceID());				
			paPK.setRupture_ID(saRupture.getRuptureID());
			paPK.setRup_Var_ID(saRupture.rupVar.variationNumber);
			paPK.setRup_Var_Scenario_ID(1);
			paPK.setSGT_Variation_ID(1);
			paPK.setIM_Type(new String("SA_Period_" + saperiods.getNextValue()));
			expectedPA.setPaPK(paPK);
			expectedPA.setIM_Value(saRupture.rupVar.geomAvgComp.periods[0]);
			expectedPA.setUnits("cm per second squared");
		}
		
		/*@Test public void useExample() {
			Criteria crit = sess.createCriteria(PeakAmplitude.class);
			crit.setMaxResults(50);
			Example examplePA = Example.create(expectedPA)
			.excludeProperty("IM_Value");
			crit.add(examplePA);
			PeakAmplitude resultsPA = (PeakAmplitude)crit.uniqueResult();
			//List peakamps = crit.list();
			
			PeakAmplitude resultsPA = (PeakAmplitude) sess.createQuery(
		    "from PeakAmplitudes as pa where pa = ?")
		    .setEntity(0, expectedPA)
		    .uniqueResult();
			
			System.out.println("expectedPA: " + expectedPA);
			
			PeakAmplitude resultantPA = (PeakAmplitude) sess.createQuery(
					"from PeakAmplitudes as pa where pa = ?")
					.setEntity(0,expectedPA)
					.uniqueResult();

			
			//PeakAmplitude resultsPA = (PeakAmplitude)peakamps.get(0);
			if (isInDebugMode) {
				System.out.println("-- results from example --");
				System.out.println(resultsPA);
			}
			
			assertEquals(expectedPA, resultsPA);
		}*/
		
		@Test public void useHQL() {
			String hqlQuery = "from PeakAmplitude where ERF_ID = 29 AND Site_ID = 18 AND Rup_Var_Scenario_ID = 1 AND Source_ID = ? AND Rupture_ID = ? AND Rup_Var_ID = ? AND SGT_Variation_ID = 1 AND IM_Type = ?";
			
			PeakAmplitude resultantPA = (PeakAmplitude) sess.createQuery(hqlQuery)
				.setInteger(0, saRupture.getSourceID())
				.setInteger(1, saRupture.getRuptureID())
				.setInteger(2, saRupture.rupVar.variationNumber)
				.setString(3, "SA_Period_" + saperiods.getCurrentValue())
				.uniqueResult();
			
			if (isInDebugMode) {
				System.out.println("-- results from hql query --");
				System.out.println(resultantPA);	
			}
			
			assertEquals(expectedPA, resultantPA);
		}
		
		@Test public void useSQLRestriction() {
			Criteria crit = sess.createCriteria(PeakAmplitude.class);
			crit.setMaxResults(50);
			String sqlRestriction = "ERF_ID = 29 AND Site_ID = 18 AND Rup_Var_Scenario_ID = 1 AND Source_ID = "  + saRupture.getSourceID() 
									+ " AND Rupture_ID = " + saRupture.getRuptureID() + " AND Rup_Var_ID = " + saRupture.rupVar.variationNumber
									+ " AND SGT_Variation_ID = 1 AND IM_Type = '" + "SA_Period_" + saperiods.getCurrentValue() + "' ";
									
			crit.add( Restrictions.sqlRestriction(sqlRestriction));
			List peakamps = crit.list();
			
			PeakAmplitude resultsPA = (PeakAmplitude)peakamps.get(0);
			if (isInDebugMode) {
				System.out.println("-- results from sqlrestriction --");
				System.out.println(resultsPA);	
			}
			
			assertEquals(expectedPA, resultsPA);
		}
		
		@AfterClass public static void runAfterAllTests() {
			sess.close();
			assertFalse(sess.isConnected());
			assertFalse(sess.isOpen());
		}
		
		public static junit.framework.Test suite() {
			return new JUnit4TestAdapter (TestPeakAmplitude.class);
		}
	}
