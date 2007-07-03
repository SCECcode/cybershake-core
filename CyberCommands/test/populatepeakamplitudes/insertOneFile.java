package test.populatepeakamplitudes;

import java.util.ArrayList;

import mapping.PeakAmplitudes;
import mapping.PeakAmplitudesPK;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.Transaction;
import org.hibernate.cfg.Configuration;

import util.SABinary2Float;
import util.SwapBytes;
import data.SARuptureVariation;
import data.SAPeriods;
import data.SARuptureFromFile;

public class insertOneFile {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		byte[] byteArrayFromFile = SwapBytes.swapBinaryFileToArrayForFloats("safiles/PeakVals_allUSC_127_6.bsa");
		ArrayList<Float> floats = SABinary2Float.convertByteArrayToArrayOfFloats(byteArrayFromFile);
		SARuptureFromFile saRupture = new SARuptureFromFile(floats);
		//printEastandNorthComponents(saRupture);
		saRupture.computeAllGeomAvgComponents();
		printGeomAvgComponents(saRupture);
		
		insertRupture(saRupture);
		
	}

	private static void insertRupture(SARuptureFromFile saRupture) {
		
		// Fire up Hibernate        
        SessionFactory sf = new Configuration().configure("surface.cfg.xml").buildSessionFactory();
        Session sess = sf.openSession();
        
        SAPeriods saPeriods = new SAPeriods();
    	
        

        
        int outerLoopMax = saRupture.rupVars.size();
		/*System.out.println("number of rupture variations: " + saRupture.rupVars.size());*/
		for (int rupVarIter=0;rupVarIter<outerLoopMax;rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);
			SARuptureVariation currRupVar = saRupture.rupVars.get(rupVarIter);
			int innerLoopMax = currRupVar.geomAvgComp.periods.length;
			for (int periodIter=0;periodIter<innerLoopMax;periodIter++) {
				/*System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
						+ ", Geometrically Averaged Component, Period " + (periodIter+1) + " : " 
						+ currRupVar.geomAvgComp.periods[periodIter] );*/
				
				// 1. Create the PeakAmplitudes Class
				PeakAmplitudes pa = new PeakAmplitudes();
				PeakAmplitudesPK paPK = new PeakAmplitudesPK();
				paPK.setERF_ID(29);
				paPK.setSite_ID(1);
				paPK.setSource_ID(127);
				paPK.setRupture_ID(6);
				paPK.setRup_Var_ID(currRupVar.variationNumber);
				paPK.setRup_Var_Scenario_ID(1);
				paPK.setIM_Type(new String("SA Period " + saPeriods.getNextValue()));
				pa.setPaPK(paPK);
				pa.setIM_Value(currRupVar.geomAvgComp.periods[periodIter]);
				pa.setUnits("meters per second squared");
		        
		
		        sess.beginTransaction();
				
		        // 3. Save and Commit PeakAmplitude instance		        
		        sess.save(pa);
				
		        sess.getTransaction().commit();
		        
		        
			}
			
		}
		
		
		sess.close();

		
		
		
	}
	
	private static void printGeomAvgComponents(SARuptureFromFile saRupture) {
		System.out.println("number of rupture variations: " + saRupture.rupVars.size());
		for (int rupVarIter=0;rupVarIter<saRupture.rupVars.size();rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);
			SARuptureVariation currRupVar = saRupture.rupVars.get(rupVarIter);
			for (int periodIter=0;periodIter<currRupVar.geomAvgComp.periods.length;periodIter++) {
				System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
						+ ", Geometrically Averaged Component, Period " + (periodIter+1) + " : " 
						+ currRupVar.geomAvgComp.periods[periodIter] );
			}
			
		}
	}

	private static void printEastandNorthComponents(SARuptureFromFile saRupture) {
		System.out.println("number of rupture variations: " + saRupture.rupVars.size());
		for (int rupVarIter=0;rupVarIter<saRupture.rupVars.size();rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);
			SARuptureVariation currRupVar = saRupture.rupVars.get(rupVarIter);
			for (int periodIter=0;periodIter<currRupVar.eastComp.periods.length;periodIter++) {
				System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
						+ ", East Component, Period " + (periodIter+1) + " : " 
						+ currRupVar.eastComp.periods[periodIter] );
			}
			for (int periodIter=0;periodIter<currRupVar.northComp.periods.length;periodIter++) {
				System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
						+ ", North Component, Period " + (periodIter+1) + " : " 
						+ currRupVar.northComp.periods[periodIter] );
			}
			
		}
	}
}
