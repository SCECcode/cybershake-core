package test.nonautomated;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import mapping.PeakAmplitude;
import mapping.PeakAmplitudePK;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;

import util.BSAFilenameFilter;
import data.SAPeriods;
import data.SARuptureFromFile;
import data.SARuptureVariation;

public class insertAllUSCFiles {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		BufferedWriter filelistwriter = new BufferedWriter(new FileWriter("filelist"));


		File safiles = new File("safiles/lbp");
		File[] safilesList =  safiles.listFiles(new BSAFilenameFilter());
		//SARuptureFromFile saRupture = new SARuptureFromFile(floats);

		//		 Fire up Hibernate        
		SessionFactory sf = new Configuration().configure("surface.cfg.xml").buildSessionFactory();
		Session sess = sf.openSession();

		for (int i=0; i<safilesList.length; i++) {

			writeFileListToFile(filelistwriter, safilesList, i);

			//System.out.println("Filename: " + safilesList[i].getName());
			prepAndExecuteRuptureInsertion(safilesList, sess, i, "LBP");
		}

		sess.close();


	}

	private static void writeFileListToFile(BufferedWriter filelistwriter, File[] safilesList, int i) throws IOException {
		filelistwriter.write("Filename: " + safilesList[i].getName());
		filelistwriter.newLine();
	}

	private static void prepAndExecuteRuptureInsertion(File[] safilesList, Session sess, int i, String siteName) {
		SARuptureFromFile saRupture = new SARuptureFromFile(safilesList[i], siteName);
		//printEastandNorthComponents(saRupture);
		saRupture.computeAllGeomAvgComponents();
		//printGeomAvgComponents(saRupture);

		insertRupture(saRupture, sess);
	}

	private static void insertRupture(SARuptureFromFile saRupture, Session sess) {

		SAPeriods saPeriods = new SAPeriods();
		
		
		
		int outerLoopMax = saRupture.rupVars.size();
		/*System.out.println("number of rupture variations: " + saRupture.rupVars.size());*/
		for (int rupVarIter=0;rupVarIter<outerLoopMax;rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);
			
			// Do garbage collection
			if ((rupVarIter+1)%250==0) System.gc();
			int currentERF_ID = 29;
			int currentSource_ID = saRupture.getSourceID();
			int currentRupture_ID = saRupture.getRuptureID();
			int currentRup_Var_Scenario_ID = 1;

			SARuptureVariation currRupVar = saRupture.rupVars.get(rupVarIter);

/*			// Initialize RuptureVariation class
			RuptureVariation rv = new RuptureVariation();
			RuptureVariationPK rupVarPK = new RuptureVariationPK();
			// Create the RuptureVariation Class
			rupVarPK.setERF_ID(currentERF_ID );
			rupVarPK.setSource_ID(currentSource_ID);
			rupVarPK.setRupture_ID(currentRupture_ID);
			rupVarPK.setRup_Var_ID(currRupVar.variationNumber);
			rupVarPK.setRup_Var_Scenario_ID(currentRup_Var_Scenario_ID);
			rv.setRupVarPK(rupVarPK);
			rv.setRup_Var_LFN("LFN.file");*/

			sess.beginTransaction();

			/*sess.save(rv);*/

			//sess.getTransaction().commit();

			

			int innerLoopMax = currRupVar.geomAvgComp.periods.length;
			for (int periodIter=0;periodIter<innerLoopMax;periodIter++) {
				/*System.out.println("SA for Rupture Variation " + currRupVar.variationNumber 
						+ ", Geometrically Averaged Component, Period " + (periodIter+1) + " : " 
						+ currRupVar.geomAvgComp.periods[periodIter] );*/

				// Initialize PeakAmplitudes class
				PeakAmplitude pa = new PeakAmplitude();
				PeakAmplitudePK paPK = new PeakAmplitudePK();
				// Set values for the PeakAmplitudes Class
				paPK.setERF_ID(currentERF_ID);
				paPK.setSite_ID(11);				
				paPK.setSource_ID(currentSource_ID);				
				paPK.setRupture_ID(currentRupture_ID);
				paPK.setRup_Var_ID(currRupVar.variationNumber);
				paPK.setRup_Var_Scenario_ID(currentRup_Var_Scenario_ID);
				paPK.setIM_Type(new String("SA_Period_" + saPeriods.getNextValue()));
				pa.setPaPK(paPK);
				pa.setIM_Value(currRupVar.geomAvgComp.periods[periodIter]);
				pa.setUnits("meters per second squared");


				// 3. Save and Commit PeakAmplitude instance

				//sess.beginTransaction();

				sess.save(pa);
				if ((rupVarIter+1)%50==0) {
					// flush a batch of inserts and release memory
					sess.flush();
					sess.clear();
				}
			}
			
			sess.getTransaction().commit();

		}


	}

	@SuppressWarnings("unused")
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

	@SuppressWarnings("unused")
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
