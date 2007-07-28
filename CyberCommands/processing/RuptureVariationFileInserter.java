package processing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import mapping.PeakAmplitude;
import mapping.PeakAmplitudePK;

import org.hibernate.Criteria;
import org.hibernate.Hibernate;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;
import org.hibernate.criterion.Example;

import util.BSAFileUtil;
import data.SAPeriods;
import data.SARuptureFromRuptureVariationFile;
import data.SARuptureVariation;

public class RuptureVariationFileInserter {
	
	public static int siteID = 18;
	private static ArrayList<File> totalFilesList;
	private String siteIDQuery;
	private SessionFactory sessFactory;
	private static String siteName;
	
	
	public RuptureVariationFileInserter(String pathName, String newSiteName) throws IOException {
		sessFactory = new Configuration().configure("intensity.cfg.xml").buildSessionFactory();
		Session retrieveSiteIDSess = sessFactory.openSession();
		
		siteName = newSiteName;
		siteIDQuery = "SELECT CS_Site_ID FROM CyberShake_Sites WHERE CS_Site_Name = '" + siteName + "'";
		//System.out.println(query);
		List siteIDList = retrieveSiteIDSess.createSQLQuery(siteIDQuery).addScalar("CS_Site_ID", Hibernate.INTEGER).list();
		Object siteIDObject = siteIDList.get(0); 
		siteID = (Integer)siteIDObject;
		
		//System.out.println("Site_ID for " + siteName + ": " + siteID);
		
		retrieveSiteIDSess.close();
		
		BSAFileUtil.totalFilenameList = new ArrayList<String>();
		BSAFileUtil.totalFileList = new ArrayList<File>();
		File saFile = new File(pathName);
		
		totalFilesList = BSAFileUtil.createTotalFileList(saFile);
		
		Session sess = sessFactory.openSession();
		insertAllRuptureVariationFiles(sess);
		sess.getTransaction().commit();
		sess.close();


	}

	private void insertAllRuptureVariationFiles(Session sess) {
		for (int i=0; i<totalFilesList.size(); i++) {

			//writeFileListToFile(filelistwriter, safilesList, i);

			//System.out.println("Filename: " + safilesList[i].getName());
			prepAndExecuteSingleRuptureVariationFileInsertion(sess, totalFilesList.get(i));
			if ((i+1)%250==0) System.gc();
			if ((i+1)%50==0) {
				// flush a batch of inserts and release memory
				sess.flush();
				sess.clear();
			}
			/*if ((i+1)%1000==0) {
				sess.getTransaction().commit();
			}*/
		}
	}

	@SuppressWarnings("unused")
	private static void writeFileListToFile(BufferedWriter filelistwriter, File[] safilesList, int i) throws IOException {
		//System.out.println("Filename: " + safilesList[i].getName());
		filelistwriter.write("Filename: " + safilesList[i].getName());
		filelistwriter.newLine();
	}

	private static void prepAndExecuteSingleRuptureVariationFileInsertion(Session sess, File bsaFile) {
		SARuptureFromRuptureVariationFile saRuptureWithSingleRupVar = new SARuptureFromRuptureVariationFile(bsaFile, siteName);
		
		//printEastandNorthComponents(saRupture);
		saRuptureWithSingleRupVar.computeAllGeomAvgComponents();
		//printGeomAvgComponents(saRupture);
		/*System.out.println("Source_ID: " + saRuptureWithSingleRupVar.getSourceID() + 
				", Rupture_ID: " + saRuptureWithSingleRupVar.getRuptureID() + 
				", Rup_Var_ID: " + saRuptureWithSingleRupVar.rupVar.variationNumber);*/

		insertRupture(saRuptureWithSingleRupVar, sess);
	}

	private static void insertRupture(SARuptureFromRuptureVariationFile saRuptureWithSingleRupVar, Session sess) {

		SAPeriods saPeriods = new SAPeriods();



		int outerLoopMax = saRuptureWithSingleRupVar.rupVars.size();
		/*System.out.println("number of rupture variations: " + saRupture.rupVars.size());*/
		for (int rupVarIter=0;rupVarIter<outerLoopMax;rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);

			// Do garbage collection
			//if ((rupVarIter+1)%250==0) System.gc();
			int currentERF_ID = 29;
			int currentSource_ID = saRuptureWithSingleRupVar.getSourceID();
			int currentRupture_ID = saRuptureWithSingleRupVar.getRuptureID();
			int currentRup_Var_Scenario_ID = 1;

			SARuptureVariation currRupVar = saRuptureWithSingleRupVar.rupVar;

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
				paPK.setSite_ID(siteID);				
				paPK.setSource_ID(currentSource_ID);				
				paPK.setRupture_ID(currentRupture_ID);
				paPK.setRup_Var_ID(currRupVar.variationNumber);
				paPK.setRup_Var_Scenario_ID(currentRup_Var_Scenario_ID);
				paPK.setIM_Type(new String("SA_Period_" + saPeriods.getNextValue()));
				pa.setPaPK(paPK);
				pa.setIM_Value(currRupVar.geomAvgComp.periods[periodIter]);
				pa.setUnits("cm per second squared");


				// 3. Save and Commit PeakAmplitude instance

				//sess.beginTransaction();
				
				/*Criteria crit = sess.createCriteria(PeakAmplitude.class);
				crit.setMaxResults(50);
				Example examplePA = Example.create(pa).excludeProperty("IM_Value").excludeProperty("Units");
				crit.add(examplePA);
				List peakamps = crit.list();
				
				if (peakamps.size() == 0)  {
					sess.save(pa);
				}*/
				
				sess.save(pa);
				
				/*if ((rupVarIter+1)%50==0) {
					// flush a batch of inserts and release memory
					sess.flush();
					sess.clear();
				}*/
			}
		}
	}
}
