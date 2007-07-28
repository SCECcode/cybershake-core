package test.incompletetest;

import static org.junit.Assert.assertEquals;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import mapping.PeakAmplitude;
import mapping.PeakAmplitudePK;

import org.hibernate.Criteria;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;
import org.hibernate.criterion.Example;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import util.BSAFileUtil;

import data.SAPeriods;
import data.SARuptureFromRuptureVariationFile;

import junit.framework.JUnit4TestAdapter;

public class CheckMultipleSAFilesOnDB {
	
	private static SessionFactory sf;
	private static Session sess;
	private static PeakAmplitude expectedPA;
	private static SAPeriods saperiods;
	private static SARuptureFromRuptureVariationFile saRupture;
	private static ArrayList<File> totalFileList;
	
	@BeforeClass public static void runBeforeAllTests() {
		sf = new Configuration().configure("intensity.cfg.xml").buildSessionFactory();
		sess = sf.openSession();
		totalFileList = BSAFileUtil.createTotalFileList(new File("safiles/loadampstest"));
		saperiods = new SAPeriods();
	}
	
	@Test public void checkEqualityForAllLoadampstestFiles() throws IOException {
		
		
		System.out.println("totalFileList.size(): " + totalFileList.size());
		for (int fileListIndex=0; fileListIndex<totalFileList.size(); fileListIndex++) {
			saRupture = new SARuptureFromRuptureVariationFile(totalFileList.get(fileListIndex),"USC");
			
			
			int innerLoopMax = saRupture.rupVar.geomAvgComp.periods.length;
			System.out.println("innerLoopMax: " + innerLoopMax);
			
			for (int periodIndex=0; periodIndex<innerLoopMax; periodIndex++) {
				expectedPA = new PeakAmplitude();
				PeakAmplitudePK paPK = new PeakAmplitudePK();
				// Set values for the PeakAmplitudes Class
				paPK.setERF_ID(29);
				paPK.setSite_ID(18);
				paPK.setSource_ID(saRupture.getSourceID());				
				paPK.setRupture_ID(saRupture.getRuptureID());
				paPK.setRup_Var_ID(saRupture.rupVar.variationNumber);
				paPK.setRup_Var_Scenario_ID(1);
				paPK.setIM_Type(new String("SA_Period_" + saperiods.getNextValue()));
				expectedPA.setPaPK(paPK);
				expectedPA.setIM_Value(saRupture.rupVar.geomAvgComp.periods[periodIndex]);
				expectedPA.setUnits("meters per second squared");

				Criteria crit = sess.createCriteria(PeakAmplitude.class);
				crit.setMaxResults(50);
				Example examplePA = Example.create(expectedPA).excludeProperty("IM_Value").excludeProperty("Units");
				crit.add(examplePA);
				List peakamps = crit.list();
				
				System.out.println("o");
				
				//for (int peakampsIndex=0; peakampsIndex< peakamps.size(); peakampsIndex++ ) {
				PeakAmplitude resultsPA = (PeakAmplitude)peakamps.get(0);
				System.out.println(resultsPA.toString());
				
				/*BufferedWriter out = new BufferedWriter(new FileWriter("test/incompletetest/check_multiple_safiles_ondb.txt",true));
				out.append("-- results from example --");
				out.newLine();

				//for (int peakampsIndex=0; peakampsIndex< peakamps.size(); peakampsIndex++ ) {
				PeakAmplitude resultsPA = (PeakAmplitude)peakamps.get(0);
				out.append(resultsPA.toString());
				out.newLine();
				out.close();*/
				
			}

		}
	}
	
	/*@Test public void testFileReadAndObjectInstantiationSpeed() throws IOException {
		//System.out.println("totalFileList.size(): " + totalFileList.size());
		for (int fileListIndex=0; fileListIndex<totalFileList.size(); fileListIndex++) {
			saRupture = new SARuptureFromRuptureVariationFile(totalFileList.get(fileListIndex),"USC");
			
			
			int innerLoopMax = saRupture.rupVar.geomAvgComp.periods.length;
			//System.out.println("innerLoopMax: " + innerLoopMax);
			
			for (int periodIndex=0; periodIndex<innerLoopMax; periodIndex++) {
				expectedPA = new PeakAmplitude();
				PeakAmplitudePK paPK = new PeakAmplitudePK();
				// Set values for the PeakAmplitudes Class
				paPK.setERF_ID(29);
				paPK.setSite_ID(18);
				paPK.setSource_ID(saRupture.getSourceID());				
				paPK.setRupture_ID(saRupture.getRuptureID());
				paPK.setRup_Var_ID(saRupture.rupVar.variationNumber);
				paPK.setRup_Var_Scenario_ID(1);
				paPK.setIM_Type(new String("SA_Period_" + saperiods.getNextValue()));
				expectedPA.setPaPK(paPK);
				expectedPA.setIM_Value(saRupture.rupVar.geomAvgComp.periods[periodIndex]);
				expectedPA.setUnits("meters per second squared");

				BufferedWriter out = new BufferedWriter(new FileWriter("test/incompletetest/check_multiple_safiles_ondb.nodbquery.txt",true));
				out.append(expectedPA.toString());
				out.newLine();
				out.close();
				
			}

		}
	}*/
	
	@AfterClass public static void runAfterAllTests() {
		sess.close();
	}
	
	
	public static junit.framework.Test suite() {
		return new JUnit4TestAdapter (CheckMultipleSAFilesOnDB.class);
	}
}
