package processing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;

import mapping.PeakAmplitude;
import mapping.PeakAmplitudePK;

import org.hibernate.Hibernate;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;

import util.BSAFileUtil;
import util.NumberHelper;
import data.DirectionalComponent;
import data.SAPeriods;
import data.SARuptureFromRuptureVariationFile;
import data.SARuptureVariation;

public class RuptureVariationFileInserter {
	
	public static int siteID = 18;
	private static ArrayList<File> totalFilesList;
	private String siteIDQuery;
	private SessionFactory sessFactory;
	
	private String siteName;
	private String pathName;
	private String hibernate_cfg_filename = "intensity.cfg.xml";
	private int currentSGT_Variation_ID;
	//added by SC
    private int currentRup_Var_Scenario_ID;
    private int currentERF_ID;
    private double[] desiredPeriods = {2.0, 3.00003, 5.0, 10.0};
    private ArrayList<Integer> desiredPeriodsIndices = null;
	//maps period indices to IM Type IDs
	HashMap<Integer, Integer> periodIndexToIDMap = null;

    private boolean zipOption;
    
	public RuptureVariationFileInserter(String newPathName, String newSiteName, String sgtVariationID, String serverName, String rupVarID, String erfID, boolean zipOpt) throws IOException {
		siteName = newSiteName;
		pathName = newPathName;
		zipOption = zipOpt;
		
		if (serverName.equals("intensity")) {
			hibernate_cfg_filename = "intensity.cfg.xml";
		}
		else if (serverName.equals("surface")) {
			hibernate_cfg_filename = "surface.cfg.xml";
		} else if (serverName.equals("focal")) {
			hibernate_cfg_filename = "focal.cfg.xml";
		}
		
		try {
		    currentSGT_Variation_ID = Integer.parseInt(sgtVariationID);
		    //added by SC
		    currentRup_Var_Scenario_ID = Integer.parseInt(rupVarID);
		    currentERF_ID = Integer.parseInt(erfID);
		} catch (NumberFormatException nfe) {
            System.out.println("SGT Variation ID and Rupture Variation ID must be positive integers.");
            System.exit(1);
        }
	}

	private void initFileList(boolean zipOpt) {
		BSAFileUtil.totalFilenameList = new ArrayList<String>();
		BSAFileUtil.totalFileList = new ArrayList<File>();
		File saFile = new File(pathName);
		totalFilesList = BSAFileUtil.createTotalFileList(saFile, zipOpt);
	}

	private void retrieveSiteIDFromDB() {
		initSessionFactory();
		Session retrieveSiteIDSess = sessFactory.openSession();
		
		siteIDQuery = "SELECT CS_Site_ID FROM CyberShake_Sites WHERE CS_Short_Name = '" + siteName + "'";
		//System.out.println(query);
		List siteIDList = retrieveSiteIDSess.createSQLQuery(siteIDQuery).addScalar("CS_Site_ID", Hibernate.INTEGER).list();
		Object siteIDObject = siteIDList.get(0); 
		siteID = (Integer)siteIDObject;
		
		//System.out.println("Site_ID for " + siteName + ": " + siteID);
		retrieveSiteIDSess.close();
	}

	private void initSessionFactory() {
		sessFactory = new Configuration().configure(hibernate_cfg_filename).buildSessionFactory();
	}

	public void performInsertions() {
		retrieveSiteIDFromDB();
		initFileList(zipOption);
		
		if (zipOption) {
			for (File f: totalFilesList) {
				System.out.println(f.getName());
			}
			Session sess = sessFactory.openSession();
			insertRuptureVariationFilesFromZip(sess);
			sess.getTransaction().commit();
			sess.close();
		} else {
			Session sess = sessFactory.openSession();
			insertAllRuptureVariationFiles(sess);
			sess.getTransaction().commit();
			sess.close();
		}
	}

	private void insertRuptureVariationFilesFromZip(Session sess) {
		for (File zf: totalFilesList) {
		try {
			System.out.println("Entering zip file " + zf.getName());
			ZipFile saZipFile = new ZipFile(zf);
			Enumeration<? extends ZipEntry> e = saZipFile.entries();
			long size;
			byte[] data;
			int counter = 0;
			ZipEntry ze;
			InputStream is;
			while (e.hasMoreElements()) {
				counter++;
				ze = e.nextElement();
				is = saZipFile.getInputStream(ze);
				size = ze.getSize();
				data = new byte[(int)size];
				if (is.read(data, 0, (int)size)!=size) {
					System.err.println("Error reading " + size + " bytes of zip entry " + ze.getName());
					System.exit(3);
				}
				
				SARuptureFromRuptureVariationFile saRuptureWithSingleRupVar = new SARuptureFromRuptureVariationFile(data, siteName, ze);
//				saRuptureWithSingleRupVar.computeAllGeomAvgComponents();
				insertRupture(saRuptureWithSingleRupVar, sess);				
				if (counter%250==0) System.gc();
				if (counter%100==0) {
					sess.flush();
					sess.clear();
				}
//				for (Float period : saRuptureWithSingleRupVar.rupVar.geomAvgComp.periods) {
//					System.out.println(period);
//				}
				is.close();
			}
		} catch (IOException ex) {
			System.err.println("Error reading from zip file " + pathName);
		}
		}
	}
	
	private void insertAllRuptureVariationFiles(Session sess) {
			for (int i=0; i<totalFilesList.size(); i++) {

			//writeFileListToFile(filelistwriter, safilesList, i);

			//System.out.println("Filename: " + safilesList[i].getName());
				try {
					prepAndExecuteSingleRuptureVariationFileInsertion(sess, totalFilesList.get(i));
					if ((i+1)%250==0) System.gc();
					if ((i+1)%50==0) {
					// flush a batch of inserts and release memory
						sess.flush();
						sess.clear();
					}
				} catch (Exception e) {
					System.out.println("Exception in file " + totalFilesList.get(i).getAbsolutePath());
					e.printStackTrace();
					System.exit(-1);
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

	private void prepAndExecuteSingleRuptureVariationFileInsertion(Session sess, File bsaFile) {
		SARuptureFromRuptureVariationFile saRuptureWithSingleRupVar = new SARuptureFromRuptureVariationFile(bsaFile, siteName);
		
		//printEastandNorthComponents(saRupture);
		saRuptureWithSingleRupVar.computeAllGeomAvgComponents();
		//printGeomAvgComponents(saRupture);
		/*System.out.println("Source_ID: " + saRuptureWithSingleRupVar.getSourceID() + 
				", Rupture_ID: " + saRuptureWithSingleRupVar.getRuptureID() + 
				", Rup_Var_ID: " + saRuptureWithSingleRupVar.rupVar.variationNumber);*/

		insertRupture(saRuptureWithSingleRupVar, sess);
	}

	private void insertRupture(SARuptureFromRuptureVariationFile saRuptureWithSingleRupVar, Session sess) {

		SAPeriods saPeriods = new SAPeriods();

		//open session
		Session imTypeIDSess = sessFactory.openSession();
		
		String imTypeIDqueryPrefix = "SELECT IM_Type_ID FROM IM_Types WHERE IM_Type_Measure = 'spectral acceleration' AND "; 
		
		//going to only insert 2.0s, 3.0s, 5.0s, 10s periods to save space
		if (desiredPeriodsIndices==null) {
			desiredPeriodsIndices = new ArrayList<Integer>();
			periodIndexToIDMap = new HashMap<Integer, Integer>();
			for (int i=0; i<saPeriods.values.length; i++) {
				for (int j=0; j<desiredPeriods.length; j++) {
					if (Math.abs(saPeriods.values[i]-desiredPeriods[j])<0.0001) {
						desiredPeriodsIndices.add(i);
						int typeID = (Integer)imTypeIDSess.createSQLQuery(imTypeIDqueryPrefix + "IM_Type_Value = " + desiredPeriods[j]).addScalar("IM_Type_ID", Hibernate.INTEGER).list().get(0);
						periodIndexToIDMap.put(i, typeID);
					}
				}
			}
		}

		imTypeIDSess.close();

		int outerLoopMax = saRuptureWithSingleRupVar.rupVars.size();
		/*System.out.println("number of rupture variations: " + saRupture.rupVars.size());*/
		for (int rupVarIter=0;rupVarIter<outerLoopMax;rupVarIter++) {
			//System.out.println("rupVarIter: " + rupVarIter);

			// Do garbage collection
			//if ((rupVarIter+1)%250==0) System.gc();
			int currentSource_ID = saRuptureWithSingleRupVar.getSourceID();
			int currentRupture_ID = saRuptureWithSingleRupVar.getRuptureID();

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
			
//			int innerLoopMax = currRupVar.geomAvgComp.periods.length;
			for (int periodIter: desiredPeriodsIndices) {
//			for (int periodIter=0;periodIter<innerLoopMax;periodIter++) {
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
				paPK.setSGT_Variation_ID(currentSGT_Variation_ID);
				paPK.setIM_Type_ID(periodIndexToIDMap.get(periodIter));
				pa.setPaPK(paPK);
				double psaValue = currRupVar.geomAvgComp.periods[periodIter];
				if (psaValue>8400 || psaValue<0.01) {
					System.err.println("Found value " + psaValue + " for source " + currentSource_ID + ", " + currentRupture_ID + ", " + currRupVar.variationNumber + ", period index " + periodIter + ", period value " + saPeriods.values[periodIter]);
					throw new IllegalArgumentException();
				}
//				System.out.println("Inserting value " + psaValue + " for source " + currentSource_ID + ", " + currentRupture_ID + ", " + currRupVar.variationNumber + ", period index " + periodIter + ", period value " + saPeriods.values[periodIter]);
				pa.setIM_Value(psaValue);

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

	public int getCurrentSGT_Variation_ID() {
		return currentSGT_Variation_ID;
	}

	public void setCurrentSGT_Variation_ID(int currentSGT_Variation_ID) {
		this.currentSGT_Variation_ID = currentSGT_Variation_ID;
	}

	public  String getSiteName() {
		return siteName;
	}

	public void setSiteName(String siteName) {
		this.siteName = siteName;
	}

	public String getPathName() {
		return pathName;
	}

	public void setPathName(String pathName) {
		this.pathName = pathName;
	}
}
