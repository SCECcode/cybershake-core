package processing;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
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
import org.hibernate.NonUniqueObjectException;
import org.hibernate.SQLQuery;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;
import org.hibernate.exception.ConstraintViolationException;

import commands.CyberLoadamps.Mode;
import util.BSAFileUtil;
import util.NumberHelper;
import util.SwapBytes;
import data.BSAHeader;
import data.DirectionalComponent;
import data.RotDEntry;
import data.RunID;
import data.SAPeriods;
import data.SARuptureFromRuptureVariationFile;
import data.SARuptureVariation;

public class RuptureVariationFileInserter {
	
	private static ArrayList<File> totalFilesList;

	private SessionFactory sessFactory;
	
	private String pathName;
	private String hibernate_cfg_filename = "intensity.cfg.xml";
	private RunID run_ID;

    private ArrayList<Double> desiredPeriods = new ArrayList<Double>();
    private ArrayList<Integer> desiredPeriodsIndices = null;
	//maps period indices to IM Type IDs
	HashMap<Integer, Integer> periodIndexToIDMap = null;
	HashMap<Float, Integer> rd100periodValueToIDMap = null;
	HashMap<Float, Integer> rd50periodValueToIDMap = null;
	
    private ArrayList<Integer> desiredPeriodsIndicesX = null;
	HashMap<Integer, Integer> periodIndexToIDMapX = null;
    private ArrayList<Integer> desiredPeriodsIndicesY = null;
	HashMap<Integer, Integer> periodIndexToIDMapY = null;

	private Mode fileMode;
    private boolean insertGeoMean = false;
    private boolean insertXY = false;
    
    private boolean convertGtoCM = false;
    private final float G_TO_CM_2 = 980.665f;

	public RuptureVariationFileInserter(String newPathName, RunID rid, String serverName, Mode m, String periods, String insertValues, boolean convert) throws IOException {
		pathName = newPathName;
		fileMode = m;
		run_ID = rid;
		if (insertValues.indexOf("gm")!=-1) {
			insertGeoMean = true;
		}
		if (insertValues.indexOf("xy")!=-1) {
			insertXY = true;
		}
		if (insertGeoMean==false && insertXY == false) { //neither option was picked.  Bad.
			System.err.println("Insertion values was " + insertValues + ", but it must be one of gm, xy, or gmxy.");
			System.exit(-1);
		}
		
		if (serverName.equals("intensity")) {
			hibernate_cfg_filename = "intensity.cfg.xml";
		} else if (serverName.equals("surface")) {
			hibernate_cfg_filename = "surface.cfg.xml";
		} else if (serverName.equals("focal")) {
			hibernate_cfg_filename = "focal.cfg.xml";
		} else {
			System.err.println("Server name was " + serverName + ", but it must be one of intensity, surface, or focal.  Exiting.");
			System.exit(-2);
		}
		String[] pieces = periods.split(",");
		for (String p: pieces) {
			desiredPeriods.add(Double.parseDouble(p));
		}
		convertGtoCM = convert;
		initSessionFactory();
	}

	private void initFileList(Mode m) {
		BSAFileUtil.totalFilenameList = new ArrayList<String>();
		BSAFileUtil.totalFileList = new ArrayList<File>();
		File saFile = new File(pathName);
		totalFilesList = BSAFileUtil.createTotalFileList(saFile, m);
	}

	
//	private void retrieveSiteIDFromDB() {
//		initSessionFactory();
//		Session retrieveSiteIDSess = sessFactory.openSession();
//		
//		siteIDQuery = "SELECT CS_Site_ID FROM CyberShake_Sites WHERE CS_Short_Name = '" + siteName + "'";
//		//System.out.println(query);
//		List siteIDList = retrieveSiteIDSess.createSQLQuery(siteIDQuery).addScalar("CS_Site_ID", Hibernate.INTEGER).list();
//		Object siteIDObject = siteIDList.get(0); 
////		siteID = (Integer)siteIDObject;
//		
//		//System.out.println("Site_ID for " + siteName + ": " + siteID);
//		retrieveSiteIDSess.close();
//	}

	private void initSessionFactory() {
		sessFactory = new Configuration().configure(hibernate_cfg_filename).buildSessionFactory();
	}

	public void performInsertions() {
//		retrieveSiteIDFromDB();
		initFileList(fileMode);
		Session sess = sessFactory.openSession();
		
		if (fileMode==Mode.ZIP) {
			for (File f: totalFilesList) {
				System.out.println(f.getName());
			}
			insertRuptureVariationFilesFromZip(sess);
		} else if (fileMode==Mode.HEAD) {
			insertRuptureVariationFilesWithHeader(sess);
		} else if (fileMode==Mode.ROTD) {
			insertRotDFiles(sess);
		} else {
			insertAllRuptureVariationFiles(sess);

		}
		sess.getTransaction().commit();
		sess.close();
	}

	private void insertRotDFiles(Session sess) {
		int counter = 0;
		for (File f: totalFilesList) {
			try {
				/*Each rotd file has
				 * <PSA header for RV1>
				 * <num of rotD entries>
				 * <RotD entry 1>
				 * <RotD entry 2>
				 * ...
				 * <RotD entry 16>
				 * <PSA header for RV2
				 */
				System.out.println("Processing file " + f.getName());
				FileInputStream stream = new FileInputStream(f);
				BSAHeader head = new BSAHeader();
				
				try {
					while (true) {
						int ret = head.parse(stream);
						if (ret==-1) break;
						//Check the site name
						if (!head.siteString.equals(run_ID.getSiteName())) {
							System.err.println("Header in " + f.getName() + " lists site name as " + head.siteString + ", but the site for run ID " + run_ID.getRunID() + " is " + run_ID.getSiteName());
							System.exit(-6);
						}
						ArrayList<RotDEntry> entries = new ArrayList<RotDEntry>();
						//Get the # of periods
						byte[] period_bytes = new byte[4];
						stream.read(period_bytes);
						byte[] out_period_bytes = SwapBytes.swapByteArrayToByteArrayForFloats(period_bytes);
						DataInputStream ds = new DataInputStream(new ByteArrayInputStream(out_period_bytes));
						int num_periods = ds.readInt();
						//16 bytes per RotD entry
						byte[] rotdEntryBytes = new byte[num_periods*16];
						stream.read(rotdEntryBytes);
						byte[] outByteArray = SwapBytes.swapByteArrayToByteArrayForFloats(rotdEntryBytes);
						ds = new DataInputStream(new ByteArrayInputStream(outByteArray));
						for (int i=0; i<num_periods; i++) {
							RotDEntry rde = new RotDEntry();
							rde.populate(ds);
							entries.add(rde);
						}
						ds.close();
						insertRotdRupture(entries, head, sess);
					}
				} catch (IOException ie) {
					ie.printStackTrace();
				}
				
				//Do this more often because there are multiple inserts per file
				if ((counter+1)%25==0) System.gc();
				if ((counter+1)%5==0) {
				// flush a batch of inserts and release memory
					sess.flush();
					sess.clear();
				}
				counter++;
			} catch (IOException ioe) {
				System.err.println("Error reading from file " + pathName);
			} catch (ConstraintViolationException ex) {
				ex.printStackTrace();
				System.err.println("Offending SQL statement was: " + ex.getSQL());
				System.exit(-2);
			}
		}
	}

	
	private void insertRuptureVariationFilesWithHeader(Session sess) {
		int counter = 0;
		//Set this now.  Do it.
		DirectionalComponent.periodsLength = SAPeriods.num_head_periods;
		for (File f: totalFilesList) {
			try {
				/*Each file consists of
				 * <PSA header for RV1>
				 * <RV1, comp 1>
				 * ...
				 * <RV1, comp n>
				 * <PSA header, RV2>
				 * <RV2, comp 1>
				 * ...
				 */
				//Can't use DataInputStream because we have to byte-swap
				System.out.println("Processing file " + f.getName());
				FileInputStream stream = new FileInputStream(f);
				BSAHeader head = new BSAHeader();
				
				try {
					while (true) { //Will exit when we reach the end of the file
						//Read the header
						int ret = head.parse(stream);
						if (ret==-1) break;
//						head.print();
						//Check the site name
						if (!head.siteString.equals(run_ID.getSiteName())) {
							System.err.println("Header in " + f.getName() + " lists site name as " + head.siteString + ", but the site for run ID " + run_ID.getRunID() + " is " + run_ID.getSiteName());
							System.exit(-6);
						}
						//Get the data
						//*4 since they're floats
						byte[] data = new byte[SAPeriods.num_head_periods * head.num_comps * 4];
						ret = stream.read(data);
						if (ret==-1) { //we didn't read anything and this is a weird place to have stopped
							System.err.println("Found header but not data for " + f.getName() + " at offset " + stream.getChannel().position());
							System.exit(-2);
						}
						SARuptureFromRuptureVariationFile saRuptureWithSingleRupVar = new SARuptureFromRuptureVariationFile(data, run_ID.getSiteName(), head);
						insertRupture(saRuptureWithSingleRupVar, sess, true);
					}
				} catch (IOException e) {
					//just means we reached the end of the file
					e.printStackTrace();
				}
				
				stream.close();
				
				//Do this more often because there are multiple inserts per file
				if ((counter+1)%25==0) System.gc();
				if ((counter+1)%5==0) {
				// flush a batch of inserts and release memory
					sess.flush();
					sess.clear();
				}
				counter++;
			} catch (IOException ex) {
				System.err.println("Error reading from file " + pathName);
			} catch (ConstraintViolationException ex) {
				ex.printStackTrace();
				System.err.println("Offending SQL statement was: " + ex.getSQL());
				System.exit(-2);
			}
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
				
				SARuptureFromRuptureVariationFile saRuptureWithSingleRupVar = new SARuptureFromRuptureVariationFile(data, run_ID.getSiteName(), ze);
//				saRuptureWithSingleRupVar.computeAllGeomAvgComponents();
				insertRupture(saRuptureWithSingleRupVar, sess, false);				
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
		} catch (ConstraintViolationException ex) {
			ex.printStackTrace();
			System.err.println("Offending SQL statement was: " + ex.getSQL());
			System.exit(-2);
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
		SARuptureFromRuptureVariationFile saRuptureWithSingleRupVar = new SARuptureFromRuptureVariationFile(bsaFile, run_ID.getSiteName());
		
		//printEastandNorthComponents(saRupture);
		saRuptureWithSingleRupVar.computeAllGeomAvgComponents();
		//printGeomAvgComponents(saRupture);
		/*System.out.println("Source_ID: " + saRuptureWithSingleRupVar.getSourceID() + 
				", Rupture_ID: " + saRuptureWithSingleRupVar.getRuptureID() + 
				", Rup_Var_ID: " + saRuptureWithSingleRupVar.rupVar.variationNumber);*/

		insertRupture(saRuptureWithSingleRupVar, sess, false);
	}

	private void insertRotdRupture(ArrayList<RotDEntry> entries, BSAHeader head, Session sess) {
		Session rotdSession = sessFactory.openSession();
		
		String rd100Prefix = "SELECT IM_Type_ID FROM IM_Types WHERE IM_Type_Measure = 'spectral acceleration' AND IM_Type_Component = 'RotD100' AND ";
		String rd50Prefix = "SELECT IM_Type_ID FROM IM_Types WHERE IM_Type_Measure = 'spectral acceleration' AND IM_Type_Component = 'RotD50' AND ";
		
		//Determine mapping from periods to IM_Type_IDs
		if (rd100periodValueToIDMap==null) {
			rd100periodValueToIDMap = new HashMap<Float, Integer>();
			for (int i=0; i<desiredPeriods.size(); i++) {
				SQLQuery query = rotdSession.createSQLQuery(rd100Prefix + "IM_Type_Value = " + desiredPeriods.get(i)).addScalar("IM_Type_ID", Hibernate.INTEGER);
				int typeID = (Integer)query.list().get(0);
				rd100periodValueToIDMap.put(desiredPeriods.get(i).floatValue(), typeID);
				System.out.println("Adding IM_Type_ID " + typeID + " to list.");
			}
		}
		if (rd50periodValueToIDMap==null) {
			rd50periodValueToIDMap = new HashMap<Float, Integer>();
			for (int i=0; i<desiredPeriods.size(); i++) {
				SQLQuery query = rotdSession.createSQLQuery(rd50Prefix + "IM_Type_Value = " + desiredPeriods.get(i)).addScalar("IM_Type_ID", Hibernate.INTEGER);
				int typeID = (Integer)query.list().get(0);
				rd50periodValueToIDMap.put(desiredPeriods.get(i).floatValue(), typeID);
				System.out.println("Adding IM_Type_ID " + typeID + " to list.");
			}
		}
		
		sess.beginTransaction();
		
		//Look for entries with matching period
		for (RotDEntry e: entries) {
			if (rd100periodValueToIDMap.containsKey(e.period)) {
//				System.out.println("Adding source " + head.source_id + ", rupture " + head.rupture_id + ", rup var " + head.rup_var_id + ", period " + e.period);
				// Initialize PeakAmplitudes class
				PeakAmplitude pa = new PeakAmplitude();
				PeakAmplitudePK paPK = new PeakAmplitudePK();
				// Set values for the PeakAmplitudes Class
				paPK.setRun_ID(run_ID.getRunID());
				paPK.setSource_ID(head.source_id);				
				paPK.setRupture_ID(head.rupture_id);
				paPK.setRup_Var_ID(head.rup_var_id);
				paPK.setIM_Type_ID(rd100periodValueToIDMap.get(e.period));
				pa.setPaPK(paPK);
				float rd100 = e.rotD100;
				if (convertGtoCM) {
					rd100 = e.rotD100 * G_TO_CM_2;
				}
				pa.setIM_Value(rd100);
				try {
					sess.save(pa);
				} catch (NonUniqueObjectException nuoe) {
					//Occurs if there's a duplicate entry in the PSA file, which can happen on rare occasions.  REport and keep going.
					System.err.println("ERROR:  found duplicate entry for run_id " + paPK.getRun_ID() + ", source " + paPK.getSource_ID() + " rupture " + paPK.getRupture_ID() + " rup_var " + paPK.getRup_Var_ID() + " IM_Type " + paPK.getIM_Type_ID() + ".  Skipping.");
				}
			}
			if (rd50periodValueToIDMap.containsKey(e.period)) {
				// Initialize PeakAmplitudes class
				PeakAmplitude pa = new PeakAmplitude();
				PeakAmplitudePK paPK = new PeakAmplitudePK();
				// Set values for the PeakAmplitudes Class
				paPK.setRun_ID(run_ID.getRunID());
				paPK.setSource_ID(head.source_id);				
				paPK.setRupture_ID(head.rupture_id);
				paPK.setRup_Var_ID(head.rup_var_id);
				paPK.setIM_Type_ID(rd50periodValueToIDMap.get(e.period));
				pa.setPaPK(paPK);
				float rd50 = e.rotD50;
				if (convertGtoCM) {
					rd50 = e.rotD50 * G_TO_CM_2;
				}
				pa.setIM_Value(rd50);
				try {
					sess.save(pa);
				} catch (NonUniqueObjectException nuoe) {
					//Occurs if there's a duplicate entry in the PSA file, which can happen on rare occasions.  REport and keep going.
					System.err.println("ERROR:  found duplicate entry for run_id " + paPK.getRun_ID() + ", source " + paPK.getSource_ID() + " rupture " + paPK.getRupture_ID() + " rup_var " + paPK.getRup_Var_ID() + " IM_Type " + paPK.getIM_Type_ID() + ".  Skipping.");
				}
			}
		}
	}

	private void insertRupture(SARuptureFromRuptureVariationFile saRuptureWithSingleRupVar, Session sess, boolean headers) {
		//open session
		Session imTypeIDSess = sessFactory.openSession();
		
		double[] ourPeriods = SAPeriods.values;
		if (headers) {
			ourPeriods = SAPeriods.head_values;
		}
		
		String imTypeIDqueryPrefix = "SELECT IM_Type_ID FROM IM_Types WHERE IM_Type_Measure = 'spectral acceleration' AND IM_Type_Component='geometric mean' AND "; 
		
		//going to only insert 2.0s, 3.0s, 5.0s, 10s periods to save space
		//Going to skip 2.0s since we don't really trust it
		if (insertGeoMean) {
			if (desiredPeriodsIndices==null) {
				desiredPeriodsIndices = new ArrayList<Integer>();
				periodIndexToIDMap = new HashMap<Integer, Integer>();
				for (int i=0; i<ourPeriods.length; i++) {
					for (int j=0; j<desiredPeriods.size(); j++) {
						if (Math.abs(ourPeriods[i]-desiredPeriods.get(j))<0.0001) {
							desiredPeriodsIndices.add(i);
							SQLQuery query = imTypeIDSess.createSQLQuery(imTypeIDqueryPrefix + "IM_Type_Value = " + ourPeriods[i]).addScalar("IM_Type_ID", Hibernate.INTEGER);
							int typeID = (Integer)query.list().get(0);
							periodIndexToIDMap.put(i, typeID);
						}
					}
				}
			}
		}
		
		if (insertXY) {
			String imTypeIDqueryPrefixX = "SELECT IM_Type_ID FROM IM_Types WHERE IM_Type_Measure = 'spectral acceleration' AND IM_Type_Component = 'X component' AND ";
			String imTypeIDqueryPrefixY = "SELECT IM_Type_ID FROM IM_Types WHERE IM_Type_Measure = 'spectral acceleration' AND IM_Type_Component = 'Y component' AND ";
			
			if (desiredPeriodsIndices==null) {
				desiredPeriodsIndices = new ArrayList<Integer>();
				periodIndexToIDMapX = new HashMap<Integer, Integer>();
				periodIndexToIDMapY = new HashMap<Integer, Integer>();
				for (int i=0; i<ourPeriods.length; i++) {
					for (int j=0; j<desiredPeriods.size(); j++) {
						if (Math.abs(ourPeriods[i]-desiredPeriods.get(j))<0.0001) {
							desiredPeriodsIndices.add(i);
							int typeID = (Integer)imTypeIDSess.createSQLQuery(imTypeIDqueryPrefixX + "IM_Type_Value = " + ourPeriods[i]).addScalar("IM_Type_ID", Hibernate.INTEGER).list().get(0);
							periodIndexToIDMapX.put(i, typeID);
							typeID = (Integer)imTypeIDSess.createSQLQuery(imTypeIDqueryPrefixY + "IM_Type_Value = " + ourPeriods[i]).addScalar("IM_Type_ID", Hibernate.INTEGER).list().get(0);
							periodIndexToIDMapY.put(i, typeID);
						}
					}
				}
			} else { //desiredPeriodIndices is already populated
				periodIndexToIDMapX = new HashMap<Integer, Integer>();
				periodIndexToIDMapY = new HashMap<Integer, Integer>();
				for (int period: desiredPeriodsIndices) {
					int typeID = (Integer)imTypeIDSess.createSQLQuery(imTypeIDqueryPrefixX + "IM_Type_Value = " + ourPeriods[period]).addScalar("IM_Type_ID", Hibernate.INTEGER).list().get(0);
					periodIndexToIDMapX.put(period, typeID);
					typeID = (Integer)imTypeIDSess.createSQLQuery(imTypeIDqueryPrefixY + "IM_Type_Value = " + ourPeriods[period]).addScalar("IM_Type_ID", Hibernate.INTEGER).list().get(0);
					periodIndexToIDMapY.put(period, typeID);
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

				//for error checking
				float mag = -1.0f;
				float rupture_dist = -1.0f;
				if (insertGeoMean) {
					// Initialize PeakAmplitudes class
					PeakAmplitude pa = new PeakAmplitude();
					PeakAmplitudePK paPK = new PeakAmplitudePK();
					// Set values for the PeakAmplitudes Class
					paPK.setRun_ID(run_ID.getRunID());
					paPK.setSource_ID(currentSource_ID);				
					paPK.setRupture_ID(currentRupture_ID);
					paPK.setRup_Var_ID(currRupVar.variationNumber);
					paPK.setIM_Type_ID(periodIndexToIDMap.get(periodIter));
					pa.setPaPK(paPK);
					double psaValue = currRupVar.geomAvgComp.periods[periodIter];
					if (psaValue>8400 || psaValue<0.01) {
						//Tiny PSA values are ok if the event is small and far away
						System.out.println("Found psaValue " + psaValue + " for source " + currentSource_ID + " rupture " + currentRupture_ID + " rup var " + currRupVar.variationNumber);
						if (mag==-1.0) {
							//Then we need to do the query for this rupture
							Session checkMagSession = sessFactory.openSession();
							String query = "select R.Mag, SR.Site_Rupture_Dist from Ruptures R, CyberShake_Site_Ruptures SR " +
							"where SR.CS_Site_ID=" + run_ID.getSiteID() + " and SR.ERF_ID=" + run_ID.getErfID() +
							" and SR.Source_ID=" + currentSource_ID + " and SR.Rupture_ID=" + currentRupture_ID +
							" and R.ERF_ID=SR.ERF_ID and R.Source_ID=SR.Source_ID and R.Rupture_ID=SR.Rupture_ID";
							System.out.println(query);
							System.out.flush();
							SQLQuery sq = checkMagSession.createSQLQuery(query);
							sq.addScalar("R.Mag", Hibernate.FLOAT);
							sq.addScalar("SR.Site_Rupture_Dist", Hibernate.FLOAT);
							List<Object[]> results = sq.list();
							mag = (Float)(results.get(0)[0]);
							//Rupture distances aren't available for all ruptures
							if (results.get(0).length>1) {
								rupture_dist = (Float)(results.get(0)[1]);
							} else {
								rupture_dist = -1;
							}
							checkMagSession.close();
						}
						if (mag>6.8 || rupture_dist<300.0 || psaValue<0.003) {
							System.err.println("Found value " + psaValue + " for source " + currentSource_ID + ", " + currentRupture_ID + ", " + currRupVar.variationNumber + ", period index " + periodIter + ", period value " + ourPeriods[periodIter]);
							System.err.println("Mag=" + mag + " rupture_dist=" + rupture_dist);
							throw new IllegalArgumentException();
						} else {
							System.out.println("Found value " + psaValue + " for source " + currentSource_ID + ", " + currentRupture_ID + ", " + currRupVar.variationNumber + ", period index " + periodIter + ", period value " + ourPeriods[periodIter]);
							System.out.println("Since mag=" + mag + " and source rupture dist=" + rupture_dist + ", permitting insertion.");
						}
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
				
					try {
						sess.save(pa);
					} catch (NonUniqueObjectException nuoe) {
						//Occurs if there's a duplicate entry in the PSA file, which can happen on rare occasions.  REport and keep going.
						System.err.println("ERROR:  found duplicate entry for run_id " + paPK.getRun_ID() + ", source " + paPK.getSource_ID() + " rupture " + paPK.getRupture_ID() + " rup_var " + paPK.getRup_Var_ID() + " IM_Type " + paPK.getIM_Type_ID() + ".  Skipping.");
					}
				}
				
				if (insertXY) {
					//do X
					PeakAmplitude pa = new PeakAmplitude();
					PeakAmplitudePK paPK = new PeakAmplitudePK();
					paPK.setRun_ID(run_ID.getRunID());
					paPK.setSource_ID(currentSource_ID);				
					paPK.setRupture_ID(currentRupture_ID);
					paPK.setRup_Var_ID(currRupVar.variationNumber);
					paPK.setIM_Type_ID(periodIndexToIDMapX.get(periodIter));
					pa.setPaPK(paPK);
					double psaValue = currRupVar.eastComp.periods[periodIter];
					if (psaValue>8400 || psaValue<0.01) {
						System.err.println("Found value " + psaValue + " for source " + currentSource_ID + ", " + currentRupture_ID + ", " + currRupVar.variationNumber + ", period index " + periodIter + ", period value " + ourPeriods[periodIter]);
						throw new IllegalArgumentException();
					}
//					System.out.println("Inserting value " + psaValue + " for source " + currentSource_ID + ", " + currentRupture_ID + ", " + currRupVar.variationNumber + ", period index " + periodIter + ", period value " + saPeriods.values[periodIter]);
					pa.setIM_Value(psaValue);
					
					sess.save(pa);
					
					//do Y
					pa = new PeakAmplitude();
					paPK = new PeakAmplitudePK();
					paPK.setRun_ID(run_ID.getRunID());
					paPK.setSource_ID(currentSource_ID);				
					paPK.setRupture_ID(currentRupture_ID);
					paPK.setRup_Var_ID(currRupVar.variationNumber);
					paPK.setIM_Type_ID(periodIndexToIDMapY.get(periodIter));
					pa.setPaPK(paPK);
					psaValue = currRupVar.northComp.periods[periodIter];
					if (psaValue>8400 || psaValue<0.01) {
						System.err.println("Found value " + psaValue + " for source " + currentSource_ID + ", " + currentRupture_ID + ", " + currRupVar.variationNumber + ", period index " + periodIter + ", period value " + ourPeriods[periodIter]);
						throw new IllegalArgumentException();
					}
//					System.out.println("Inserting value " + psaValue + " for source " + currentSource_ID + ", " + currentRupture_ID + ", " + currRupVar.variationNumber + ", period index " + periodIter + ", period value " + saPeriods.values[periodIter]);
					pa.setIM_Value(psaValue);
					
					sess.save(pa);
				}
				
				/*if ((rupVarIter+1)%50==0) {
					// flush a batch of inserts and release memory
					sess.flush();
					sess.clear();
				}*/
			}
		}
	}

	public String getPathName() {
		return pathName;
	}

	public void setPathName(String pathName) {
		this.pathName = pathName;
	}
}
