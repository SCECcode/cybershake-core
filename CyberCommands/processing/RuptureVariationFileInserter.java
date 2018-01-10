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
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;
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
import org.hibernate.StatelessSession;
import org.hibernate.cfg.Configuration;
import org.hibernate.exception.ConstraintViolationException;

import commands.CyberLoadamps.Mode;
import util.BSAFileUtil;
import util.NumberHelper;
import util.SwapBytes;
import data.BSAHeader;
import data.DirectionalComponent;
import data.DurationEntry;
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

	private HashMap<String, Integer> durationToTypeID = null;
	
	private Mode fileMode;
    private boolean insertGeoMean = false;
    private boolean insertXY = false;
    
    private boolean convertGtoCM = false;
    private final float G_TO_CM_2 = 980.665f;
    
    private boolean forceInsert = false;
    
    private Connection conn = null;

	public RuptureVariationFileInserter(String newPathName, RunID rid, String serverName, Mode m, String periods, String insertValues, boolean convert, boolean forceInsert) throws IOException {
		pathName = newPathName;
		fileMode = m;
		run_ID = rid;
		this.forceInsert = forceInsert;
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
		} else if (serverName.equals("csep-x")) {
			hibernate_cfg_filename = "csep-x.cfg.xml";
		} else if (serverName.equals("moment")) {
			hibernate_cfg_filename = "moment.cfg.xml";
		} else {
			System.err.println("Server name was " + serverName + ", but it must be one of intensity, surface, focal, or csep-x.  Exiting.");
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
		} else if (fileMode==Mode.DURATION) {
			insertDurationFiles(sess);
		} else {
			insertAllRuptureVariationFiles(sess);

		}		
		sess.getTransaction().commit();
		sess.close();
	}

	private void insertDurationFiles(Session sess) {
		int counter = 0;
		for (File f: totalFilesList) {
			try {
				/*Each duration file has
				 * <Header for RV1>
				 * <number of duration entries per component
				 * <X comp, duration entry 1>
				 * <X comp, duration entry 2>
				 * ...
				 * <X comp, duration entry N>
				 * <Y comp, duration entry 1>
				 * ...
				 * <Y comp, duration entry N>
				 * <Header for RV2>
				 * ...
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
						ArrayList<DurationEntry> entries = new ArrayList<DurationEntry>();
						//Get the # of durations per component
						byte[] num_durations_bytes = new byte[4];
						stream.read(num_durations_bytes);
						byte[] out_durations_bytes = SwapBytes.swapByteArrayToByteArrayForFloats(num_durations_bytes);
						DataInputStream ds = new DataInputStream(new ByteArrayInputStream(out_durations_bytes));
						int num_durations = ds.readInt();
						//16 bytes per duration entry (type, type_value, component, value)
						//Times 2 because 2 components
						byte[] durationEntryBytes = new byte[num_durations*2*16];
						stream.read(durationEntryBytes);
						byte[] outByteArray = SwapBytes.swapByteArrayToByteArrayForFloats(durationEntryBytes);
						ds = new DataInputStream(new ByteArrayInputStream(outByteArray));
						for (int i=0; i<2*num_durations; i++) {
							DurationEntry de = new DurationEntry();
							de.populate(ds);
							entries.add(de);
						}
						ds.close();
						insertDurationRupture(entries, head, sess);
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
				ioe.printStackTrace();
			} catch (ConstraintViolationException ex) {
				ex.printStackTrace();
				System.err.println("Offending SQL statement was: " + ex.getSQL());
				System.exit(-2);
			}
		}
	}
		
	
	private void insertRotDFiles(Session sess) {
		//Track time
		long start, end;
		long fileReading = 0;
		long insertionSetup = 0;
		long dbCommit = 0;
		long garbageCollection = 0;
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
				start = System.currentTimeMillis();
				System.out.println("Processing file " + f.getName());
				FileInputStream stream = new FileInputStream(f);
				BSAHeader head = new BSAHeader();
				end = System.currentTimeMillis();
				fileReading += (end-start);
				
				try {
					while (true) {
						start = System.currentTimeMillis();
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
						end = System.currentTimeMillis();
						fileReading += (end-start);
						start = System.currentTimeMillis();
//						insertRotdRuptureSQL(entries, head);
						insertRotdRupture(entries, head, sess);
						end = System.currentTimeMillis();
						insertionSetup += (end-start);
					}
				} catch (IOException ie) {
					ie.printStackTrace();
				}
				
				//Do this more often because there are multiple inserts per file
				if ((counter+1)%25==0) {
					start = System.currentTimeMillis();
					System.gc();
					end = System.currentTimeMillis();
					garbageCollection += (end-start);
				}
				if ((counter+1)%25==0) {
				// flush a batch of inserts and release memory
					start = System.currentTimeMillis();
					sess.flush();
					sess.clear();
					end = System.currentTimeMillis();
					dbCommit += (end-start);
				}
				counter++;
			} catch (IOException ioe) {
				System.err.println("Error reading from file " + pathName);
				ioe.printStackTrace();
			} catch (ConstraintViolationException ex) {
				ex.printStackTrace();
				System.err.println("Offending SQL statement was: " + ex.getSQL());
				System.exit(-2);
			}
		}
		System.out.println("File reading (sec): " + (fileReading/1000));
		System.out.println("Insertion setup (sec): " + (insertionSetup/1000));
		System.out.println("Garbage collection (sec): " + (garbageCollection/1000));
		System.out.println("DB Commit (sec): " + (dbCommit/1000));
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
//				System.out.println("Processing file " + f.getName());
//				System.out.println("Prepping for loop.");
				FileInputStream stream = new FileInputStream(f);
				BSAHeader head = new BSAHeader();
				
				boolean isFileEmpty = true;
				
				try {
					while (true) { //Will exit when we reach the end of the file
						//Read the header
//						System.out.println("Reading header.");
						int ret = head.parse(stream);
						if (ret==-1) break;
						//Unset if we read something from the file
						isFileEmpty = false;
//						head.print();
						//Check the site name
						if (!head.siteString.equals(run_ID.getSiteName())) {
							System.err.println("Header in " + f.getName() + " lists site name as " + head.siteString + ", but the site for run ID " + run_ID.getRunID() + " is " + run_ID.getSiteName());
							System.exit(-6);
						}
						//Get the data
						//*4 since they're floats
//						System.out.println("Reading " + SAPeriods.num_head_periods * head.num_comps * 4 +  " bytes.");
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
				
				//Make sure there was something in the file
				if (isFileEmpty) {
					System.err.println("File " + f.getName() + " was empty, aborting.");
					System.exit(-3);
				}
				
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

	private void insertDurationRupture(ArrayList<DurationEntry> entries, BSAHeader head, Session sess) {
		Session durationSession = sessFactory.openSession();
		
		String durationPrefix = "SELECT IM_Type_ID FROM IM_Types WHERE IM_Type_Measure = ";
		
		if (durationToTypeID==null) {
			durationToTypeID = new HashMap<String, Integer>();
			//Determine mapping to IM_Type_IDs
			//For now, we will insert for both X and Y components, for both acceleration and velocity, 5-75% and 5-95% cutoffs
			String[] type = {"acceleration", "velocity"};
			String[] measure = {"5% to 75%", "5% to 95%"};
			String[] component = {"X", "Y"};
			for (String c: component) {
				for (String t: type) {
					for (String m: measure) {
						String typeMeasureString = t + " significant duration, " + m;
						String query = durationPrefix + "'" + typeMeasureString + "' AND IM_Type_Component='" + c + "'";
						SQLQuery q = durationSession.createSQLQuery(query);
						int typeID = (Integer)(q.list().get(0));
						String keyString = DurationEntry.getKeyString(typeMeasureString, c);
						durationToTypeID.put(keyString, typeID);
					}
				}
			}
		}
		
		sess.beginTransaction();
		
		for (DurationEntry e: entries) {
			String keyString = "" + e.type;
			if (e.type_value!=-1) {
				keyString += "_" + e.type_value;
			}
			keyString += "_" + e.component;
			
			if (durationToTypeID.containsKey(keyString)) {
				PeakAmplitude pa = new PeakAmplitude();
				PeakAmplitudePK paPK = new PeakAmplitudePK();
				paPK.setRun_ID(run_ID.getRunID());
				paPK.setSource_ID(head.source_id);				
				paPK.setRupture_ID(head.rupture_id);
				paPK.setRup_Var_ID(head.rup_var_id);
				paPK.setIM_Type_ID(durationToTypeID.get(keyString));
				pa.setPaPK(paPK);
				pa.setIM_Value(e.value);
				try {
					sess.save(pa);
				} catch (NonUniqueObjectException nuoe) {
					//Occurs if there's a duplicate entry in the PSA file, which can happen on rare occasions.  Because of the Study 15.4 issues, abort.
					System.err.println("ERROR:  found duplicate entry in file for run_id " + paPK.getRun_ID() + ", source " + paPK.getSource_ID() + " rupture " + paPK.getRupture_ID() + " rup_var " + paPK.getRup_Var_ID() + " IM_Type " + paPK.getIM_Type_ID() + ".  Aborting.");
					System.exit(2);
				}
			}
		}
	}
	
	//Insert rupture data using SQL queries directly, not hibernate
	private void insertRotdRuptureSQL(ArrayList<RotDEntry> entries, BSAHeader head) {
		String rd100Prefix = "SELECT IM_Type_ID FROM IM_Types WHERE IM_Type_Measure = 'spectral acceleration' AND IM_Type_Component = 'RotD100' AND ";
		String rd50Prefix = "SELECT IM_Type_ID FROM IM_Types WHERE IM_Type_Measure = 'spectral acceleration' AND IM_Type_Component = 'RotD50' AND ";
		
		try {
			if (conn==null) {
				Properties props = new Properties();
				props.put("user", "cybershk");
				props.put("password", "re@lStil1");
				conn = DriverManager.getConnection("jdbc:mysql://focal.usc.edu:3306/CyberShake", props);
			}
			Statement stmt = conn.createStatement();

			//Determine mapping from periods to IM_Type_IDs
			if (rd100periodValueToIDMap==null) {
				rd100periodValueToIDMap = new HashMap<Float, Integer>();
				for (int i=0; i<desiredPeriods.size(); i++) {
					String query = rd100Prefix + "IM_Type_Value = " + desiredPeriods.get(i);
					System.out.println(query);
					ResultSet rs = stmt.executeQuery(query);
					try {
						rs.first();
					 	if (rs.getRow()==0) {
				      	    System.err.println("No IM_Type_ID found for RotD100, value " + desiredPeriods.get(i) + ".");
				      	    System.exit(1);
				      	}
					} catch (SQLException e) {
						e.printStackTrace();
						System.exit(2);
					}
					int typeID = rs.getInt(1);
					rd100periodValueToIDMap.put(desiredPeriods.get(i).floatValue(), typeID);
					System.out.println("Adding IM_Type_ID " + typeID + " to list.");
				}
			}
			if (rd50periodValueToIDMap==null) {
				rd50periodValueToIDMap = new HashMap<Float, Integer>();
				for (int i=0; i<desiredPeriods.size(); i++) {
					String query = rd50Prefix + "IM_Type_Value = " + desiredPeriods.get(i);
					ResultSet rs = stmt.executeQuery(query);
					try {
						rs.first();
					 	if (rs.getRow()==0) {
				      	    System.err.println("No IM_Type_ID found for RotD50, value " + desiredPeriods.get(i) + ".");
				      	    System.exit(1);
				      	}
					} catch (SQLException e) {
						e.printStackTrace();
						System.exit(2);
					}
					int typeID = rs.getInt(1);
					rd50periodValueToIDMap.put(desiredPeriods.get(i).floatValue(), typeID);
					System.out.println("Adding IM_Type_ID " + typeID + " to list.");
				}
			}
			String insertString = "insert into PeakAmplitudes " + 
					"(Run_ID, Source_ID, Rupture_ID, Rup_Var_ID, IM_Type_ID, IM_Value) " +
					"values (" + run_ID.getRunID() + ", ?, ?, ?, ?, ?)";
			PreparedStatement ps = conn.prepareStatement(insertString);
			int counter = 0;
			for (RotDEntry e: entries) {
				if (rd100periodValueToIDMap.containsKey(e.period)) {
//					System.out.println("Adding source " + head.source_id + ", rupture " + head.rupture_id + ", rup var " + head.rup_var_id + ", period " + e.period);
					// Initialize PeakAmplitudes class
					ps.setString(1, "" + head.source_id);
					ps.setString(2, "" + head.rupture_id);
					ps.setString(3, "" + head.rup_var_id);
					ps.setString(4, "" + rd100periodValueToIDMap.get(e.period));
					float rd100 = e.rotD100;
					if (convertGtoCM) {
						rd100 = e.rotD100 * G_TO_CM_2;
					}
					ps.setString(5, "" + rd100);
					ps.addBatch();
				}
				if (rd50periodValueToIDMap.containsKey(e.period)) {
					// Initialize PeakAmplitudes class
					ps.setString(1, "" + head.source_id);
					ps.setString(2, "" + head.rupture_id);
					ps.setString(3, "" + head.rup_var_id);
					ps.setString(4, "" + rd50periodValueToIDMap.get(e.period));
					float rd50 = e.rotD50;
					if (convertGtoCM) {
						rd50 = e.rotD50 * G_TO_CM_2;
					}
					ps.setString(5, "" + rd50);
					ps.addBatch();
				}
				counter++;
				if (counter % 1000==0) {
					ps.executeBatch();
				}
			}
			if ((counter%1000)!=0) {
				//do the stragglers
				ps.executeBatch();
			}
			ps.close();
			stmt.close();
		} catch (SQLException sqe) {
			sqe.printStackTrace();
			System.exit(2);
		}
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
				System.out.println(query.getQueryString());
				int typeID = (Integer)(query.list().get(0));
				rd100periodValueToIDMap.put(desiredPeriods.get(i).floatValue(), typeID);
				System.out.println("Adding IM_Type_ID " + typeID + " to list.");
			}
		}
		if (rd50periodValueToIDMap==null) {
			rd50periodValueToIDMap = new HashMap<Float, Integer>();
			for (int i=0; i<desiredPeriods.size(); i++) {
				SQLQuery query = rotdSession.createSQLQuery(rd50Prefix + "IM_Type_Value = " + desiredPeriods.get(i)).addScalar("IM_Type_ID", Hibernate.INTEGER);
				int typeID = (Integer)(query.list().get(0));
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
//					sess.insert(pa);
				} catch (NonUniqueObjectException nuoe) {
					//Occurs if there's a duplicate entry in the PSA file, which can happen on rare occasions.  Because of the Study 15.4 issues, abort if this happens.
					System.err.println("ERROR:  found duplicate entry in file for run_id " + paPK.getRun_ID() + ", source " + paPK.getSource_ID() + " rupture " + paPK.getRupture_ID() + " rup_var " + paPK.getRup_Var_ID() + " IM_Type " + paPK.getIM_Type_ID() + ".  Aborting.");
					System.exit(2);
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
//					sess.insert(pa);
				} catch (NonUniqueObjectException nuoe) {
					//Occurs if there's a duplicate entry in the PSA file, which can happen on rare occasions.  Because of the Study 15.4 issues, abort if this happens.
					System.err.println("ERROR:  found duplicate entry for run_id " + paPK.getRun_ID() + ", source " + paPK.getSource_ID() + " rupture " + paPK.getRupture_ID() + " rup_var " + paPK.getRup_Var_ID() + " IM_Type " + paPK.getIM_Type_ID() + ".  Aborting.");
					System.exit(2);
				}
			}
		}
	}

	private void insertRupture(SARuptureFromRuptureVariationFile saRuptureWithSingleRupVar, Session sess, boolean headers) {
		//open session
//		System.out.println("Determining type IDs.");
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

//		System.out.println("Looping over rupture variations.");
		
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
//			System.out.println("Looping over periods for RV " + currRupVar.variationNumber);
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
					if (psaValue>8400 || psaValue<0.008) {
						//If force insert is on, we don't care
						if (forceInsert) {
							System.out.println("Found psaValue " + psaValue + " for source " + currentSource_ID + " rupture " + currentRupture_ID + " rup var " + currRupVar.variationNumber);
							System.out.println("Force-insert option was selected, inserting anyway.");
						} else {
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
							//We feel OK with values < 0.008 for M6.85 @ 300 km, and M6.05 @ 200 km
							//Based on that, we made a line with some tweaks
							double okDist = 125.0*mag - 576.25;
							if (rupture_dist<okDist || psaValue<0.003) {
								System.err.println("Found value " + psaValue + " for source " + currentSource_ID + ", " + currentRupture_ID + ", " + currRupVar.variationNumber + ", period index " + periodIter + ", period value " + ourPeriods[periodIter]);
								System.err.println("Mag=" + mag + " rupture_dist=" + rupture_dist);
								throw new IllegalArgumentException();
							} else {
								System.out.println("Found value " + psaValue + " for source " + currentSource_ID + ", " + currentRupture_ID + ", " + currRupVar.variationNumber + ", period index " + periodIter + ", period value " + ourPeriods[periodIter]);
								System.out.println("Since mag=" + mag + " and source rupture dist=" + rupture_dist + ", permitting insertion.");
							}
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
						//Occurs if there's a duplicate entry in the PSA file, which can happen on rare occasions.  Because of the issue in Study 15.4, abort.
						System.err.println("ERROR:  found duplicate entry for run_id " + paPK.getRun_ID() + ", source " + paPK.getSource_ID() + " rupture " + paPK.getRupture_ID() + " rup_var " + paPK.getRup_Var_ID() + " IM_Type " + paPK.getIM_Type_ID() + ".  Aborting.");
						System.exit(2);
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
