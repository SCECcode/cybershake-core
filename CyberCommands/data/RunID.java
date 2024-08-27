package data;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

public class RunID {
	private int runID;
	private int erfID;
	private int siteID = -1;
	private String siteName;
	private int sgtVarID;
	private int ruptVarScenID;
	private double cutoffDist;
	private Connection connection = null;
	
	private String HOSTNAME = "focal.usc.edu";
	private boolean isSQLite = false;
	private final int PORT = 3306;
	private final String DB_NAME = "CyberShake";
	private final String USER = "cybershk";
	private final String PASS = "phy$ic@1St3ady";
	
	public RunID(int runID, String server) {
		if (server.equals("moment")) {
			HOSTNAME = "moment.usc.edu";
		} else if (server.equals("focal")) {
			HOSTNAME = "focal.usc.edu";
		} else if (server.equals("moment_carc")) {
			HOSTNAME = "moment";
		} else if (server.indexOf("sqlite:")==0) {
			HOSTNAME = server;
			isSQLite = true;
		}
		this.runID = runID;
		createConnection();
		populateRunIDInfo();
		populateSiteInfo();
		closeConnection();
	}
	
	public RunID(int runID) {
		this.runID = runID;
		createConnection();
		populateRunIDInfo();
		populateSiteInfo();
		closeConnection();
	}

	private void closeConnection() {
		if (connection!=null) {
			try {
				connection.close();
			} catch (SQLException e) {
				e.printStackTrace();
				System.exit(-3);
			}
		}
	}

	private void createConnection() {
		if (connection==null) {
			String drivers = "com.mysql.cj.jdbc.Driver";
			//Have to add serverTimezone now; assume 
			String url = "jdbc:mysql://"+HOSTNAME+":"+PORT+"/"+DB_NAME + "?serverTimezone=America/Los_Angeles";
			//System.out.println("Connecing with url " + url);
			if (isSQLite) {
				//Use SQLite drivers
				url = "jdbc:" + HOSTNAME;
			}
	    
			try {
				Class.forName(drivers).newInstance();
			} catch (Exception ex) {
				ex.printStackTrace();
				System.exit(-1);
			}
			try {
				if (isSQLite) {
					connection = DriverManager.getConnection(url);
				} else {
					connection = DriverManager.getConnection(url,USER,PASS);
				}
			} catch (SQLException ex) {
				ex.printStackTrace();
				System.exit(-2);
			}
		}
	}

	private void populateSiteInfo() {
		try {
			Statement stat = connection.createStatement();
			String query;
			if (siteID==-1) {
				populateRunIDInfo();
			}
			query = "SELECT S.CS_Short_Name, R.Cutoff_Dist FROM CyberShake_Sites S, CyberShake_Site_Regions R WHERE R.ERF_ID=" + erfID + " and R.CS_Site_ID=S.CS_Site_ID and S.CS_Site_ID=" + siteID;
			System.out.println(query);
			ResultSet res = stat.executeQuery(query);
			if (!isSQLite) {
				res.first();
			} else {
				res.next();
			}
			if (res.getRow()==0 || res.isClosed()) {
				System.err.println("Couldn't find a site name for site ID " + siteID);
				System.exit(3);
			}
			String shortNameField = "S.CS_Short_Name";
			String cutoffDistField = "R.Cutoff_Dist";
			if (isSQLite) {
				shortNameField = "CS_Short_Name";
				cutoffDistField = "Cutoff_Dist";
			}
			siteName = res.getString(shortNameField);
			cutoffDist = res.getDouble(cutoffDistField);
    		boolean nextRow = res.next();
    		if ((isSQLite && nextRow) || (!isSQLite && !res.isAfterLast())) {
    			System.err.println("More than one Run_ID matched Run_ID "  + runID + ", aborting.");
    			System.exit(3);
    		}
    		res.close();
    		stat.close();
		} catch (Exception ex) {
			ex.printStackTrace();
			System.exit(2);
		}
	}

	private void populateRunIDInfo() {
		try {
			Statement stat = connection.createStatement();
			String query = "SELECT Site_ID, ERF_ID, SGT_Variation_ID, Rup_Var_Scenario_ID FROM CyberShake_Runs WHERE Run_ID=" + runID;
			ResultSet res = stat.executeQuery(query);
			if (!isSQLite) {
				res.first();
			} else {
				res.next();
			}
			
    		if (res.getRow()==0 || res.isClosed()) {
    			System.err.println("Couldn't find run " + runID + ".");
    			System.exit(1);
    		}
    		siteID = res.getInt("Site_ID");
    		erfID = res.getInt("ERF_ID");
    		sgtVarID = res.getInt("SGT_Variation_ID");
    		ruptVarScenID = res.getInt("Rup_Var_Scenario_ID");
    		boolean nextRow = res.next();
    		if ((isSQLite && nextRow) || (!isSQLite && !res.isAfterLast())) {
    			System.err.println("More than one Run_ID matched Run_ID "  + runID + ", aborting.");
    			System.exit(3);
    		}
    		res.close();
    		stat.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	public int getRunID() {
		return runID;
	}

	public int getErfID() {
		return erfID;
	}

	public int getSiteID() {
		return siteID;
	}

	public String getSiteName() {
		return siteName;
	}

	public int getSgtVarID() {
		return sgtVarID;
	}

	public int getRuptVarScenID() {
		return ruptVarScenID;
	}

	public double getCutoffDist() {
		return cutoffDist;
	}
	
}
