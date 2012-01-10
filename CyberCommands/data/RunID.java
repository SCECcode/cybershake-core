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
	
	private final String HOSTNAME = "focal.usc.edu";
	private final int PORT = 3306;
	private final String DB_NAME = "CyberShake";
	private final String USER = "cybershk_ro";
	private final String PASS = "CyberShake2007";
	
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
			String drivers = "com.mysql.jdbc.Driver";
			String url = "jdbc:mysql://"+HOSTNAME+":"+PORT+"/"+DB_NAME;
	    
			try {
				Class.forName(drivers).newInstance();
			} catch (Exception ex) {
				ex.printStackTrace();
				System.exit(-1);
			}
			try {
				connection = DriverManager.getConnection(url,USER,PASS);
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
			query = "SELECT S.CS_Short_Name, R.Cutoff_Dist FROM CyberShake_Sites S, CyberShake_Site_Regions R WHERE R.CS_Site_ID=S.CS_Site_ID and S.CS_Site_ID=" + siteID;
			ResultSet res = stat.executeQuery(query);
			res.first();
			if (res.getRow()==0) {
				System.err.println("Couldn't find a site name for site ID " + siteID);
				System.exit(3);
			}
			siteName = res.getString("S.CS_Short_Name");
			cutoffDist = res.getDouble("R.Cutoff_Dist");
    		res.next();
    		if (!res.isAfterLast()) {
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
			res.first();
    		if (res.getRow()==0) {
    			System.err.println("Couldn't find run " + runID + ".");
    			System.exit(1);
    		}
    		siteID = res.getInt("Site_ID");
    		erfID = res.getInt("ERF_ID");
    		sgtVarID = res.getInt("SGT_Variation_ID");
    		ruptVarScenID = res.getInt("Rup_Var_Scenario_ID");
    		res.next();
    		if (!res.isAfterLast()) {
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
