package faultlist;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Calendar;


public class CreateFaultList {
	private static final String DB_SERVER = "intensity.usc.edu";
	private static final String DB = "CyberShake";

	public static void main(String[] args) {
		if (args.length<4) {
			System.out.println("Usage:  CreateFaultList <site> <erf_id> <path to rup vars> <outputFile>");
			System.exit(1);
		}
		createList(args[0], Integer.parseInt(args[1]), args[2], args[3]);
		System.exit(0);
	}

	private static void createList(String site, int erf_id, String pathToVars, String output) {
		DBConnect dbc = new DBConnect(DB_SERVER, DB);
		String query = "select R.Source_ID, R.Rupture_ID " +
		"from CyberShake_Site_Ruptures R, CyberShake_Sites S " +
		"where S.CS_Short_Name=\"" + site + "\" " +
		"and R.CS_Site_ID=S.CS_Site_ID " +
		"and R.ERF_ID=" + erf_id + " " + 
		"order by R.Source_ID, R.Rupture_ID";
		System.out.println(query);
		ResultSet ruptures = dbc.selectData(query);
		try {
			ruptures.first();
			if (ruptures.getRow()==0) {
				System.err.println("No ruptures found for site " + site + ", exiting.");
				System.exit(2);
			}
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(output));
			int count = 0;
			while (!ruptures.isAfterLast()) {
				//write path/sourceID/ruptureID/sourceID_ruptureID.txt
				int sourceID = ruptures.getInt("Source_ID");
				int ruptureID = ruptures.getInt("Rupture_ID");
				bw.write(pathToVars + File.separator + sourceID + File.separator + ruptureID + File.separator + sourceID + "_" + ruptureID + ".txt nheader=6 latfirst=1\n");
				ruptures.next();
				count++;
				if (count%100==0) System.out.println("Processed " + count + " ruptures.");
			}
			bw.flush();
			bw.close();
			ruptures.close();
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
