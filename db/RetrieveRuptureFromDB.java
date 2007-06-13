import java.io.File;
import java.sql.ResultSet;
import java.sql.SQLException;



public class RetrieveRuptureFromDB {
	private static DBConnect dbc;
	private static final String DB_SERVER = "intensity.usc.edu";
	private static final String DB = "CyberShake";
	
	public static void main(String[] args) {
		dbc = new DBConnect(DB_SERVER,DB);
		if (args.length==1) {
			if (args[0].equals("all")) {
				getAllRuptures("");
				dbc.closeConnection();
				return;
			}
		} else if (args.length==2) {
			if (args[0].equals("all")) {
				getAllRuptures(args[1]);
				dbc.closeConnection();
				return;
			} else {
				//TODO:  read a file with erf, source, rupture, rupture file entries
			}
		} else if (args.length<3) {
			dbc.closeConnection();
			System.out.println("Usage: RetrieveRuptureFromDB <ERF_ID> <Source_ID> <Rupture_ID> [filename]");
			System.exit(0);
		} else {
			String erf = args[0];
			String source = args[1];
			String rupture = args[2];
			String filename = null;
			if (args.length==4) {
				filename = args[3];
			}
		
			getOneRupture(erf, source, rupture, filename);
			dbc.closeConnection();

		}
	}
	
	private static void getOneRupture(String erf, String source, String rupture, String filename) {
		RuptureFileData rfd = getRuptureFileData(erf, source, rupture);
		if (rfd==null) {
			System.exit(0);
		}
		getPointsData(rfd, erf, source, rupture);
		if (filename!=null) {
			rfd.printToFile(filename);
		} else {
			rfd.print();
		}
	}

	private static void getAllRuptures(String path) {
//		dbc = new DBConnect(DB_SERVER,DB);
		String select = "select ERF_ID, Source_ID, Rupture_ID from Ruptures";
		ResultSet rs = dbc.selectData(select);
		int totalRows = 0;
		
		try {
			rs.last();
			totalRows = rs.getRow();
			rs.first();
			if (rs.getRow()==0) {
				System.err.println("No ruptures in table.");
				System.exit(0);
			} else {
				while (!rs.isAfterLast()) {
					if (rs.getRow()%1000==0) {
						System.gc();
					}
					String erf_id = rs.getString("ERF_ID");
					String source_id = rs.getString("Source_ID");
					String rupture_id = rs.getString("Rupture_ID");
					String filename;
					if (path=="") {
						filename = source_id + "_" + rupture_id + ".txt";
					} else {
						filename = path + File.separatorChar + source_id + "_" + rupture_id + ".txt";
					}
					System.out.println("Creating file " + filename + ", " + rs.getRow() + " of " + totalRows);
					getOneRupture(erf_id, source_id, rupture_id, filename);
					rs.next();
				}
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	private static RuptureFileData getRuptureFileData(String erf, String source, String rupture) {
		String select = "select Prob, Mag, Grid_Spacing, Num_Rows, Num_Columns " + 
		"from Ruptures " +
		"where ERF_ID = " + erf +
		" and Source_ID = " + source +
		" and Rupture_ID = " + rupture;
		ResultSet rs = dbc.selectData(select);
		RuptureFileData rfd = null;
		try {
			rs.first();
			if (rs.getRow()==0) {
				System.out.println("No hits in the db for ERF_ID " + erf + ", Source_ID " + source + ", Rupture_ID " + rupture + ", exiting.");
				rs.close();
			} else {
//				rfd.setProbability(rs.getDouble("Prob"));
//				rfd.setMagnitude(rs.getDouble("Mag"));
//				rfd.setGridSpacing(rs.getDouble("Grid_Spacing"));
//				rfd.setNumRows(rs.getInt("Num_Rows"));
//				rfd.setNumCols(rs.getInt("Num_Columns"));				
				rfd = new RuptureFileData(rs.getDouble("Prob"), rs.getDouble("Mag"), rs.getDouble("Grid_Spacing"), rs.getInt("Num_Rows"), rs.getInt("Num_Columns"));
				rs.close();
				return rfd;
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return null;
	}

	private static void getPointsData(RuptureFileData rfd, String erf, String source, String rupture) {
		ResultSet rs = dbc.selectData("select Lat, Lon, Depth, Rake, Dip, Strike " +
				"from Points " +
				"where ERF_ID = " + erf +
				" and Source_ID = " + source +
				" and Rupture_ID = " + rupture);
		RupturePointData rpd;
//		rfd.clearPoints();
		try {
			rs.first();
			if (rs.getRow()==0) {
				System.err.println("No points in the db for ERF_ID " + erf + ", Source_ID " + source + ", Rupture_ID " + rupture + ", exiting.");
				System.exit(1);
			} else {
				while (!rs.isAfterLast()) {
					rpd = new RupturePointData(rs.getDouble("Lat"), rs.getDouble("Lon"), rs.getDouble("Depth"), rs.getDouble("Rake"), rs.getDouble("Dip"), rs.getDouble("Strike"));
					rfd.addPoint(rpd);
					rs.next();
				}
			}
			rs.close();
		} catch (SQLException e) {
			e.printStackTrace();
			System.exit(0);
		}
	}
}

