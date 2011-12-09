import java.io.File;
import java.io.FileNotFoundException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Scanner;
import java.io.PrintStream;



public class RetrieveRuptureFromDB {
	private static DBConnect dbc;
	private static final String DB_SERVER = "focal.usc.edu";
	private static final String DB = "CyberShake";
	private static String ERF_ID;
        private static final String RuptureFile = "ruptures.list";

	public static void main(String[] args) {
		dbc = new DBConnect(DB_SERVER,DB);
		if (args.length==2) {
			if (args[1].equals("all")) {
				ERF_ID = args[0];
				getAllRuptures("");
				dbc.closeConnection();
				return;
			}
		} else if (args.length==3) {
			if (args[1].equals("all")) {
				ERF_ID = args[0];
				getAllRuptures(args[2]);
				dbc.closeConnection();
				return;
			} else if (args[1].equals("-f")) {
				ERF_ID = args[0];
				getSomeRuptures(args[2]);
				dbc.closeConnection();
				return;
			} else {
				ERF_ID = args[0];
				String source = args[1];
				String rupture = args[2];
				getOneRupture(source, rupture, null, null, null);
				dbc.closeConnection();
				return;
			}
		} else if (args.length==4) {
			if (args[1].equals("-f")) {
				ERF_ID = args[0];
				String path = args[3];
				getSomeRuptures(args[2], path);
				dbc.closeConnection();
				return;
			} else {
				ERF_ID = args[0];
				String source = args[1];
				String rupture = args[2];
				String filename = args[3];
				getOneRupture(source, rupture, filename, null, null);
				dbc.closeConnection();
				return;
			}
		} else {
			dbc.closeConnection();
			System.out.println("Usage: RetrieveRuptureFromDB <ERF_ID> (<all> | (<Source_ID> <Rupture_ID>) | -f <filename>) [path])");
			System.out.println("File format:\n\tERF_ID\n\tsource1 rupture1\n\tsource2 rupture2\n\t...");
			System.exit(0);
		}
	}
	
	private static void getSomeRuptures(String inputFile, String path) {
		/* file format:
		 * 
		 * ERF_ID
		 * source1 rupture1
		 * source2 rupture2
		 * ...
		 */
		try {
			Scanner input = new Scanner(new File(inputFile));
			ERF_ID = input.next();
			String source, rupture, filename;
			while (input.hasNext()) {
				source = input.next();
				rupture = input.next();
				if (path=="") {
					filename = source + "_" + rupture + ".txt";
				} else {
					filename = path + File.separatorChar + source + "_" + rupture + ".txt";
				}
				getOneRupture(source, rupture, filename, null, null);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
	}

	private static void getSomeRuptures(String inputFile) {
		getSomeRuptures(inputFile, ".");
	}

	private static void getOneRupture(String source, String rupture, String filename, PrintStream index, 
					  String indexentry) {
		RuptureFileData rfd = getRuptureFileData(source, rupture);
		if (rfd==null) {
			System.exit(0);
		}
		getPointsData(rfd, source, rupture);
		if (filename!=null) {
			rfd.printToFile(filename);
			if ((index != null) && (indexentry != null)) {
			    index.println(source + " " + rupture + " " + rfd.getMagnitude() + " " + indexentry);
			}
		} else {
			rfd.print();
		}
	}

	private static void getAllRuptures(String path) {
//		dbc = new DBConnect(DB_SERVER,DB);
		String select = "select Source_ID, Rupture_ID from Ruptures where ERF_ID=" + ERF_ID;
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
			        // Create index file
			        try {
			        PrintStream index = new PrintStream(new File(path + File.separatorChar + RuptureFile));
				while (!rs.isAfterLast()) {
					if (rs.getRow()%1000==0) {
						System.gc();
					}
					String source_id = rs.getString("Source_ID");
					String rupture_id = rs.getString("Rupture_ID");
					String filename;
					String indexentry = null;
					if (path=="") {
						filename = source_id + "_" + rupture_id + ".txt";
					} else {
						File file = new File(path + File.separatorChar + source_id + 
								     File.separatorChar + rupture_id);
						file.mkdirs();
						indexentry = File.separatorChar + source_id + File.separatorChar + rupture_id + 
						    File.separatorChar + source_id + "_" + rupture_id + ".txt";
						filename = path + indexentry;
					}
					System.out.println("Creating file " + filename + ", " + rs.getRow() + " of " + totalRows);
					if (path=="") {
					    getOneRupture(source_id, rupture_id, filename, null, null);
					} else {
					    getOneRupture(source_id, rupture_id, filename, index, indexentry);
					}
					rs.next();
				}
				index.close();
		                } catch (FileNotFoundException e) {
				        e.printStackTrace();
					System.out.println(e.getMessage());
		                }

			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	private static RuptureFileData getRuptureFileData(String source, String rupture) {
		String select = "select Prob, Mag, Grid_Spacing, Num_Rows, Num_Columns, Num_Points " + 
		"from Ruptures " +
		"where ERF_ID = " + ERF_ID +
		" and Source_ID = " + source +
		" and Rupture_ID = " + rupture;
		ResultSet rs = dbc.selectData(select);
		RuptureFileData rfd = null;
		try {
			rs.first();
			if (rs.getRow()==0) {
				System.out.println("No hits in the db for ERF_ID " + ERF_ID + ", Source_ID " + source + ", Rupture_ID " + rupture + ", exiting.");
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

	private static void getPointsData(RuptureFileData rfd, String source, String rupture) {
		ResultSet rs = dbc.selectData("select Lat, Lon, Depth, Rake, Dip, Strike " +
				"from Points " +
				"where ERF_ID = " + ERF_ID +
				" and Source_ID = " + source +
				" and Rupture_ID = " + rupture);
		RupturePointData rpd;
//		rfd.clearPoints();
		try {
			rs.first();
			if (rs.getRow()==0) {
				System.err.println("No points in the db for ERF_ID " + ERF_ID + ", Source_ID " + source + ", Rupture_ID " + rupture + ", exiting.");
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

