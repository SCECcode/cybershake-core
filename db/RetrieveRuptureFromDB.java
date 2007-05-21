import java.sql.ResultSet;
import java.sql.SQLException;



public class RetrieveRuptureFromDB {

	public static void main(String[] args) {
		if (args.length<3) {
			System.out.println("Usage: RetrieveRuptureFromDB ERF_ID Source_ID Rupture_ID");
			System.exit(0);
		}
		String erf = args[0];
		String source = args[1];
		String rupture = args[2];
		
		DBConnect dbc = new DBConnect("intensity.usc.edu","CyberShake");
		ResultSet rs = dbc.selectData("select Prob, Mag, Grid_Spacing, Num_Rows, Num_Cols " + 
				"from Ruptures " +
				"where ERF_ID = " + erf +
				" and Source_ID = " + source +
				" and Rupture_ID = " + rupture);
		RuptureFileData rfd;
		try {
			rs.first();
			while (!rs.isAfterLast()) {
				rfd = new RuptureFileData(rs.getDouble("Prob"), rs.getDouble("Mag"), rs.getDouble("Grid_Spacing"), rs.getInt("Num_Rows"), rs.getInt("Num_Cols"));
				
				rs.next();
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}
}

