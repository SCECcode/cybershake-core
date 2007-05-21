import java.sql.ResultSet;
import java.sql.SQLException;



public class TestDBConnect {

	public static void main(String[] args) {
		DBConnect dbc = new DBConnect("intensity.usc.edu","CyberShake");
		ResultSet rs = dbc.selectData("SHOW TABLES");
		try {
			System.out.println(rs.getMetaData().getColumnCount());
			rs.first();
			while (!rs.isAfterLast()) {
				System.out.println(rs.getString(1));
				rs.next();
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
	} 
}
