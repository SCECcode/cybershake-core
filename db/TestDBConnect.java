import java.sql.ResultSet;
import java.sql.SQLException;



public class TestDBConnect {

	public static void main(String[] args) {
		DBConnect dbc = new DBConnect("surface.usc.edu","CyberShake");
		ResultSet rs = dbc.selectData("SHOW TABLES");
		try {
			System.out.println(rs.getMetaData().getColumnCount());
			rs.first();
			System.out.println(rs.getRow());
			while (!rs.isAfterLast()) {
				System.out.println(rs.getString(1));
				rs.next();
			}
            rs.close();
            rs = dbc.selectData("select count(*) from CyberShake_Sites");
            rs.getInt("count(*)");
		} catch (SQLException e) {
			e.printStackTrace();
		}
	} 
}
