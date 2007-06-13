import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;



public class AddPointsToDB {
	private static String erf;
	private static String source;
	private static String rupture;
	private static String filename;
	private static DBConnect dbc;
	private static final String DB_SERVER = "intensity.usc.edu";
	private static final String DB = "CyberShake";
	
	public static void main(String[] args) {
		if (args.length<4) {
			System.out.println("Usage: AddPointsToDB <ERF_ID> <Source_ID> <Rupture_ID> <filename>");
			System.exit(0);
		}
		
		erf = args[0];
		source = args[1];
		rupture = args[2];
		filename = args[3];
		
		dbc = new DBConnect(DB_SERVER,DB);
		
		populate();
	}
	
	private static void populate() {
		BufferedReader br;
		String insertString = "insert into Points (ERF_ID, Source_ID, Rupture_ID, Lat, Lon, Depth, Rake, Dip, Strike) ";
		String valueString;
		try {
			br = new BufferedReader(new FileReader(filename));
			br.readLine(); //prob
			br.readLine(); //mag
			br.readLine(); //spacing
			br.readLine(); //rows
			br.readLine(); //cols
			br.readLine(); //comment line
			String line = br.readLine();
			int i=1;
			while (line!=null) {
				String[] pieces = line.split("\\s+");
				//Point_ID, ERF_ID, Source_ID, Rupture_ID, Lat, Lon, Depth, Rake, Dip, Strike
				valueString = "values (" + erf + ", " + source + ", " + rupture + ", " +
				pieces[0] + ", " + pieces[1] + ", " + pieces[2] + ", " + pieces[3] + ", " + pieces[4] + ", " + pieces[5] + ")";
				if (!dbc.insertData(insertString + valueString)) {
					System.err.println("Error inserting data with query " + insertString + valueString);
				}
				System.out.println("inserting point " + i);
				line = br.readLine();
				i++;
			}
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
