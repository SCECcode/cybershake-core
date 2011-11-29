import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;



public class AddPointsToDB {
	private static String erf;
	private static DBConnect dbc;
	private static final String DB_SERVER = "intensity.usc.edu";
	private static final String DB = "CyberShake";
	
	public static void main(String[] args) {
		if (args.length==2) {
			erf = args[0];
			String directory = args[1];
			populate(directory);
			return;
		} else if (args.length<4) {
			System.out.println("Usage: AddPointsToDB <ERF_ID> <Source_ID> <Rupture_ID> <filename>");
			System.exit(0);
		}
		
		erf = args[0];
		String source = args[1];
		String rupture = args[2];
		String filename = args[3];
		
		dbc = new DBConnect(DB_SERVER,DB);
		
		populate(source, rupture, filename);
		dbc.closeConnection();
	}
	
	private static void populate(String directory) {
		File dir = new File(directory);
		if (!dir.isDirectory()) {
			System.err.println(directory + " is not a directory.");
			System.exit(1);
		}
		
		File[] files = dir.listFiles();
		for (int i=0; i<files.length; i++) {
			String name = files[i].getName();
			String[] pieces = name.split("\\.");
			String prefix = pieces[0];
			pieces = prefix.split("_");
			System.out.println("Inserting source " + pieces[0] + ", rupture " + pieces[1] + " filename " + name);
			populate(pieces[0], pieces[1], files[i].getAbsolutePath());
		}
	}
	
	private static void populate(String source, String rupture, String filename) {
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
//				System.out.println("inserting point " + i);
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
