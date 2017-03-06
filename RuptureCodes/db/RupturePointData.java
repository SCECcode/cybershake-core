import java.io.PrintStream;



public class RupturePointData {
	private double lat;
	private double lon;
	private double depth;
	private double rake;
	private double dip;
	private double strike;
	
	public RupturePointData(double la, double lo, double de, double r, double di, double s) {
		lat = la;
		lon = lo;
		depth = de;
		rake = r;
		dip = di;
		strike = s;
	}

	public void print(PrintStream ps) {
		ps.println(lat + "    " + lon + "    " + depth + "   " + rake + "    " + dip + "   " + strike);
	}
	
	
}
