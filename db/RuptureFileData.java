import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;

public class RuptureFileData {
	private double probability;
	private double magnitude;
	private double gridSpacing;
	private int numRows;
	private int numCols;
	private ArrayList<RupturePointData> points = new ArrayList<RupturePointData>();
	
	public RuptureFileData(double pr, double ma, double gr, int nr, int nc) {
		probability = pr;
		magnitude = ma;
		gridSpacing = gr;
		numRows = nr;
		numCols = nc;
	}
	
	public RuptureFileData(double pr, double ma, double gr, int nr, int nc, ArrayList<RupturePointData> pts) {
		this(pr, ma, gr, nr, nc);
		points = pts;
	}
	
	public void addPoint(RupturePointData rpd) {
		points.add(rpd);
	}
	
	public void print() {
		print(System.out);
	}
	
	public void printToFile(String filename) {
		try {
			PrintStream ps = new PrintStream(new File(filename));
			print(ps);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	private void print(PrintStream ps) {
		ps.println("Probability = " + probability);
		ps.println("Magnitude = " + magnitude);
		ps.println("GridSpacing = " + gridSpacing);
		ps.println("NumRows = " + numRows);
		ps.println("NumCols = " + numCols);
		ps.println("#   Lat         Lon         Depth      Rake    Dip     Strike");
		for (RupturePointData data: points) {
			data.print(ps);
		}
	}
}
