package test.fitnesse;

import data.DirectionalComponent;
import data.SARuptureFromRuptureVariationFile;
import fit.ColumnFixture;


public class RuptureVariationFixture extends ColumnFixture {
	private SARuptureFromRuptureVariationFile saRupture;
	public DirectionalComponent northComp;
	public DirectionalComponent eastComp;
	
	public DirectionalComponent geomAvgComp() {
		saRupture = new SARuptureFromRuptureVariationFile();
		saRupture.rupVar.northComp = northComp;
		saRupture.rupVar.eastComp = eastComp;
		try {
			saRupture.rupVar.computeGeomAvg();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return saRupture.rupVar.geomAvgComp;
		/*return new DirectionalComponent();*/
	}
}
