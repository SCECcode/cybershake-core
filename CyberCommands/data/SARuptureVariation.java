package data;

public class SARuptureVariation {
	public DirectionalComponent northComp = new DirectionalComponent();
	public DirectionalComponent eastComp = new DirectionalComponent();
	public DirectionalComponent geomAvgComp = new DirectionalComponent();
	public int variationNumber;
	
	public void computeGeomAvg() throws Exception {
		if (northComp.periods.length != eastComp.periods.length) {
			throw new Exception("Directional Component has north and east components of different lengths");
		}
		for (int i=0; i<northComp.periods.length; i++) {
			geomAvgComp.periods[i] = (float) Math.sqrt(eastComp.periods[i]*northComp.periods[i]);
			//System.out.println("Rupture Variation " + variationNumber + ", Period " + (i + 1) + " : Square root of (" + eastComp.periods[i] + "*" + eastComp.periods[i] + "+" + northComp.periods[i] + "*" + northComp.periods[i] + ") = " + geomAvgComp.periods[i]);
		}
	}
	
	
}
