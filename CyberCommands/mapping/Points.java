package mapping;

// Generated Aug 9, 2007 9:53:47 PM by Hibernate Tools 3.2.0.b9

/**
 * Points generated by hbm2java
 */
public class Points implements java.io.Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -8560105081734393786L;

	private long pointId;

	private ErfIds erfIds;

	private Ruptures ruptures;

	private double lat;

	private double lon;

	private double depth;

	private double rake;

	private double dip;

	private double strike;

	public Points() {
	}

	public Points(long pointId, ErfIds erfIds, Ruptures ruptures, double lat,
			double lon, double depth, double rake, double dip, double strike) {
		this.pointId = pointId;
		this.erfIds = erfIds;
		this.ruptures = ruptures;
		this.lat = lat;
		this.lon = lon;
		this.depth = depth;
		this.rake = rake;
		this.dip = dip;
		this.strike = strike;
	}

	public long getPointId() {
		return this.pointId;
	}

	public void setPointId(long pointId) {
		this.pointId = pointId;
	}

	public ErfIds getErfIds() {
		return this.erfIds;
	}

	public void setErfIds(ErfIds erfIds) {
		this.erfIds = erfIds;
	}

	public Ruptures getRuptures() {
		return this.ruptures;
	}

	public void setRuptures(Ruptures ruptures) {
		this.ruptures = ruptures;
	}

	public double getLat() {
		return this.lat;
	}

	public void setLat(double lat) {
		this.lat = lat;
	}

	public double getLon() {
		return this.lon;
	}

	public void setLon(double lon) {
		this.lon = lon;
	}

	public double getDepth() {
		return this.depth;
	}

	public void setDepth(double depth) {
		this.depth = depth;
	}

	public double getRake() {
		return this.rake;
	}

	public void setRake(double rake) {
		this.rake = rake;
	}

	public double getDip() {
		return this.dip;
	}

	public void setDip(double dip) {
		this.dip = dip;
	}

	public double getStrike() {
		return this.strike;
	}

	public void setStrike(double strike) {
		this.strike = strike;
	}

}