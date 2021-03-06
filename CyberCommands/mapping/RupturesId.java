package mapping;

// Generated Aug 9, 2007 9:53:47 PM by Hibernate Tools 3.2.0.b9

/**
 * RupturesId generated by hbm2java
 */
public class RupturesId implements java.io.Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -380203655321709518L;

	private int erfId;

	private int sourceId;

	private int ruptureId;

	public RupturesId() {
	}

	public RupturesId(int erfId, int sourceId, int ruptureId) {
		this.erfId = erfId;
		this.sourceId = sourceId;
		this.ruptureId = ruptureId;
	}

	public int getErfId() {
		return this.erfId;
	}

	public void setErfId(int erfId) {
		this.erfId = erfId;
	}

	public int getSourceId() {
		return this.sourceId;
	}

	public void setSourceId(int sourceId) {
		this.sourceId = sourceId;
	}

	public int getRuptureId() {
		return this.ruptureId;
	}

	public void setRuptureId(int ruptureId) {
		this.ruptureId = ruptureId;
	}

	public boolean equals(Object other) {
		if ((this == other))
			return true;
		if ((other == null))
			return false;
		if (!(other instanceof RupturesId))
			return false;
		RupturesId castOther = (RupturesId) other;

		return (this.getErfId() == castOther.getErfId())
				&& (this.getSourceId() == castOther.getSourceId())
				&& (this.getRuptureId() == castOther.getRuptureId());
	}

	public int hashCode() {
		int result = 17;

		result = 37 * result + this.getErfId();
		result = 37 * result + this.getSourceId();
		result = 37 * result + this.getRuptureId();
		return result;
	}

}
