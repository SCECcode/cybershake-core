package mapping;

public class PeakAmplitudes {
	private PeakAmplitudesPK paPK;
	private double IM_Value;
	private String Units;
	
	public double getIM_Value() {
		return IM_Value;
	}
	public void setIM_Value(double value) {
		IM_Value = value;
	}
	public String getUnits() {
		return Units;
	}
	public void setUnits(String units) {
		Units = units;
	}
	public PeakAmplitudesPK getPaPK() {
		return paPK;
	}
	public void setPaPK(PeakAmplitudesPK paPK) {
		this.paPK = paPK;
	}
	
	public boolean equals(Object obj) {
		if(this == obj)
			return true;
		if((obj == null) || (obj.getClass() != this.getClass()))
			return false;
		// object must be PeakAmplitudesPK at this point
		PeakAmplitudes pa = (PeakAmplitudes)obj;
		//return (paPK.equals(pa.getPaPK())  && IM_Value == pa.getIM_Value() && Units == pa.getUnits());
		return (paPK.equals(pa.getPaPK()));
	}
	
	public int hashCode() {
		return 0;
	}
	
	public String toString() {
		return paPK.toString() + ", IM_Value: " + IM_Value + ", Units: " + Units;
	}
	
}
