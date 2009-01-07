package mapping;

public class PeakAmplitude {
	private PeakAmplitudePK paPK;
	private double IM_Value;
	
	public double getIM_Value() {
		return IM_Value;
	}
	public void setIM_Value(double value) {
		IM_Value = value;
	}
	public PeakAmplitudePK getPaPK() {
		return paPK;
	}
	public void setPaPK(PeakAmplitudePK paPK) {
		this.paPK = paPK;
	}
	
	public boolean equals(Object obj) {
		if(this == obj)
			return true;
		if((obj == null) || (obj.getClass() != this.getClass()))
			return false;
		// object must be PeakAmplitudesPK at this point
		PeakAmplitude pa = (PeakAmplitude)obj;
		//return (paPK.equals(pa.getPaPK())  && IM_Value == pa.getIM_Value() && Units == pa.getUnits());
		return (paPK.equals(pa.getPaPK()));
	}
	
	public int hashCode() {
		return 0;
	}
	
	public String toString() {
		return paPK.toString() + ", IM_Value: " + IM_Value;
	}
	
}
