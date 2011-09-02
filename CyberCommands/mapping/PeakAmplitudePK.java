package mapping;

import java.io.Serializable;

public class PeakAmplitudePK implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int Source_ID;
	private int Rupture_ID;
	private long Rup_Var_ID;
	private int Run_ID;
	private int IM_Type_ID;

	public int getIM_Type_ID() {
		return IM_Type_ID;
	}
	public void setIM_Type_ID(int type_ID) {
		IM_Type_ID = type_ID;
	}
	public long getRup_Var_ID() {
		return Rup_Var_ID;
	}
	public void setRup_Var_ID(long rup_Var_ID) {
		Rup_Var_ID = rup_Var_ID;
	}
	
	public int getRupture_ID() {
		return Rupture_ID;
	}
	public void setRupture_ID(int rupture_ID) {
		Rupture_ID = rupture_ID;
	}
	public int getSource_ID() {
		return Source_ID;
	}
	public void setSource_ID(int source_ID) {
		Source_ID = source_ID;
	}
	
	public int getRun_ID() {
		return Run_ID;
	}
	public void setRun_ID(int run_ID) {
		Run_ID = run_ID;
	}
	
	public boolean equals(Object obj)
	{
		if(this == obj)
			return true;
		if((obj == null) || (obj.getClass() != this.getClass()))
			return false;
		// object must be PeakAmplitudesPK at this point
		PeakAmplitudePK papk = (PeakAmplitudePK)obj;
		return (Source_ID == papk.getSource_ID() &&
				Rupture_ID == papk.getRupture_ID() && 
				Rup_Var_ID == papk.getRup_Var_ID() &&
				Run_ID == papk.getRun_ID() &&
				IM_Type_ID == papk.getIM_Type_ID() );
	}
	
	public int hashCode() {
		return 0;
	}
	
	public String toString() {
		return "Run_ID: " + Run_ID + ", Source_ID: " + Source_ID + ", Rupture_ID: " + Rupture_ID
		+ ", Rup_Var_ID: " + Rup_Var_ID + ", IM_Type_ID: " + IM_Type_ID;
	}

}
