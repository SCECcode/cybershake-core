package mapping;

import java.io.Serializable;

public class PeakAmplitudePK implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int ERF_ID;
	private int Site_ID; 
	private int Source_ID;
	private int Rupture_ID;
	private long Rup_Var_ID;
	private int Rup_Var_Scenario_ID;
	private String IM_Type;
	
	public int getERF_ID() {
		return ERF_ID;
	}
	public void setERF_ID(int erf_id) {
		ERF_ID = erf_id;
	}
	public String getIM_Type() {
		return IM_Type;
	}
	public void setIM_Type(String type) {
		IM_Type = type;
	}
	public long getRup_Var_ID() {
		return Rup_Var_ID;
	}
	public void setRup_Var_ID(long rup_Var_ID) {
		Rup_Var_ID = rup_Var_ID;
	}
	public int getRup_Var_Scenario_ID() {
		return Rup_Var_Scenario_ID;
	}
	public void setRup_Var_Scenario_ID(int rup_Var_Scenario_ID) {
		Rup_Var_Scenario_ID = rup_Var_Scenario_ID;
	}
	public int getRupture_ID() {
		return Rupture_ID;
	}
	public void setRupture_ID(int rupture_ID) {
		Rupture_ID = rupture_ID;
	}
	public int getSite_ID() {
		return Site_ID;
	}
	public void setSite_ID(int site_ID) {
		Site_ID = site_ID;
	}
	public int getSource_ID() {
		return Source_ID;
	}
	public void setSource_ID(int source_ID) {
		Source_ID = source_ID;
	}
	
	public boolean equals(Object obj)
	{
		if(this == obj)
			return true;
		if((obj == null) || (obj.getClass() != this.getClass()))
			return false;
		// object must be PeakAmplitudesPK at this point
		PeakAmplitudePK papk = (PeakAmplitudePK)obj;
		return (ERF_ID == papk.getERF_ID() &&
				Site_ID == papk.getSite_ID() &&
				Source_ID == papk.getSource_ID() &&
				Rupture_ID == papk.getRupture_ID() && 
				Rup_Var_ID == papk.getRup_Var_ID() &&
				Rup_Var_Scenario_ID == papk.getRup_Var_Scenario_ID() &&
				IM_Type.equals(papk.getIM_Type()) );
	}
	
	public int hashCode() {
		return 0;
	}
	
	public String toString() {
		return "ERF_ID: " + ERF_ID + ", Site_ID: " + Site_ID + ", Source_ID: " + Source_ID + ", Rupture_ID: " + Rupture_ID
		+ ", Rup_Var_ID: " + Rup_Var_ID + ", Rup_Var_Scenario_ID: " + Rup_Var_Scenario_ID + ", IM_Type: " + IM_Type;
	}
}
