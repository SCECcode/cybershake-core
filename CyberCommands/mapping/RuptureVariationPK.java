package mapping;

import java.io.Serializable;

public class RuptureVariationPK implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private long Rup_Var_ID;
	private int Rup_Var_Scenario_ID;
	private int ERF_ID;
	private int Source_ID;
	private int Rupture_ID;
	
	public int getERF_ID() {
		return ERF_ID;
	}
	public void setERF_ID(int erf_id) {
		ERF_ID = erf_id;
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
	public int getSource_ID() {
		return Source_ID;
	}
	public void setSource_ID(int source_ID) {
		Source_ID = source_ID;
	}
	
	public boolean equals(Object obj) {
		if(this == obj)
			return true;
		if((obj == null) || (obj.getClass() != this.getClass()))
			return false;
		// object must be RuptureVariationPK at this point
		RuptureVariationPK papk = (RuptureVariationPK)obj;
		return (ERF_ID == papk.getERF_ID() &&
				Source_ID == papk.getSource_ID() &&
				Rupture_ID == papk.getRupture_ID() && 
				Rup_Var_ID == papk.getRup_Var_ID() &&
				Rup_Var_Scenario_ID == papk.getRup_Var_Scenario_ID());
	}
	
	public int hashCode() {
		return 0;
	}
	

}
