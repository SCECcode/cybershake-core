package mapping;

import java.io.Serializable;
import java.util.Date;

public class CyberShake_Runs implements Serializable {
	private static final long serialVersionUID = 1L;
	private int Run_ID;
	
	private CyberShakeSites cyberShakeSites;
	private ErfIds erfIds;
	private SgtVariationIds sgtVariationIds;
	private RuptureVariationScenarioIds ruptureVariationScenarioIds;
	
	private String SGT_Host;
	private String PP_Host;
	private Date SGT_Time;
	private Date PP_Time;
	private String Status;
	private Date Status_Time;
	private String Last_User;
	private String Comment;
	private String Job_ID;
	private String Submit_Dir;
	private String Notify_User;
	
	public boolean equals(Object obj)
	{
		if(this == obj)
			return true;
		if((obj == null) || (obj.getClass() != this.getClass()))
			return false;
		// object must be PeakAmplitudesPK at this point
		CyberShake_Runs csr = (CyberShake_Runs)obj;
		return (Run_ID == csr.getRun_ID() &&
				cyberShakeSites == csr.getCyberShakeSites() &&
				erfIds == csr.getErfIds() &&
				sgtVariationIds == csr.getSgtVariationIds()&&
				ruptureVariationScenarioIds == csr.getRuptureVariationScenarioIds() &&
				SGT_Host.equals(csr.getSGT_Host()) &&
				PP_Host.equals(csr.getPP_Host()) &&
				SGT_Time.equals(csr.getSGT_Time()) &&
				PP_Time.equals(csr.getPP_Time()) &&
				Status.equals(csr.getStatus()) &&
				Status_Time.equals(csr.getStatus_Time()) && 
				Last_User.equals(csr.getLast_User()) &&
				Comment.equals(csr.getComment()) &&
				Job_ID.equals(csr.getJob_ID()) &&
				Submit_Dir.equals(csr.getSubmit_Dir()) &&
				Notify_User.equals(csr.getNotify_User()) );
	}
	
	public int hashCode() {
		return 0;
	}
	
	public String toString() {
		return "Run_ID: " + Run_ID + ", Site_ID: " + cyberShakeSites.getCsSiteId() + ", ERF_ID: " + erfIds.getErfId() + ", SGT_Variation_ID: " + sgtVariationIds.getSgtVariationId() +
		", Rup_Var_Scenario_ID: " + ruptureVariationScenarioIds.getRupVarScenarioId() + ", SGT_Host: " + SGT_Host + ", PP_Host: " + PP_Host +
		", SGT_Time: " + SGT_Time.toString() + ", PP_Time: " + PP_Time.toString() + ", Status: " + Status +
		", Status_Time: " + Status_Time.toString() + ", Last_User: " + Last_User + ", Comment:" + Comment +
		", Job_ID: " + Job_ID + ", Submit_Dir: " + Submit_Dir + ", Notify_User: " + Notify_User;
	}

	public int getRun_ID() {
		return Run_ID;
	}

	public void setRun_ID(int run_ID) {
		Run_ID = run_ID;
	}

	public String getSGT_Host() {
		return SGT_Host;
	}

	public void setSGT_Host(String host) {
		SGT_Host = host;
	}

	public String getPP_Host() {
		return PP_Host;
	}

	public void setPP_Host(String host) {
		PP_Host = host;
	}

	public Date getSGT_Time() {
		return SGT_Time;
	}

	public void setSGT_Time(Date time) {
		SGT_Time = time;
	}

	public Date getPP_Time() {
		return PP_Time;
	}

	public void setPP_Time(Date time) {
		PP_Time = time;
	}

	public String getStatus() {
		return Status;
	}

	public void setStatus(String status) {
		Status = status;
	}

	public Date getStatus_Time() {
		return Status_Time;
	}

	public void setStatus_Time(Date status_Time) {
		Status_Time = status_Time;
	}

	public String getLast_User() {
		return Last_User;
	}

	public void setLast_User(String last_User) {
		Last_User = last_User;
	}

	public String getComment() {
		return Comment;
	}

	public void setComment(String comment) {
		Comment = comment;
	}

	public String getJob_ID() {
		return Job_ID;
	}

	public void setJob_ID(String job_ID) {
		Job_ID = job_ID;
	}

	public String getSubmit_Dir() {
		return Submit_Dir;
	}

	public void setSubmit_Dir(String submit_Dir) {
		Submit_Dir = submit_Dir;
	}

	public String getNotify_User() {
		return Notify_User;
	}

	public void setNotify_User(String notify_User) {
		Notify_User = notify_User;
	}

	public CyberShakeSites getCyberShakeSites() {
		return cyberShakeSites;
	}

	public void setCyberShakeSites(CyberShakeSites cyberShakeSites) {
		this.cyberShakeSites = cyberShakeSites;
	}

	public ErfIds getErfIds() {
		return erfIds;
	}

	public void setErfIds(ErfIds erfIds) {
		this.erfIds = erfIds;
	}

	public SgtVariationIds getSgtVariationIds() {
		return sgtVariationIds;
	}

	public void setSgtVariationIds(SgtVariationIds sgtVariationIds) {
		this.sgtVariationIds = sgtVariationIds;
	}

	public RuptureVariationScenarioIds getRuptureVariationScenarioIds() {
		return ruptureVariationScenarioIds;
	}

	public void setRuptureVariationScenarioIds(
			RuptureVariationScenarioIds ruptureVariationScenarioIds) {
		this.ruptureVariationScenarioIds = ruptureVariationScenarioIds;
	}

}

