package mapping;

// Generated Aug 9, 2007 9:53:47 PM by Hibernate Tools 3.2.0.b9

import java.util.HashSet;
import java.util.Set;

/**
 * RuptureVariationScenarioIds generated by hbm2java
 */
public class RuptureVariationScenarioIds implements java.io.Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 8113619046795108510L;

	private int rupVarScenarioId;

	private ErfIds erfIds;

	private String rupVarScenarioName;

	private String rupVarScenarioDescription;

	private Set<RuptureVariationScenarioMetadata> ruptureVariationScenarioMetadatas = new HashSet<RuptureVariationScenarioMetadata>(
			0);

	private Set<PeakAmplitudes> peakAmplitudeses = new HashSet<PeakAmplitudes>(
			0);

	private Set<RuptureVariations> ruptureVariationses = new HashSet<RuptureVariations>(
			0);

	public RuptureVariationScenarioIds() {
	}

	public RuptureVariationScenarioIds(int rupVarScenarioId, ErfIds erfIds,
			String rupVarScenarioName) {
		this.rupVarScenarioId = rupVarScenarioId;
		this.erfIds = erfIds;
		this.rupVarScenarioName = rupVarScenarioName;
	}

	public RuptureVariationScenarioIds(
			int rupVarScenarioId,
			ErfIds erfIds,
			String rupVarScenarioName,
			String rupVarScenarioDescription,
			Set<RuptureVariationScenarioMetadata> ruptureVariationScenarioMetadatas,
			Set<PeakAmplitudes> peakAmplitudeses,
			Set<RuptureVariations> ruptureVariationses) {
		this.rupVarScenarioId = rupVarScenarioId;
		this.erfIds = erfIds;
		this.rupVarScenarioName = rupVarScenarioName;
		this.rupVarScenarioDescription = rupVarScenarioDescription;
		this.ruptureVariationScenarioMetadatas = ruptureVariationScenarioMetadatas;
		this.peakAmplitudeses = peakAmplitudeses;
		this.ruptureVariationses = ruptureVariationses;
	}

	public int getRupVarScenarioId() {
		return this.rupVarScenarioId;
	}

	public void setRupVarScenarioId(int rupVarScenarioId) {
		this.rupVarScenarioId = rupVarScenarioId;
	}

	public ErfIds getErfIds() {
		return this.erfIds;
	}

	public void setErfIds(ErfIds erfIds) {
		this.erfIds = erfIds;
	}

	public String getRupVarScenarioName() {
		return this.rupVarScenarioName;
	}

	public void setRupVarScenarioName(String rupVarScenarioName) {
		this.rupVarScenarioName = rupVarScenarioName;
	}

	public String getRupVarScenarioDescription() {
		return this.rupVarScenarioDescription;
	}

	public void setRupVarScenarioDescription(String rupVarScenarioDescription) {
		this.rupVarScenarioDescription = rupVarScenarioDescription;
	}

	public Set<RuptureVariationScenarioMetadata> getRuptureVariationScenarioMetadatas() {
		return this.ruptureVariationScenarioMetadatas;
	}

	public void setRuptureVariationScenarioMetadatas(
			Set<RuptureVariationScenarioMetadata> ruptureVariationScenarioMetadatas) {
		this.ruptureVariationScenarioMetadatas = ruptureVariationScenarioMetadatas;
	}

	public Set<PeakAmplitudes> getPeakAmplitudeses() {
		return this.peakAmplitudeses;
	}

	public void setPeakAmplitudeses(Set<PeakAmplitudes> peakAmplitudeses) {
		this.peakAmplitudeses = peakAmplitudeses;
	}

	public Set<RuptureVariations> getRuptureVariationses() {
		return this.ruptureVariationses;
	}

	public void setRuptureVariationses(
			Set<RuptureVariations> ruptureVariationses) {
		this.ruptureVariationses = ruptureVariationses;
	}

}
