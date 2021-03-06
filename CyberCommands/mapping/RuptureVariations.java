package mapping;

// Generated Aug 9, 2007 9:53:47 PM by Hibernate Tools 3.2.0.b9

import java.util.HashSet;
import java.util.Set;

/**
 * RuptureVariations generated by hbm2java
 */
public class RuptureVariations implements java.io.Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -7958570295601555096L;

	private RuptureVariationsId id;

	private RuptureVariationScenarioIds ruptureVariationScenarioIds;

	private ErfIds erfIds;

	private Ruptures ruptures;

	private String rupVarLfn;

	private Set<PeakAmplitudes> peakAmplitudeses = new HashSet<PeakAmplitudes>(
			0);

	public RuptureVariations() {
	}

	public RuptureVariations(RuptureVariationsId id,
			RuptureVariationScenarioIds ruptureVariationScenarioIds,
			ErfIds erfIds, Ruptures ruptures, String rupVarLfn) {
		this.id = id;
		this.ruptureVariationScenarioIds = ruptureVariationScenarioIds;
		this.erfIds = erfIds;
		this.ruptures = ruptures;
		this.rupVarLfn = rupVarLfn;
	}

	public RuptureVariations(RuptureVariationsId id,
			RuptureVariationScenarioIds ruptureVariationScenarioIds,
			ErfIds erfIds, Ruptures ruptures, String rupVarLfn,
			Set<PeakAmplitudes> peakAmplitudeses) {
		this.id = id;
		this.ruptureVariationScenarioIds = ruptureVariationScenarioIds;
		this.erfIds = erfIds;
		this.ruptures = ruptures;
		this.rupVarLfn = rupVarLfn;
		this.peakAmplitudeses = peakAmplitudeses;
	}

	public RuptureVariationsId getId() {
		return this.id;
	}

	public void setId(RuptureVariationsId id) {
		this.id = id;
	}

	public RuptureVariationScenarioIds getRuptureVariationScenarioIds() {
		return this.ruptureVariationScenarioIds;
	}

	public void setRuptureVariationScenarioIds(
			RuptureVariationScenarioIds ruptureVariationScenarioIds) {
		this.ruptureVariationScenarioIds = ruptureVariationScenarioIds;
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

	public String getRupVarLfn() {
		return this.rupVarLfn;
	}

	public void setRupVarLfn(String rupVarLfn) {
		this.rupVarLfn = rupVarLfn;
	}

	public Set<PeakAmplitudes> getPeakAmplitudeses() {
		return this.peakAmplitudeses;
	}

	public void setPeakAmplitudeses(Set<PeakAmplitudes> peakAmplitudeses) {
		this.peakAmplitudeses = peakAmplitudeses;
	}

}
