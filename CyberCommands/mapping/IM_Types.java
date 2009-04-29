package mapping;

import java.io.Serializable;

public class IM_Types implements Serializable {
	private static final long serialVersionUID = 1L;
	private int IM_Type_ID;
	private String IM_Type_Measure;
	private double IM_Type_Value;
	private String Units;

	public boolean equals(Object obj)
	{
		if(this == obj)
			return true;
		if((obj == null) || (obj.getClass() != this.getClass()))
			return false;

		IM_Types imt = (IM_Types)obj;
		return (IM_Type_ID == imt.getIM_Type_ID() &&
				IM_Type_Measure.equals(imt.getIM_Type_Measure()) && 
				IM_Type_Value == imt.getIM_Type_Value() &&
				Units.equals(imt.getUnits()) );
	}
	
	public int hashCode() {
		return 0;
	}
	
	public String toString() {
		return "IM_Type_ID: " + IM_Type_ID + ", IM_Type_Measure: " + IM_Type_Measure + ", IM_Type_Value: " + IM_Type_Value + ", Units: " + Units;
	}

	public int getIM_Type_ID() {
		return IM_Type_ID;
	}

	public void setIM_Type_ID(int type_ID) {
		IM_Type_ID = type_ID;
	}

	public String getIM_Type_Measure() {
		return IM_Type_Measure;
	}

	public void setIM_Type_Measure(String type_Measure) {
		IM_Type_Measure = type_Measure;
	}

	public double getIM_Type_Value() {
		return IM_Type_Value;
	}

	public void setIM_Type_Value(double type_Value) {
		IM_Type_Value = type_Value;
	}

	public String getUnits() {
		return Units;
	}

	public void setUnits(String units) {
		Units = units;
	}
}
