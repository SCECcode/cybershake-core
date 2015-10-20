package data;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.Serializable;

public class DurationEntry implements Serializable {
	private static final long serialVersionUID = 1L;

	public int type;
	public int type_value;
	public int component;
	public float value;

	//Constants
	public static final int ARIAS_INTENSITY = 0;
	public static final int ENERGY_INTEGRAL = 1;
	public static final int CAV = 2;
	public static final int DV = 3;
	public static final int DA = 4;
	public static final int D5_75 = 5;
	public static final int D5_95 = 6;
	public static final int D20_80 = 7;
	
	public static final int X_COMP = 0;
	public static final int Y_COMP = 1;
	
	public DurationEntry() {
	}

	public int populate(DataInputStream stream) throws IOException {
		type = stream.readInt();
		type_value = stream.readInt();
		component = stream.readInt();
		value = stream.readFloat();
		return 0;
	}
	
	public static String getKeyString(String databaseString, String component) {
		StringBuffer keyString = new StringBuffer("");
		if (databaseString.contains("Arias intensity")) {
			keyString.append(ARIAS_INTENSITY);
		} else if (databaseString.contains("energy integral")) {
			keyString.append(ENERGY_INTEGRAL);
		} else if (databaseString.contains("Cumulative absolute velocity")) {
			keyString.append(CAV);
		} else if (databaseString.contains("velocity")) {
			keyString.append(DV);
		} else if (databaseString.contains("acceleration")) {
			keyString.append(DA);
		}
		if (databaseString.contains("5% to 75%")) {
			keyString.append("_" + D5_75);
		} else if (databaseString.contains("5% to 95%")) {
			keyString.append("_" + D5_95);
		} else if (databaseString.contains("20% to 80%")) {
			keyString.append("_" + D20_80);
		}
		if (component.equals("X")) {
			keyString.append("_" + X_COMP);
		} else if (component.equals("Y")) {
			keyString.append("_" + Y_COMP);
		}
		return keyString.toString();
	}
}
