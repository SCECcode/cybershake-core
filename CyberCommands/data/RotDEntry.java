package data;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.Serializable;

import util.SwapBytes;

public class RotDEntry implements Serializable {
	private static final long serialVersionUID = 1L;
	
	public float period;
	public float rotD100;
	public int rotD100_angle;
	public float rotD50;
	
	public RotDEntry() {
	}
	
	public int populate(DataInputStream stream) throws IOException {
		period = stream.readFloat();
		rotD100 = stream.readFloat();
		rotD100_angle = stream.readInt();
		rotD50 = stream.readFloat();
		return 0;
	}
	
	
}
