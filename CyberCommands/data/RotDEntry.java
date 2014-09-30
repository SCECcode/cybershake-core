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
	
	public int parse(FileInputStream stream) throws IOException {
		byte[] byteArray = new byte[16];
		stream.read(byteArray);
		byte[] outByteArray = SwapBytes.swapByteArrayToByteArrayForFloats(byteArray);
		
		DataInputStream dr = new DataInputStream(new ByteArrayInputStream(outByteArray));
		period = dr.readFloat();
		rotD100 = dr.readFloat();
		rotD100_angle = dr.readInt();
		rotD50 = dr.readFloat();
		
		dr.close();
		return 0;
	}
	
	
}
