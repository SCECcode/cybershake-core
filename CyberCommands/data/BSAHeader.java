package data;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;

import util.IncorrectNumberOfBytesException;
import util.NumberHelper;
import util.SwapBytes;

public class BSAHeader {
	private char[] version = new char[8];
	public String versionString;
	private char[] site_name = new char[8];
	public String siteString;
	private char[] padding = new char[8];
	public int source_id;
	public int rupture_id;
	public int rup_var_id;
	public float dt;
	public int nt;
	private int comps;
	public float det_max_freq;
	public float stoch_max_freq;
	
	public int num_comps = 0;
	public boolean x_comp = false;
	public boolean y_comp = false;
	public boolean z_comp = false;
	
	private final int X_COMP_FLAG = 1;
	private final int Y_COMP_FLAG = 2;
	private final int Z_COMP_FLAG = 4;
	
	public BSAHeader() {}
	
	public int parse(ByteArrayInputStream bStream) throws IOException {
		num_comps = 0;
		x_comp = false;
		y_comp = false;	
		z_comp = false;
		int i, ret;
		byte[] byteArray = new byte[version.length];
		ret = bStream.read(byteArray);
		if (ret==-1) { //EOF 
			return -1;
		}
			
		StringBuffer sb = new StringBuffer();
		//This will stop short if we find a null terminator
		for (i=0; i<version.length; i++) {
			version[i] = (char) byteArray[i];
			if (version[i]!='\0') {
				sb.append(version[i]);
			} else {
				break;
			}
		}
		versionString = sb.toString();
		if (!versionString.equals("12.10")) {
			System.err.println("Version of this header is " + versionString + ", don't know how to parse it.");
			System.exit(-5);
		}
		
		byteArray = new byte[site_name.length];
		bStream.read(byteArray);
		sb = new StringBuffer();
		for (i=0; i<site_name.length; i++) {
			site_name[i] = (char) byteArray[i];
			if (site_name[i]!='\0') {
				sb.append(site_name[i]);
			} else {
				break;
			}
		}
		siteString = sb.toString();
		
		bStream.skip((long)(padding.length));
		//read the remaining bytes into an array, swap them, then read them out with a stream
		byteArray = new byte[32];
		bStream.read(byteArray);
		byte[] outByteArray = SwapBytes.swapByteArrayToByteArrayForFloats(byteArray);
		
		DataInputStream dr = new DataInputStream(new ByteArrayInputStream(outByteArray));
		source_id = dr.readInt();
		rupture_id = dr.readInt();
		rup_var_id = dr.readInt();
		dt = dr.readFloat();
		nt = dr.readInt();
		comps = dr.readInt();
		det_max_freq = dr.readFloat();
		stoch_max_freq = dr.readFloat();
		
		if ((comps & X_COMP_FLAG)!=0) {
			x_comp = true;
			num_comps++;
		}
		if ((comps & Y_COMP_FLAG)!=0) {
			y_comp = true;
			num_comps++;
		}
		if ((comps & Z_COMP_FLAG)!=0) {
			z_comp = true;
			num_comps++;
		}
		
		dr.close();
		
		return 0;
	}
	
	public int parse(FileInputStream stream) throws IOException {
		num_comps = 0;
		x_comp = false;
		y_comp = false;	
		z_comp = false;
		int i, ret;
		byte[] byteArray = new byte[version.length];
		ret = stream.read(byteArray);
		if (ret==-1) { //EOF 
			return -1;
		}
			
		StringBuffer sb = new StringBuffer();
		//This will stop short if we find a null terminator
		for (i=0; i<version.length; i++) {
			version[i] = (char) byteArray[i];
			if (version[i]!='\0') {
				sb.append(version[i]);
			} else {
				break;
			}
		}
		versionString = sb.toString();
		if (!versionString.equals("12.10")) {
			System.err.println("Version of this header is " + versionString + ", don't know how to parse it.");
			System.err.println("Stream location: " + stream.getChannel().position());
			System.exit(-5);
		}
		
		byteArray = new byte[site_name.length];
		stream.read(byteArray);
		sb = new StringBuffer();
		for (i=0; i<site_name.length; i++) {
			site_name[i] = (char) byteArray[i];
			if (site_name[i]!='\0') {
				sb.append(site_name[i]);
			} else {
				break;
			}
		}
		siteString = sb.toString();
		
		stream.skip((long)(padding.length));
		//read the remaining bytes into an array, swap them, then read them out with a stream
		byteArray = new byte[32];
		stream.read(byteArray);
		byte[] outByteArray = SwapBytes.swapByteArrayToByteArrayForFloats(byteArray);
		
		DataInputStream dr = new DataInputStream(new ByteArrayInputStream(outByteArray));
		source_id = dr.readInt();
		rupture_id = dr.readInt();
		rup_var_id = dr.readInt();
		dt = dr.readFloat();
		nt = dr.readInt();
		comps = dr.readInt();
		det_max_freq = dr.readFloat();
		stoch_max_freq = dr.readFloat();
		
		if ((comps & X_COMP_FLAG)!=0) {
			x_comp = true;
			num_comps++;
		}
		if ((comps & Y_COMP_FLAG)!=0) {
			y_comp = true;
			num_comps++;
		}
		if ((comps & Z_COMP_FLAG)!=0) {
			z_comp = true;
			num_comps++;
		}
		
		dr.close();
		
		return 0;
	}
	
	public void print() {
		System.out.println("version: " + versionString);
		System.out.println("site: " + siteString);
		System.out.println("source_id: " + source_id);
		System.out.println("rupture_id: " + rupture_id);
		System.out.println("rup_var_id: " + rup_var_id);
		System.out.println("dt: " + dt);
		System.out.println("nt: " + nt);
		System.out.println("num_comps: " + num_comps);
		System.out.println("det_max_freq: " + det_max_freq);
		System.out.println("stoch_max_freq: " + stoch_max_freq);
	}
}
