package util;

import java.io.*;
import java.util.ArrayList;
import java.util.StringTokenizer;

/*
 * Created on Aug 26, 2005
 *
 */

/**
 * @author francoeu
 *	Simple Output utility to go from Davids binary format to space deliminated floats
 */
public class SABinary2Float {
	
	public static void convertFileToFile(String file) throws FileNotFoundException{
		try {
			DataInputStream dr = new DataInputStream(new FileInputStream(file));
			BufferedWriter out = new BufferedWriter(new FileWriter(getOutputFileName(file)));
			float val;
			int counter=0;
			while(dr.available()>=4){
				val=dr.readFloat();
				System.out.print(val+ " ");
				out.write(Float.toString(val) + "\n");
				counter++;
			}
			System.out.println();
			System.out.println("Total: "+counter);
			dr.close();
			out.flush();
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static ArrayList<Float> convertFileToArrayOfFloats(String file) throws FileNotFoundException {
		ArrayList<Float> floats = new ArrayList<Float>();
		try {
			DataInputStream dr = new DataInputStream(new FileInputStream(file));
			float val;
			while(dr.available()>=4){
				val=dr.readFloat();
				floats.add(val);
			}
			dr.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return floats;
	}
	
	public static ArrayList<Float> convertByteArrayToArrayOfFloats(byte[] byteArray) {
		ArrayList<Float> floats = new ArrayList<Float>();
		try {
			DataInputStream dr = new DataInputStream(new ByteArrayInputStream(byteArray));
			float val;
			while(dr.available()>=4){
				val=dr.readFloat();
				floats.add(val);
			}
			dr.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return floats;
	}
	
	public static ArrayList<Integer> convertByteArrayToArrayOfInts(byte[] byteArray) {
		ArrayList<Integer> ints = new ArrayList<Integer>();
		try {
			DataInputStream dr = new DataInputStream(new ByteArrayInputStream(byteArray));
			int val;
			while(dr.available()>=4){
				val=dr.readInt();
				ints.add(val);
			}
			dr.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return ints;
	}
	
	public static String getOutputFileName(String infile){
		String outname = new String();
		String text = null;
		StringTokenizer st = new StringTokenizer(infile,"/");
 	    while (st.hasMoreTokens()) {
 	         text=st.nextToken();
 	    }
 	    text = text.replaceAll(".bsa","");
 	    outname= text+ "_floats.txt";
		return(outname);
	}
}

