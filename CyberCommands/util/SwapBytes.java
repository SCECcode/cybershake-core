package util;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class SwapBytes {

	public static void swapFileToFile(String filename) {
		ArrayList<Float> floats = new ArrayList<Float>();
		read(filename, floats);
		commonSwap(floats);
		String outputFilename = constructOutputFilename(filename);
		writeDoubles(outputFilename, floats);
	}

	public static byte[] swapBinaryFileToArrayForFloats(String filename) {
		File inputFile = new File(filename);
		long fileSize = inputFile.length();
		ArrayList<byte[]> byteArray = new ArrayList<byte[]>();
		byte[] returnBytes = new byte[(int)fileSize];
		byte[] tempArray = new byte[NumberHelper.FLOAT_LENGTH];
		int counter = 0;
		try {
			FileInputStream input = new FileInputStream(inputFile);
			while (counter<fileSize) {
				input.read(tempArray);
				byteArray.add(NumberHelper.swapBytes(tempArray, NumberHelper.FLOAT_LENGTH));
				counter+=NumberHelper.FLOAT_LENGTH;
			}
			input.close();
			for (int i=0; i<byteArray.size(); i++) {
				for (int j=0; j<NumberHelper.FLOAT_LENGTH; j++) {
					returnBytes[NumberHelper.FLOAT_LENGTH*i + j] = byteArray.get(i)[j];
				}
			}
			return returnBytes;
		} catch (FileNotFoundException e) {
			System.err.println("File " + inputFile + " is missing, continuing.");
		} catch (IOException e) {
			e.printStackTrace();
		} catch (IncorrectNumberOfBytesException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	private static void commonSwap(ArrayList<Float> floats) {
		ArrayList<byte[]> byteArray = new ArrayList<byte[]>();
		makeByteArray(floats, byteArray);
		ArrayList<byte[]> swappedByteArray = new ArrayList<byte[]>();
		swapBytes(byteArray, swappedByteArray);
		floats.clear();
		makeDoubles(swappedByteArray, floats);
	}
	
	private static String constructOutputFilename(String filename) {
		return "swapped_" + filename;
	}

	public static void swapArrayOfFloatsToFile(ArrayList<Float> floats, String outputFilename) {
		commonSwap(floats);
		writeDoubles(outputFilename, floats);
	}
	
	public static ArrayList<Float> swapFileToArrayOfFloats(String filename) {
		ArrayList<Float> floats = new ArrayList<Float>();
		read(filename, floats);
		commonSwap(floats);
		return floats;
	}
	
	public static ArrayList<Float> swapArrayOfFloatsToArrayOfFloats(ArrayList<Float> floats) {
		commonSwap(floats);
		return floats;
	}
	
	private static void writeDoubles(String outputFilename, ArrayList<Float> floats) {
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilename)));
			for (float f: floats) {
				pw.println(f);
			}
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	private static void makeDoubles(ArrayList<byte[]> swappedByteArray, ArrayList<Float> floats) {
		ByteArrayInputStream byteIn;
		DataInputStream dataIn;
		for (byte[] b: swappedByteArray) {
			dataIn = new DataInputStream(byteIn = new ByteArrayInputStream(b));
			try {
				floats.add(dataIn.readFloat());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	private static void swapBytes(ArrayList<byte[]> byteArray, ArrayList<byte[]> swappedByteArray) {
		for (byte[] b: byteArray) {
			try {
				swappedByteArray.add(NumberHelper.swapBytes(b, NumberHelper.FLOAT_LENGTH));
			} catch (IncorrectNumberOfBytesException e) {
				e.printStackTrace();
			}
		}
	}

	private static void makeByteArray(ArrayList<Float> floats, ArrayList<byte[]> byteArray) {
		ByteArrayOutputStream byteOut;
		DataOutputStream dataOut = new DataOutputStream(byteOut = new ByteArrayOutputStream());
		try {
			for (float f : floats) {
				dataOut.writeFloat(f);
				byteArray.add(byteOut.toByteArray());
				byteOut.reset();
			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	private static void read(String filename, ArrayList<Float> floats) {
		File f = new File(filename);
		try {
			BufferedReader br = new BufferedReader(new FileReader(f));
			String line = br.readLine();
			while (line!=null) {
				floats.add(Float.parseFloat(line));
				line = br.readLine();
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}

