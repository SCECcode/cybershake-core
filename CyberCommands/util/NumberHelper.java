package util;

/*
 * Created on Mar 3, 2005
 *
 */


/**
 * @author wpun
 *
 */
public class NumberHelper {

	public static final int FLOAT_LENGTH = 4; 
	public static final int DOUBLE_LENGTH = 8;
	
	// bytes[0] is the least significant bytes
	public static float littleEndianByteToFloat(byte [] bytes) throws IncorrectNumberOfBytesException  {
		checkSufficientBytes(bytes, FLOAT_LENGTH); 
        int i = (bytes[0]&0xff) + ((bytes[1]&0xff) << 8) + ((bytes[2]&0xff) << 16) + (bytes[3] << 24);
        return Float.intBitsToFloat(i);
      }

	// bytes[0] is the most significant bytes
	public static float bigEndianByteToFloat(byte [] bytes) throws IncorrectNumberOfBytesException  {
		checkSufficientBytes(bytes, FLOAT_LENGTH); 
        int i = (bytes[3]&0xff) + ((bytes[2]&0xff) << 8) + ((bytes[1]&0xff) << 16) + (bytes[0] << 24);
//        int i = (bytes[3]) + ((bytes[2]) << 8) + ((bytes[1]) << 16) + (bytes[0] << 24);
		return Float.intBitsToFloat(i);
      }

	// bytes[0] is the least significant bytes
	public static double littleEndianByteToDouble(byte [] bytes) throws IncorrectNumberOfBytesException  {
		checkSufficientBytes(bytes, DOUBLE_LENGTH); 
		long value = 0;
        for (int i = 0; i < 8; i++) {
        	value += (long)((long)bytes[i] << (long)(i*8));
        }
        return Double.longBitsToDouble(value);
      }


	// bytes[0] is the most significant bytes
	public static double bigEndianByteToDouble(byte [] bytes) throws IncorrectNumberOfBytesException {
		checkSufficientBytes(bytes, DOUBLE_LENGTH); 
		long value = 0;
        for (int i = 0; i < 8; i++) {
        	value += (long)((long)bytes[7-i] << (long)(i*8));
        }
        return Double.longBitsToDouble(value);
      }
	private static void checkSufficientBytes(byte[] bytes, int typeLength) throws IncorrectNumberOfBytesException {
		int length = bytes.length;
		if (length  < typeLength) throw new IncorrectNumberOfBytesException(typeLength+" bytes required but "+length+" provided.");
	}
	public static byte[] swapBytes(byte[] bytes, int typeLength) throws IncorrectNumberOfBytesException {
		checkSufficientBytes(bytes, typeLength); 
		byte [] out = new byte[typeLength];
		for (int i = 0; i < typeLength; i++) {
			out[i] = bytes[(typeLength-1)-i];
		}
		return out; 
	}
}
