package data;

public class SAPeriods {
	public static double[] values = {10.0, 9.500286, 9.00009, 8.499788, 8.0, 7.500187, 6.99986, 6.499837, 5.99988, 5.499945, 5.0,
			4.800076, 4.600028, 4.400053, 4.199916, 4.0, 3.79997, 3.599971, 3.399973, 3.2, 3.00003, 2.800022, 2.599969,
			2.399981, 2.199978, 2.0, 1.666667, 1.428571, 1.25, 1.111111, 1.0, 0.6666667, 0.5, 0.4, 0.3333333, 0.2857143,
			0.25, 0.2222222, 0.2, 0.1666667, 0.1428571, 0.125, 0.1111111, 0.1, 0.01};
//	public static double[] values = {10.0, 9.500286, 9.00009, 8.499788, 8.0, 7.500187, 6.99986, 6.499837, 5.99988, 5.499945, 5.0,
//	4.800076, 4.600028, 4.400053, 4.199916, 4.0, 3.79997, 3.599971, 3.399973, 3.2, 3.00003, 2.800022, 2.599969,
//	2.399981, 2.199978, 2.0, 0.01};
	private int currentValueIndex = 0;
	public static int num_periods = values.length;
	
	public double getNextValue() {
		if (currentValueIndex == num_periods) 
			currentValueIndex = 0;
		currentValueIndex++;
		//System.out.println("SAPeriods.values[" + (currentValueIndex-1) + "]: " + values[currentValueIndex-1]);
		return values[currentValueIndex-1];
	}

	public double getCurrentValue() {
		return values[currentValueIndex-1];
	}
}
