package data;

public class SAPeriods {
	public double[] values = {10.0, 9.500286, 9.00009, 8.499788, 8.0, 7.500187, 6.99986, 6.499837, 5.99988, 5.499945, 5.0, 4.800076, 4.600028, 4.400053, 4.199916, 4.0, 3.79997, 3.599971, 3.399973, 3.2, 3.00003, 2.800022, 2.599969, 2.399981, 2.199978, 2.0, 0.01};
	private int currentValueIndex = 0;	
	
	public double getNextValue() {
		if (currentValueIndex == 27) 
			currentValueIndex = 0;
		currentValueIndex++;
		//System.out.println("SAPeriods.values[" + (currentValueIndex-1) + "]: " + values[currentValueIndex-1]);
		return values[currentValueIndex-1];
	}
}
