package test.fitnesse;

import fit.ColumnFixture;

public class TwoItemedFixture extends ColumnFixture {
	public TwoItemed twoItemed;
	
	public TwoItemed multiply() {
		twoItemed.multiply();
		return twoItemed;
	}
}
