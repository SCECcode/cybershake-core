package test.fitnesse;

public class TwoItemed {
	public int firstInt;
	public int secondInt;
	
	public TwoItemed(int newFirstInt, int newSecondInt) {
		firstInt = newFirstInt;
		secondInt = newSecondInt;
	}

	public TwoItemed() {
		firstInt = 0;
		secondInt = 0;
	}

	public void multiply() {
		int origFirstInt = firstInt;
		firstInt *= secondInt;
		secondInt *= origFirstInt;
	}
	
	public String toString() {
		return Integer.toString(firstInt) + " " + Integer.toString(secondInt);
	}
	
	public static TwoItemed parse(String s) {
		String[] result = s.split(" ");
		TwoItemed ti = new TwoItemed();
		ti.firstInt = Integer.parseInt(result[0]);
		ti.secondInt = Integer.parseInt(result[1]);
		return ti;
	}
	
	public boolean equals(Object obj) {
		if (obj == this) {
			return true;
		}
		if ((obj == null) || (this.getClass() != obj.getClass() ) ) {
			return false;
		}
		TwoItemed ti = (TwoItemed) obj;
		return (ti.firstInt == this.firstInt && ti.secondInt == this.secondInt);
		
		
	}
	
	public int hashCode() {
		return 0;
	}
}
