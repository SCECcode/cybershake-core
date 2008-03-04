package data;

public class DirectionalComponent {
	public static final int periodsLength = 27;
	public Float[] periods = new Float[periodsLength];
	
	
	public DirectionalComponent() {
		for (int i=0;i<periods.length;i++) {
			periods[i] = new Float(0);
		}
	}
	
	public String toString() {
		String periodsString = new String();
		
		for (int i=0; i<periods.length; i++) {
			periodsString += periods[i].toString();
			if (i+1<periods.length) {
				periodsString += " ";
			}
		}
		
		return periodsString;
	}
	
	public static DirectionalComponent parse(String s) {
		String[] result = s.split(" ");
		DirectionalComponent dc = new DirectionalComponent();
		if (result.length == periodsLength) {
			for (int i=0;i<result.length;i++) {
				dc.periods[i] = Float.parseFloat(result[i]);
			}			
		}
		return dc;

	}
	
	public boolean equals (Object obj) {
		if (obj == this) {
			return true;
		}
		if ((obj == null) || (this.getClass() != obj.getClass() ) ) {
			return false;
		}
		DirectionalComponent dc = (DirectionalComponent) obj;
		for (int i=0;i<periods.length;i++) {
			if (!(dc.periods[i].equals(this.periods[i])) ) {
				return false;
			}
		}
		return true;
	}
	
	public int hashCode() {
		return 0;
	}
}
