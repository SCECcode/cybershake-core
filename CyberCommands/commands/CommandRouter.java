package commands;

public class CommandRouter {
	public static void main(String[] args) {
		if (args[0].equals("CyberLoadamps")) {
			String[] cmdArgs = createCmdArgs(args);
			CyberLoadamps.main(cmdArgs);
		}
		else if (args[0].equals("CyberRegion")) {
			String[] cmdArgs = createCmdArgs(args);
			CyberRegion.main(cmdArgs);
		}
	}

	public static String[] createCmdArgs(String[] args) {
		String[] cmdArgs = new String[args.length-1];
		int j=0;
		for (int i=1;i<args.length;i++) {
			cmdArgs[j] = args[i];
			j++;
		}
		return cmdArgs;
	}
}
