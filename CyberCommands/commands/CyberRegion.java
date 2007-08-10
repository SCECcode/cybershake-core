package commands;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class CyberRegion {
	
	private static final String NO_SITES_SPECIFIED_MESSAGE = "Please use one of the -all, -site, or -id options to specifiy one site or all sites";
	private static final String NO_SERVER_OPTION_MESSAGE = "Please use -server to specify a database server";

	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		Options options = new Options();

		Option all = new Option("all","list all sites and region for each site");
		Option id = OptionBuilder.withArgName("ID").hasArg().withDescription("give region for the specified site id").create("id");
		Option site = OptionBuilder.withArgName("name").hasArg().withDescription("gives region for the specified site name").create("site");
		Option server = OptionBuilder.withArgName("name").hasArg().withDescription("server name (surface or intensity) - this options is required").create("server");

		Option help = new Option("help", "print this message");
		
		options.addOption(all);
		options.addOption(id);
		options.addOption(site);
		options.addOption(server);
		options.addOption(help);
		
		CommandLineParser parser = new GnuParser();
		
		try {
			CommandLine cmd = parser.parse( options, args);

			if (cmd.hasOption("help")) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp( "CyberLoadamps", options, true );	
			}
			else {
				if (!cmd.hasOption("server")) {
					System.out.println(NO_SERVER_OPTION_MESSAGE);
					return;
				}
				
				if (!cmd.hasOption("all") && !cmd.hasOption("site") && !cmd.hasOption("id")) {
					System.out.println(NO_SITES_SPECIFIED_MESSAGE);
					return;
				}
			}
		} catch (ParseException e) {
			if (e instanceof  MissingArgumentException) {
				System.out.println(e.getMessage());
			}
			else {
				e.printStackTrace();
			}
		}
		
		
	}

	public static String getNO_SERVER_OPTION_MESSAGE() {
		return NO_SERVER_OPTION_MESSAGE;
	}

	public static String getNO_SITES_SPECIFIED_MESSAGE() {
		return NO_SITES_SPECIFIED_MESSAGE;
	}
}
