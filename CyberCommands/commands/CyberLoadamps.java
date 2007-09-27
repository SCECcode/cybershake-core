package commands;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.UnrecognizedOptionException;

import processing.RuptureVariationFileInserter;

public class CyberLoadamps {

	private static final String NO_SGT_OPTION_MESSAGE = "Please use -sgt to specify an SGT Variation ID";
	private static final String NO_SITE_OPTION_MESSAGE = "Please use -site to specify a site name";
	private static final String NO_P_OPTION_MESSAGE = "Please use -p to set the path with the spectral acceleration files";
	private static final String NO_SERVER_OPTION_MESSAGE = "Please use -server to specify a database server";
	//added by SC
    private static final String NO_RUP_VAR_ID_OPTION_MESSAGE = "Please use -rvid to specify a rupture variation ID";
    
	/**
	 * @param args
	 */
	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		Options options = new Options();

		Option sgt = OptionBuilder.withArgName("ID").hasArg().withDescription("SGT Variation ID - this option is required").create("sgt");
		Option p = OptionBuilder.withArgName("directory").hasArg().withDescription("file path with spectral acceleration files - this option is required").create("p");
		Option site = OptionBuilder.withArgName("name").hasArg().withDescription("site name for spectral acceleration files - this option is required").create("site");
		Option server = OptionBuilder.withArgName("name").hasArg().withDescription("server name (surface or intensity) - this options is required").create("server");
		//added by SC
        Option ruptureVariationID = OptionBuilder.withArgName("ruptureVariationID").hasArg().withDescription("Rupture Variation ID - this option is required").create("rvid");
        
		Option help = new Option("help", "print this message");

		
		options.addOption(site);
		options.addOption(sgt);
		options.addOption(p);
		options.addOption(help);
		options.addOption(server);
		//added by SC
        options.addOption(ruptureVariationID);

		CommandLineParser parser = new GnuParser();

		try {
			CommandLine cmd = parser.parse( options, args);

			if (cmd.hasOption("help")) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp( "CyberLoadAmps", options, true );	
			}
			else {
				if (!cmd.hasOption("p")) {
					System.out.println(NO_P_OPTION_MESSAGE);
					return;
				}

				if (!cmd.hasOption("site")) {
					System.out.println(NO_SITE_OPTION_MESSAGE);
					return;
				}
				
				if (!cmd.hasOption("sgt")) {
					System.out.println(NO_SGT_OPTION_MESSAGE);
					return;
				}
				
				if (!cmd.hasOption("server")) {
					System.out.println(NO_SERVER_OPTION_MESSAGE);
					return;
				}
                //added by SC
                if (!cmd.hasOption("rupVarID")) {
                    System.out.println(NO_RUP_VAR_ID_OPTION_MESSAGE);
                    return;
                }

				System.out.println("Running loadamps using directory: " + cmd.getOptionValue("p") + " as site: " + cmd.getOptionValue("site") + " and SGT Variation ID: " + cmd.getOptionValue("sgt"));
				RuptureVariationFileInserter rvfi = new RuptureVariationFileInserter(cmd.getOptionValue("p"), cmd.getOptionValue("site"), cmd.getOptionValue("sgt"), cmd.getOptionValue("server"), cmd.getOptionValue("rvid"));
				rvfi.performInsertions();
			}

		} catch (IOException e) {
			e.printStackTrace();
		} catch (UnrecognizedOptionException e) {
			System.out.println(e.getMessage());
			return;
		} catch (ParseException e) {
			if (e instanceof  MissingArgumentException) {
				System.out.println(e.getMessage());
			}
			else {
				e.printStackTrace();
			}
		}
	}

	public static String getNO_P_OPTION_MESSAGE() {
		return NO_P_OPTION_MESSAGE;
	}

	public static String getNO_SGT_OPTION_MESSAGE() {
		return NO_SGT_OPTION_MESSAGE;
	}

	public static String getNO_SITE_OPTION_MESSAGE() {
		return NO_SITE_OPTION_MESSAGE;
	}

	public static String getNO_SERVER_OPTION_MESSAGE() {
		return NO_SERVER_OPTION_MESSAGE;
	}

}
