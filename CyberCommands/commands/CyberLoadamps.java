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

import data.RunID;

import processing.RuptureVariationFileInserter;

public class CyberLoadamps {

	private static final String NO_P_OPTION_MESSAGE = "Please use -p to set the path with the spectral acceleration files";
	private static final String NO_SERVER_OPTION_MESSAGE = "Please use -server to specify a database server";
	private static final String NO_RUNID_OPTION_MESSAGE = "Please use -run to specify the RunID";

	/**
	 * @param args
	 */
	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		Options options = new Options();

		Option sgt = OptionBuilder.withArgName("RunID").hasArg().withDescription("Run ID - this option is required").create("run");
		Option path = OptionBuilder.withArgName("directory").hasArg().withDescription("file path with spectral acceleration files, either top-level directory or zip file - this option is required").create("p");
		Option server = OptionBuilder.withArgName("name").hasArg().withDescription("server name (focal, surface or intensity) - this options is required").create("server");
        Option zip = new Option("z", "Read zip files instead of bsa.");
        Option help = new Option("help", "print this message");

		options.addOption(sgt);
		options.addOption(path);
		options.addOption(help);
		options.addOption(server);
        options.addOption(zip);
        
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

				if (!cmd.hasOption("server")) {
					System.out.println(NO_SERVER_OPTION_MESSAGE);
					return;
				}
				
                if (!cmd.hasOption("run")) {
                    System.out.println(NO_RUNID_OPTION_MESSAGE);
                    return;
                }

				System.out.println("Running loadamps using directory: " + cmd.getOptionValue("p") + " with Run ID: " + cmd.getOptionValue("run"));
				RunID rid = new RunID(Integer.parseInt(cmd.getOptionValue("run")));
				RuptureVariationFileInserter rvfi = new RuptureVariationFileInserter(cmd.getOptionValue("p"), rid.getSiteName(), rid.getSgtVarID(), cmd.getOptionValue("server"), rid.getRuptVarScenID(), rid.getErfID(), cmd.hasOption("z"));
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

	public static String getNO_RUNID_OPTION_MESSAGE() {
		return NO_RUNID_OPTION_MESSAGE;
	}

	public static String getNO_SERVER_OPTION_MESSAGE() {
		return NO_SERVER_OPTION_MESSAGE;
	}

}
