package commands;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.UnrecognizedOptionException;

import data.RunID;

import processing.RuptureVariationFileInserter;

public class CyberLoadamps {

	private static final String NO_P_OPTION_MESSAGE = "Please use -p to set the path with the spectral acceleration files";
	private static final String NO_SERVER_OPTION_MESSAGE = "Please use -server to specify a database server";
	private static final String NO_RUNID_OPTION_MESSAGE = "Please use -run to specify the RunID";
	private static final String NO_PERIODS_OPTION_MESSAGE = "Please use -periods to specify the periods to insert";

	public enum Mode {BSA, ZIP, HEAD, ROTD, DURATION};
	/**
	 * @param args
	 */
	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		Options options = new Options();

		Option sgt = OptionBuilder.withArgName("RunID").hasArg().withDescription("Run ID - this option is required").create("run");
		Option path = OptionBuilder.withArgName("directory").hasArg().withDescription("file path with spectral acceleration files, either top-level directory or zip file - this option is required").create("p");
		Option server = OptionBuilder.withArgName("name").hasArg().withDescription("server name (focal, surface, intensity, moment, or csep-x), or sqlite:<filename> - this option is required").create("server");
        Option zip = new Option("z", "Read zip files instead of bsa.");
        Option header = new Option("d", "Assume one BSA file per rupture, with embedded header information.");
        Option valuesToInsert = OptionBuilder.withArgName("insertion_values").hasArg().withDescription("Which values to insert -\ngm:\tgeometric mean PSA data (default)\nxy:\tX and Y component PSA data\ngmxy:  Geometric mean and X and Y components").create("i");
        Option periods = OptionBuilder.withArgName("periods").hasArg().withDescription("Comma-delimited periods to insert").create("periods");
        Option rotd = new Option("r", "Read rotd files (instead of bsa.)");
        Option duration = new Option("u", "Read duration files (instead of bsa.)");
        Option convert = new Option("c", "Convert values from g to cm/sec^2");
        Option force = new Option("f", "Don't apply value checks to insertion values; use with care!.");
        Option verboseOption = new Option("v", "Print more messages about what is going on.");
        Option help = new Option("help", "print this message");

        OptionGroup fileGroup = new OptionGroup();
        fileGroup.addOption(zip);
        fileGroup.addOption(header);
        fileGroup.addOption(rotd);
        fileGroup.addOption(duration);
        
		options.addOption(sgt);
		options.addOption(path);
		options.addOption(help);
		options.addOption(server);
        options.addOption(valuesToInsert);
        options.addOption(periods);
        options.addOption(convert);
        options.addOption(force);
        options.addOption(verboseOption);
        options.addOptionGroup(fileGroup);
        
		CommandLineParser parser = new GnuParser();

		boolean forceInsert = false;
		boolean verbose = false;
		
		try {
			if (args.length==0) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp( "CyberLoadAmps", options, true );
				return;
			}
			
			CommandLine cmd = parser.parse( options, args);

			if (cmd.hasOption("help")) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp( "CyberLoadAmps", options, true );	
			}
			else {
				if (!cmd.hasOption("p")) {
					System.out.println(NO_P_OPTION_MESSAGE);
					System.exit(-1);
				}

				if (!cmd.hasOption("server")) {
					System.out.println(NO_SERVER_OPTION_MESSAGE);
					System.exit(-1);
				}
				
                if (!cmd.hasOption("run")) {
                    System.out.println(NO_RUNID_OPTION_MESSAGE);
                    System.exit(-1);
                }
                
                if (!cmd.hasOption("periods")) {
                	System.out.println(NO_PERIODS_OPTION_MESSAGE);
                	System.exit(-1);
                }

                if (cmd.hasOption("f")) {
                	forceInsert = true;
                }
                
                if (cmd.hasOption("v")) {
                	verbose = true;
                }
                
				System.out.println("Running loadamps using directory: " + cmd.getOptionValue("p") + " with Run ID: " + cmd.getOptionValue("run"));
				RunID rid = new RunID(Integer.parseInt(cmd.getOptionValue("run")), cmd.getOptionValue("server"));
//				RuptureVariationFileInserter rvfi = new RuptureVariationFileInserter(cmd.getOptionValue("p"), rid.getSiteName(), rid.getSgtVarID(), cmd.getOptionValue("server"), rid.getRuptVarScenID(), rid.getErfID(), cmd.hasOption("z"));
				String insertValues;
				if (!cmd.hasOption("i")) {
					//Geometric mean is default
					insertValues = "gm";
				} else {
					insertValues = cmd.getOptionValue("i");
				}
					
				Mode m = Mode.BSA;
				
				if (cmd.hasOption("z")) {
					m = Mode.ZIP;
				} else if (cmd.hasOption("d")) {
					m = Mode.HEAD;
				} else if (cmd.hasOption("r")) {
					m = Mode.ROTD;
				} else if (cmd.hasOption("u")) {
					m = Mode.DURATION;
				}
				
				boolean convertGtoCM = false;
				if (cmd.hasOption("c")) {
					convertGtoCM = true;
				}
				
				RuptureVariationFileInserter rvfi = new RuptureVariationFileInserter(cmd.getOptionValue("p"), rid, cmd.getOptionValue("server"), m, cmd.getOptionValue("periods"), insertValues, convertGtoCM, forceInsert, verbose);
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
