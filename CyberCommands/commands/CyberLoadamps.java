package commands;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.UnrecognizedOptionException;

import processing.RuptureVariationFileInserter;

public class CyberLoadamps {

	/**
	 * @param args
	 */
	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		Options options = new Options();

		Option p = OptionBuilder.withArgName("directory").hasArg().withDescription("file path with spectral acceleration files").create("p");
		Option site = OptionBuilder.withArgName("name").hasArg().withDescription("site name for spectral acceleration files").create("site");

		Option help = new Option("help", "print this message");

		options.addOption(p);
		options.addOption(site);
		options.addOption(help);

		CommandLineParser parser = new GnuParser();

		try {
			CommandLine cmd = parser.parse( options, args);

			if (cmd.hasOption("help")) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp( "CyberLoadamps", options );	
			}
			else {
				if (!cmd.hasOption("p")) {
					System.out.println("Please use -p to set the path with the spectral acceleration files");
					return;
				}

				if (!cmd.hasOption("site")) {
					System.out.println("Please use -site to specify a site name");
					return;
				}

				System.out.println("Running loadamps using directory: " + cmd.getOptionValue("p") + " and site: " + cmd.getOptionValue("site"));
				RuptureVariationFileInserter rvfi = new RuptureVariationFileInserter(cmd.getOptionValue("p"), cmd.getOptionValue("site"));
				rvfi.performInsertions();
			}

		} catch (IOException e) {
			e.printStackTrace();
		} catch (UnrecognizedOptionException e) {
			System.out.println(e.getMessage());
			return;
		} catch (ParseException e) {
			e.printStackTrace();
		}
		




	}

}
