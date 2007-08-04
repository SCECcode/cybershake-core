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

		Option sgt = OptionBuilder.withArgName("ID").hasArg().withDescription("SGT Variation ID - this option is required").create("sgt");
		Option p = OptionBuilder.withArgName("directory").hasArg().withDescription("file path with spectral acceleration files - this option is required").create("p");
		Option site = OptionBuilder.withArgName("name").hasArg().withDescription("site name for spectral acceleration files - this option is required").create("site");
		Option server = OptionBuilder.withArgName("name").hasArg().withDescription("server name (surface or intensity) - this options is required").create("server");

		Option help = new Option("help", "print this message");

		
		options.addOption(site);
		options.addOption(sgt);
		options.addOption(p);
		options.addOption(help);
		options.addOption(server);
		

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
				
				if (!cmd.hasOption("sgt")) {
					System.out.println("Please use -sgt to specify an SGT Variation ID");
					return;
				}

				System.out.println("Running loadamps using directory: " + cmd.getOptionValue("p") + " as site: " + cmd.getOptionValue("site") + " and SGT Variation ID: " + cmd.getOptionValue("sgt"));
				RuptureVariationFileInserter rvfi = new RuptureVariationFileInserter(cmd.getOptionValue("p"), cmd.getOptionValue("site"), cmd.getOptionValue("sgt"), cmd.getOptionValue("server"));
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
