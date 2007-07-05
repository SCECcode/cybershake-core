package test.commands;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import processing.RuptureVariationFileInserter;

public class Cyber {

	/**
	 * @param args
	 */
	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		Options options = new Options();
		
		Option loadamps = new Option("loadamps",false,"load spectral acceleration files (requires the -p and -site options)");
		Option p = OptionBuilder.withArgName("directory").hasArg().withDescription("file path with spectral acceleration files").create("p");
		Option site = OptionBuilder.withArgName("name").hasArg().withDescription("site name for spectral acceleration files").create("site");
		
		Option help = new Option("help", "print this message");
		
		options.addOption(loadamps);
		options.addOption(p);
		options.addOption(site);
		options.addOption(help);
		
		CommandLineParser parser = new GnuParser();
		
		try {
			CommandLine cmd = parser.parse( options, args);
			if(cmd.hasOption("loadamps")) {
				
				if (!cmd.hasOption("p")) {
					System.out.println("Please use -p to set the path with the spectral acceleration files");
					return;
				}
				
				if (!cmd.hasOption("site")) {
					System.out.println("Please use -site to specify a site name");
					return;
				}
				
				System.out.println("Got loadamps with directory: " + cmd.getOptionValue("p") + " and site: " + cmd.getOptionValue("site"));
				RuptureVariationFileInserter rvfi = new RuptureVariationFileInserter(cmd.getOptionValue("p"), cmd.getOptionValue("site"));
			}
			
			if (cmd.hasOption("help")) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp( "Cyber", options );	
			}
			else {
			}
			
		} catch (ParseException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		

		

	}

}
