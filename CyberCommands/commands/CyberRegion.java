package commands;

import java.util.List;

import mapping.CyberShakeSiteRegions;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;
import org.hibernate.criterion.Restrictions;

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
				formatter.printHelp( "CyberRegion", options, true );	
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


				peformOperations(cmd);
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

	private static void peformOperations(CommandLine cmd) {
		String intensityConfigFilename = "intensity.cfg.xml";
		String surfaceConfigFilename = "surface.cfg.xml";
		String defaultConfigFilename = intensityConfigFilename;
		String configFilename;

		if (cmd.getOptionValue("server").equals("intensity")) {
			configFilename = intensityConfigFilename;
		}
		else if (cmd.getOptionValue("server").equals("surface")) {
			configFilename = surfaceConfigFilename;
		}
		else {
			configFilename = defaultConfigFilename;
		}

		SessionFactory sessFactory = new Configuration().configure(configFilename).buildSessionFactory();
		Session sess = sessFactory.openSession();

		Integer siteID = new Integer(0);

		if(sess.isOpen() && sess.isConnected()) {
			if (cmd.hasOption("id")) {
				siteID = Integer.parseInt(cmd.getOptionValue("id"));

				List regions = sess.createCriteria(CyberShakeSiteRegions.class).
				add(Restrictions.eq("id.csSiteId", siteID)).
				list();

				for (int i=0; i<regions.size(); i++) {
					CyberShakeSiteRegions region = (CyberShakeSiteRegions)regions.get(i);
					System.out.println("CS_Site Name: " + region.getCyberShakeSites().getCsSiteName() + 
										", CS_Site_ID: " + region.getId().getCsSiteId() + 
										", ERF_ID: " + region.getId().getErfId() +
										", Cutoff_Dist: " + region.getId().getCutoffDist() + 
										", Max_Lat: " + region.getMaxLat() +
										", Max_Lat_Source_ID: " + region.getMaxLatSourceId() +
										", Max_Lat_Rupture_ID: " + region.getMaxLatRuptureId() +
										", Max_Lon: " + region.getMaxLon() +
										", Max_Lon_Source_ID: " + region.getMaxLonSourceId() +
										", Max_Lon_Rupture_ID" + region.getMaxLonRuptureId() +
										", Min_Lat" + region.getMinLat() +
										", Min_Lat_Source_ID: " + region.getMinLatSourceId() + 
										", Min_Lat_Rupture_ID: " + region.getMinLatRuptureId() + 
										", Min_Lon: " + region.getMinLon() +
										", Min_Lon_Source_ID: " + region.getMinLonSourceId() +
										", Min_Lon_Rupture_ID: " + region.getMinLonRuptureId() );
										
				}
			}
			else if (cmd.hasOption("site")) {
				List regions = sess.createCriteria(CyberShakeSiteRegions.class).
								add(Restrictions.eq("cyberShakeSites.csShortName", "USC")).
								list();
				for (int i=0; i<regions.size(); i++) {
					CyberShakeSiteRegions region = (CyberShakeSiteRegions)regions.get(i);
					System.out.println("CS_Site Name: " + region.getCyberShakeSites().getCsSiteName() + 
										", CS_Site_ID: " + region.getId().getCsSiteId() + 
										", ERF_ID: " + region.getId().getErfId() +
										", Cutoff_Dist: " + region.getId().getCutoffDist() + 
										", Max_Lat: " + region.getMaxLat() +
										", Max_Lat_Source_ID: " + region.getMaxLatSourceId() +
										", Max_Lat_Rupture_ID: " + region.getMaxLatRuptureId() +
										", Max_Lon: " + region.getMaxLon() +
										", Max_Lon_Source_ID: " + region.getMaxLonSourceId() +
										", Max_Lon_Rupture_ID" + region.getMaxLonRuptureId() +
										", Min_Lat" + region.getMinLat() +
										", Min_Lat_Source_ID: " + region.getMinLatSourceId() + 
										", Min_Lat_Rupture_ID: " + region.getMinLatRuptureId() + 
										", Min_Lon: " + region.getMinLon() +
										", Min_Lon_Source_ID: " + region.getMinLonSourceId() +
										", Min_Lon_Rupture_ID: " + region.getMinLonRuptureId() );
										
				}
				
			}
		}

		sess.close();

		if(sess.isOpen() || sess.isConnected()) {
			System.out.println("Error: Hibernate session failed to close");
		}
	}

	public static String getNO_SERVER_OPTION_MESSAGE() {
		return NO_SERVER_OPTION_MESSAGE;
	}

	public static String getNO_SITES_SPECIFIED_MESSAGE() {
		return NO_SITES_SPECIFIED_MESSAGE;
	}
}
