package commands;

import java.util.List;

import mapping.CyberShakeSites;

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

public class CyberSites {
	
	private static final String NO_OPTION_FOR_RETRIEVING_DATA = "Please specify an option that retrieves data or use -help";
	private static final String NO_SERVER_OPTION_MESSAGE = "Please use -server to specify a database server";
	
	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		Options options = new Options();

		Option all = new Option("all","list all information for all sites");
		Option location = OptionBuilder.withArgName("short name").hasArg().withDescription("provides location data for the specified site").create("l");
		Option shortname = OptionBuilder.withArgName("short name").hasArg().withDescription("provides all the details for the specified site").create("n");
		Option id = OptionBuilder.withArgName("ID").hasArg().withDescription("provides all the details for the specified site").create("id");
		Option server = OptionBuilder.withArgName("name").hasArg().withDescription("server name (surface or intensity) - this options is required").create("server");

		Option summary = new Option("summary","provides all the short names for all the sites");
		Option help = new Option("help", "print this message");

		options.addOption(all);
		options.addOption(location);
		options.addOption(shortname);
		options.addOption(id);
		options.addOption(server);

		options.addOption(summary);
		options.addOption(help);

		CommandLineParser parser = new GnuParser();

		try {
			CommandLine cmd = parser.parse( options, args);

			if (cmd.hasOption("help")) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp( "CyberSites", options, true );	
			}
			else {
				if (!cmd.hasOption("server")) {
					System.out.println(NO_SERVER_OPTION_MESSAGE);
					return;
				}

				if (!cmd.hasOption("all") && 
					!cmd.hasOption("l") && 
					!cmd.hasOption("id") && 
					!cmd.hasOption("n") &&
					!cmd.hasOption("summary")) {
					System.out.println(NO_OPTION_FOR_RETRIEVING_DATA);
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

		if(sess.isOpen() && sess.isConnected()) {
			if (cmd.hasOption("l")) {
				getLocation(cmd, sess);
			}
			else if (cmd.hasOption("n")) {
				getSiteUsingShortName(cmd, sess);
			}
			else if (cmd.hasOption("id")) {
				getSiteUsingID(cmd,sess);
			}
			else if (cmd.hasOption("summary")) {
				getSummary(sess);
			}
			else if (cmd.hasOption("all")) {
				getAllSites(sess);
			}
		}

		sess.close();

		if(sess.isOpen() || sess.isConnected()) {
			System.out.println("Error: Hibernate session failed to close");
		}
		
	}

	private static void getAllSites(Session sess) {
		List sites = sess.createQuery("from CyberShakeSites").list();
		for (int i=0; i<sites.size(); i++) {
			CyberShakeSites site = (CyberShakeSites)sites.get(i);
			System.out.println(site);
		}
		
	}

	private static void getSummary(Session sess) {
		List sites = sess.createQuery("from CyberShakeSites").list();
		for (int i=0; i<sites.size(); i++) {
			CyberShakeSites site = (CyberShakeSites)sites.get(i);
			System.out.println("CS_Short_Name: " + site.getCsShortName());
		}
		
	}

	private static void getSiteUsingID(CommandLine cmd, Session sess) {
		Integer siteID;
		siteID = Integer.parseInt(cmd.getOptionValue("id"));

		List sites = sess.createCriteria(CyberShakeSites.class).
			add(Restrictions.eq("csSiteId", siteID)).
			list();

		for (int i=0; i<sites.size(); i++) {
			CyberShakeSites site = (CyberShakeSites)sites.get(i);
			System.out.println(site);
		}
	}

	private static void getSiteUsingShortName(CommandLine cmd, Session sess) {
		List sites = sess.createCriteria(CyberShakeSites.class).
			add(Restrictions.eq("csShortName", cmd.getOptionValue("n"))).
			list();
		for (int i=0; i<sites.size(); i++) {
			CyberShakeSites site = (CyberShakeSites)sites.get(i);
			System.out.println(site);
		}		
	}

	private static void getLocation(CommandLine cmd, Session sess) {
		List sites = sess.createCriteria(CyberShakeSites.class).
		add(Restrictions.eq("csShortName", cmd.getOptionValue("l"))).
		list();
		for (int i=0; i<sites.size(); i++) {
			CyberShakeSites site = (CyberShakeSites)sites.get(i);
			System.out.println("CS_Short_Name: " + site.getCsShortName() + 
					" , CS_Site_Lat: " + site.getCsSiteLat() + 
					" , CS_Site_Lon: " + site.getCsSiteLon());
		}		

	}
}
