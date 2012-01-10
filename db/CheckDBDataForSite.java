import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;

/*
 * Created on Jun 30, 2007
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

/**
 * @author Scott
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class CheckDBDataForSite {
    private static final String DB_SERVER = "focal.usc.edu";
    private static final String DB = "CyberShake";
    
    public static void main(String[] args) {
        if (args.length<2) {
            System.out.println("Usage:  CheckDBDataForSite <runID> <outputFile>");
            System.exit(1);
        }
        
        String runID = args[0];
        String outputFile = args[1];
        
        checkDBData(runID, outputFile);
    }

    private static void checkDBData(String runID, String outputFile) {
        DBConnect dbc = new DBConnect(DB_SERVER, DB);
        
        
        String query = "select count(*) " +
        "from Rupture_Variations V, CyberShake_Runs U, CyberShake_Site_Ruptures R " +
        "where U.Run_ID=" + runID + " " + 
        "and R.CS_Site_ID=U.Site_ID " +
        "and R.ERF_ID=U.ERF_ID " +
        "and V.ERF_ID=U.ERF_ID " +
        "and R.Source_ID=V.Source_ID " +
        "and R.Rupture_ID=V.Rupture_ID " +
        "and V.Rup_Var_Scenario_ID=U.Rup_Var_Scenario_ID";
        
        System.out.println(query);
        
        ResultSet rupVarSet = dbc.selectData(query);
        
        try {
            rupVarSet.first();
            if (rupVarSet.getRow()==0) {
                System.err.println("No rup vars in DB.");
                System.exit(2);
            }
            
            int rupVarSetNum = rupVarSet.getInt("count(*)");
            
            query = "select count(*) " +
            "from PeakAmplitudes A " +
            "where A.Run_ID=" + runID + " ";

            System.out.println(query);
            
            ResultSet ampSet = dbc.selectData(query);
            ampSet.first();
            if (ampSet.getRow()==0) {
                System.err.println("No rup vars in DB.");
                System.exit(2);
            }
            int ampSetNum = ampSet.getInt("count(*)");
            
            if (rupVarSetNum!=ampSetNum/4) {
                System.out.println(rupVarSetNum + " variations for run " + runID + " in RupVar table, but " + (ampSetNum/4) + " variations in PeakAmp table.");
                rupVarSet.close();
                ampSet.close();
                findDifferences2(dbc, runID, outputFile);
                System.exit(3);
            } else {
            	System.out.println("All ruptures are registered.");
            }
            
        } catch (SQLException e) {
            e.printStackTrace();
            System.exit(-3);
        }
        
        
    }

    private static class SourceRuptureVar implements Comparable {
    	private int Source_ID;
    	private int Rupture_ID;
    	private int Rupture_Variation_ID;
    	
    	public SourceRuptureVar(int sid, int rid, int vid) {
    		Source_ID = sid;
    		Rupture_ID = rid;
    		Rupture_Variation_ID = vid;
    	}
    	
    	public boolean equals(Object o) {
    		if (! (o instanceof SourceRuptureVar)) {
    			return false;
    		}
    		SourceRuptureVar sr = (SourceRuptureVar)o;
    		if (sr.Rupture_Variation_ID==Rupture_Variation_ID && sr.Rupture_ID==Rupture_ID && sr.Source_ID==Source_ID) {
    			return true;
    		}
    		return false;
    	}
    	
    	public int hashCode() {
    		return (Source_ID*10000+Rupture_ID)%Rupture_Variation_ID;
    	}

		public int compareTo(Object o) {
			if (! (o instanceof SourceRuptureVar)) {
				throw new ClassCastException();
			}
			SourceRuptureVar srv = (SourceRuptureVar)o;
//			System.out.println("Comparing " + srv.Source_ID + " " + srv.Rupture_ID + " " + srv.Rupture_Variation_ID +
//					" to " + Source_ID + " " + Rupture_ID + " " + Rupture_Variation_ID);
			if (Source_ID!=srv.Source_ID) {
				return (Source_ID - srv.Source_ID);
			} else {
				if (Rupture_ID!=srv.Rupture_ID) {
					return (Rupture_ID - srv.Rupture_ID);
				} else {
					return (Rupture_Variation_ID - srv.Rupture_Variation_ID);
				}
			}
		}

		public int getSource_ID() {
			return Source_ID;
		}

		public int getRupture_ID() {
			return Rupture_ID;
		}

		public int getRupture_Variation_ID() {
			return Rupture_Variation_ID;
		}
    }
    
    private static void findDifferences2(DBConnect dbc, String runID, String outputFile) {
    	System.out.println(System.currentTimeMillis());
    	
    	String rupVarQuery = "select V.Source_ID, V.Rupture_ID, V.Rup_Var_ID " +
    	"from Rupture_Variations V, CyberShake_Site_Ruptures R, CyberShake_Runs U " +
    	"where U.Run_ID=" + runID + " " +
    	"and U.Site_ID=R.CS_Site_ID " +
    	"and U.ERF_ID=R.ERF_ID " +
    	"and R.ERF_ID=V.ERF_ID " +
    	"and U.Rup_Var_Scenario_ID=V.Rup_Var_Scenario_ID " +
    	"and R.Source_ID=V.Source_ID " +
    	"and R.Rupture_ID=V.Rupture_ID " +
    	"order by V.Source_ID asc, V.Rupture_ID asc, V.Rup_Var_ID asc";
    
    	ResultSet rupVarSet = dbc.selectData(rupVarQuery);
    	TreeSet<SourceRuptureVar> rupVarTreeSet = null;
    	
    	try {
    		rupVarSet.first();
    		if (rupVarSet.getRow()==0) {
    			System.err.println("No rup var entries in DB.");
    			System.exit(3);
    		}
    	   	rupVarTreeSet = new TreeSet<SourceRuptureVar>();
        	while (!rupVarSet.isAfterLast()) {
        		SourceRuptureVar sr = new SourceRuptureVar(rupVarSet.getInt("V.Source_ID"), rupVarSet.getInt("V.Rupture_ID"), rupVarSet.getInt("V.Rup_Var_ID"));
        		rupVarTreeSet.add(sr);
        		rupVarSet.next();
        	}
        	rupVarSet.close();
    	} catch (SQLException ex) {
    		ex.printStackTrace();
    		System.exit(-3);
    	}
    	

    	
    	System.out.println(System.currentTimeMillis());
    	
    	String peakAmpsQuery = "select distinct Source_ID, Rupture_ID, Rup_Var_ID " +
    	"from PeakAmplitudes " +
    	"where Run_ID=" + runID + " " + 
    	"order by Source_ID asc, Rupture_ID asc, Rup_Var_ID asc";
    	
    	ResultSet peakAmpsSet = dbc.selectData(peakAmpsQuery);
    	TreeSet<SourceRuptureVar> peakAmpsTreeSet = null;
    	
    	System.out.println(System.currentTimeMillis());
    	
    	try {
    		peakAmpsSet.first();
    		if (peakAmpsSet.getRow()==0) {
    			System.err.println("No peak amps entries in DB.");
    			System.exit(3);
    		}
    		peakAmpsTreeSet = new TreeSet<SourceRuptureVar>();
        	while (!peakAmpsSet.isAfterLast()) {
        		SourceRuptureVar sr = new SourceRuptureVar(peakAmpsSet.getInt("Source_ID"), peakAmpsSet.getInt("Rupture_ID"), peakAmpsSet.getInt("Rup_Var_ID"));
        		peakAmpsTreeSet.add(sr);
        		peakAmpsSet.next();
        	}
        	peakAmpsSet.close();
    	} catch (SQLException ex) {
    		ex.printStackTrace();
    		System.exit(-3);
    	}
    	
    	BufferedWriter bw;
    	
		int count = 0;
		int missing = 0;
    	
    	try {
    		bw = new BufferedWriter(new FileWriter(outputFile));
    		int lastSource = -1;
    		int lastRupture = -1;
    		
    		Iterator<SourceRuptureVar> rupVarIterator = rupVarTreeSet.iterator();
    		ArrayList<String> listForRupture = new ArrayList<String>();
    		while (rupVarIterator.hasNext()) {
    			count++;
    			if (count%10000==0) {
    				System.out.println("Processed " + count + " ruptures.");
    			}
    			SourceRuptureVar srv = rupVarIterator.next();
    			if (!peakAmpsTreeSet.contains(srv)) {
    				missing++;
    				if (srv.getRupture_ID()!=lastRupture || srv.getSource_ID()!=lastSource) {
    					if (lastSource!=-1) {
    						bw.write(listForRupture.get(0) + " " + (listForRupture.size()-1) + "\n");
    						for (int i=1; i<listForRupture.size(); i++) {
    							bw.write(listForRupture.get(i) + "\n");
    						}
    					}
    					listForRupture = new ArrayList<String>();
    					listForRupture.add(srv.getSource_ID() + " " + srv.getRupture_ID());
    					lastSource = srv.getSource_ID();
    					lastRupture = srv.getRupture_ID();
    				}
					listForRupture.add(srv.getRupture_Variation_ID() + "");
    			}
    		}
    		
    		//write any stragglers
    		bw.write(listForRupture.get(0) + " " + (listForRupture.size()-1) + "\n");
			for (int i=1; i<listForRupture.size(); i++) {
				bw.write(listForRupture.get(i) + "\n");
			}
    		
    		bw.flush();
    		bw.close();
    		
    	} catch (IOException ex) {
    		ex.printStackTrace();
    		System.exit(-4);
    	}
    	System.out.println("Total missing: " + missing);
    }
    /**
     * @param site
     * @param erf_id
     * @param rup_var_id
     */
    private static void findDifferences(DBConnect dbc, String runID, String outputFile) {
        String query = "select V.Source_ID, V.Rupture_ID, V.Rup_Var_ID " +
        "from Rupture_Variations V, CyberShake_Site_Ruptures R, CyberShake_Runs U " +
        "where U.Run_ID=" + runID + " " +
        "and U.Site_ID=R.CS_Site_ID " +
        "and U.ERF_ID=R.ERF_ID " +
        "and R.ERF_ID=V.ERF_ID " +
        "and U.Rup_Var_Scenario_ID=V.Rup_Var_Scenario_ID " +
        "and R.Source_ID=V.Source_ID " +
        "and R.Rupture_ID=V.Rupture_ID "; 
        
        ResultSet rupVarSet = dbc.selectData(query);
        try {
            rupVarSet.first();
            if (rupVarSet.getRow()==0) {
                System.err.println("No entries in DB.");
                System.exit(3);
            }
            String prefix = "select * " +
            "from PeakAmplitudes A " +
            "where A.Run_ID=" + runID + " ";
            
            BufferedWriter bw = null;
            try {
                bw = new BufferedWriter(new FileWriter(outputFile));
                
                int count = 0;
                int lastSource = -1;
                int lastRupt = -1;
                ArrayList<String> listForRupture = new ArrayList<String>();
                while (!rupVarSet.isAfterLast()) {
                    count++;
                    int source = rupVarSet.getInt("Source_ID");
                    int rupture = rupVarSet.getInt("Rupture_ID");
                    
                    if (count%100==0) {
                        System.out.println("Completed " + count + " variations.");
                    }
                    
                    query = prefix + "and A.Source_ID=" + source + " " +
                    "and A.Rupture_ID=" + rupture + " " +
                    "and A.Rup_Var_ID=" + rupVarSet.getInt("Rup_Var_ID");
//                    System.out.println(query);
                    ResultSet psaSet = dbc.selectData(query);
                    psaSet.first();
                    if (psaSet.getRow()==0) { //no entries, missing
                        if (lastSource!=source || lastRupt!=rupture) {
                            if (lastSource!=-1) {
                                bw.write(listForRupture.get(0) + " " + (listForRupture.size()-1) + "\n");
                                for (int i=1; i<listForRupture.size(); i++) {
                                    bw.write(listForRupture.get(i) + "\n");
                                }
                            }
                            listForRupture.clear();
                            listForRupture.add(source + " " + rupture);
                            listForRupture.add("" + rupVarSet.getInt("Rup_Var_ID"));
                            lastSource = source;
                            lastRupt = rupture;
                            System.out.println(query);
                        } else {
                            listForRupture.add("" + rupVarSet.getInt("Rup_Var_ID"));
                        }
                    }
                    
                    psaSet.close();
                    rupVarSet.next();
                }

                bw.write(listForRupture.get(0) + " " + (listForRupture.size()-1) + "\n");
                for (int i=1; i<listForRupture.size(); i++) {
                    bw.write(listForRupture.get(i) + "\n");
                }

                bw.flush();
                bw.close();
                rupVarSet.close();
            } catch (IOException e1) {
                e1.printStackTrace();
                System.exit(-1);
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-2);
        }        
    }
}
