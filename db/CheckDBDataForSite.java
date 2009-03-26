import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;

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
    	
    	try {
    		rupVarSet.first();
    		if (rupVarSet.getRow()==0) {
    			System.err.println("No rup var entries in DB.");
    			System.exit(3);
    		}
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
    	
    	try {
    		peakAmpsSet.first();
    		if (peakAmpsSet.getRow()==0) {
    			System.err.println("No peak amps entries in DB.");
    			System.exit(3);
    		}
    	} catch (SQLException ex) {
    		ex.printStackTrace();
    		System.exit(-3);
    	}
    	
    	System.out.println(System.currentTimeMillis());
    	
    	BufferedWriter bw;
    	int rupValSource, rupValRupture, rupValRupVar;
    	int peakAmpsSource, peakAmpsRupture, peakAmpsRupVar;
    	
		int count = 0;
		int missing = 0;
    	
    	try {
    		bw = new BufferedWriter(new FileWriter(outputFile));
    		
    		while (!rupVarSet.isAfterLast()) {
        		if (count % 10000 == 0){
        			System.out.println("Processed " + count + " rupture variations.");
        		}
        		
    			rupValSource = rupVarSet.getInt("V.Source_ID");
    			rupValRupture = rupVarSet.getInt("V.Rupture_ID");
    			rupValRupVar = rupVarSet.getInt("V.Rup_Var_ID");
    			
    			peakAmpsSource = peakAmpsSet.getInt("Source_ID");
    			peakAmpsRupture = peakAmpsSet.getInt("Rupture_ID");
    			peakAmpsRupVar = peakAmpsSet.getInt("Rup_Var_ID");
    			
    			if (rupValSource != peakAmpsSource || rupValRupture != peakAmpsRupture || rupValRupVar != peakAmpsRupVar) {
    				bw.write("Mismatch between " + rupValSource + ", " + rupValRupture + ", " + rupValRupVar +
    						" and " + peakAmpsSource + ", " + peakAmpsRupture + ", " + peakAmpsRupVar + "\n");
    				bw.flush();
    				rupVarSet.next();
    				missing++;
//    				while (rupValSource != peakAmpsSource || rupValRupture != peakAmpsRupture || rupValRupVar != peakAmpsRupVar) {
//    					rupVarSet.next();
//    	    			rupValSource = rupVarSet.getInt("V.Source_ID");
//    	    			rupValRupture = rupVarSet.getInt(RUPTURE_COL_ID);
//    	    			rupValRupVar = rupVarSet.getInt(RUPT_VAR_COL_ID);
//    				}
    			} else {
    				peakAmpsSet.next();
    				rupVarSet.next();
    			}
    			
//    			rupVarSet.next();
    			count++;
    		}
    		System.out.println(count);
    		bw.flush();
    		bw.close();
    		
    	} catch (IOException ex) {
    		ex.printStackTrace();
    		System.exit(-4);
    	} catch (SQLException ex) {
    		System.out.println("line " + count);
    		ex.printStackTrace();
    		System.exit(-5);
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
