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
        if (args.length<5) {
            System.out.println("Usage:  CheckDBDataForSite <site> <erf_id> <rup_var_scen_id> <sgt_variation_id> <outputFile>");
            System.exit(1);
        }
        
        String site = args[0];
        String erf_id = args[1];
        String rup_var_scen_id = args[2];
        String sgt_variation_id = args[3];
        String outputFile = args[4];
        
        checkDBData(site, erf_id, rup_var_scen_id, sgt_variation_id, outputFile);
    }

    private static void checkDBData(String site, String erf_id, String rup_var_scen_id, String sgt_variation_id, String outputFile) {
        DBConnect dbc = new DBConnect(DB_SERVER, DB);
        
        String query = "select count(*) " +
        "from Rupture_Variations V, CyberShake_Site_Ruptures R, CyberShake_Sites S " +
        "where S.CS_Short_Name='" + site + "' " +
        "and S.CS_Site_ID=R.CS_Site_ID " +
        "and R.ERF_ID=" + erf_id + " " +
        "and V.ERF_ID=" + erf_id + " " +
        "and R.Source_ID=V.Source_ID " +
        "and R.Rupture_ID=V.Rupture_ID " +
        "and V.Rup_Var_Scenario_ID=" + rup_var_scen_id + " ";
        
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
            "from PeakAmplitudes A, CyberShake_Sites S " +
            "where S.CS_Short_Name='" + site + "' " +
            "and S.CS_Site_ID=A.Site_ID " +
            "and A.ERF_ID=" + erf_id + " " +
            "and A.Rup_Var_Scenario_ID=" + rup_var_scen_id + " " +
            "and A.SGT_Variation_ID=" + sgt_variation_id;

            System.out.println(query);
            
            ResultSet ampSet = dbc.selectData(query);
            ampSet.first();
            if (ampSet.getRow()==0) {
                System.err.println("No rup vars in DB.");
                System.exit(2);
            }
            int ampSetNum = ampSet.getInt("count(*)");
            
            if (rupVarSetNum!=ampSetNum/4) {
                System.out.println(rupVarSetNum + " variations for site " + site + " in RupVar table, but " + (ampSetNum/4) + " variations in PeakAmp table.");
                rupVarSet.close();
                ampSet.close();
                findDifferences(dbc, site, erf_id, rup_var_scen_id, sgt_variation_id, outputFile);
            } else {
            	System.out.println("All ruptures are registered.");
            }
            
        } catch (SQLException e) {
            e.printStackTrace();
            System.exit(-3);
        }
        
        
    }

    /**
     * @param site
     * @param erf_id
     * @param rup_var_id
     */
    private static void findDifferences(DBConnect dbc, String site, String erf_id, String rup_var_scen_id, String sgt_variation_id, String outputFile) {
        String query = "select V.Source_ID, V.Rupture_ID, V.Rup_Var_ID " +
        "from Rupture_Variations V, CyberShake_Site_Ruptures R, CyberShake_Sites S " +
        "where S.CS_Short_Name='" + site + "' " +
        "and S.CS_Site_ID=R.CS_Site_ID " +
        "and R.ERF_ID=" + erf_id + " " +
        "and V.ERF_ID=" + erf_id + " " +
        "and R.Source_ID=V.Source_ID " +
        "and R.Rupture_ID=V.Rupture_ID " + 
        "and V.Rup_Var_Scenario_ID=" + rup_var_scen_id;
        
        ResultSet rupVarSet = dbc.selectData(query);
        try {
            rupVarSet.first();
            if (rupVarSet.getRow()==0) {
                System.err.println("No entries in DB.");
                System.exit(3);
            }
            String prefix = "select * " +
            "from PeakAmplitudes A, CyberShake_Sites S " +
            "where S.CS_Short_Name='" + site + "' " +
            "and S.CS_Site_ID=A.Site_ID " +
            "and A.ERF_ID=" + erf_id + " " +
            "and A.Rup_Var_Scenario_ID=" + rup_var_scen_id + " " +
            "and A.SGT_Variation_ID=" + sgt_variation_id + " ";

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
