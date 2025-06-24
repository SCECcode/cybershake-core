#!/usr/bin/env python3

'''This script is used to investigate the impact of modifying rupture variation probabilities on CyberShake hazard curves.

It takes in a probability modification configuration file, and outputs a list of modified probabilities for rupture variations for use with OpenSHA.

The format of the probability modification configuration file is a series of:

#Comments
Source_ID: Rupture_IDs
    Rupture IDs can be specified as a single value, with commas, and with dashes for ranges
Direction: Percentage
    Direction indicates which part of the fault we want to increase probabilities for; the options are 'N', 'C', or 'S' for north, center, or south.  The percentage is what percent of the hypocenters in this region to select.  For example, S: 25 means to modify probabilities for the southernmost 25% of hypocenters.
Probability
    Increase probability by this much over the default, which is uniform.  Can either be given as a percentage ("50%"), or as "all", in which case all the probability will be assined to these hypocenters.
Direction to decrease probability: Percentage
    The total probability for the rupture must remain constant, so if the probability has been increased somewhere it must be decreased somewhere else.  Possible values are 'N', 'C', and 'S', and the percentage is which hypocenters to decrease probability over.

For example:
90:0-6
S:20
350%
N:80

means that for source 90, ruptures 0-6, select the southernmost 20% of hypocenters and increase their probability by 350% (4.5x).  Then, reduce the probability of the northernmost 80% of hypocenters to preserve the overall rupture probability.
'''

import sys
import os
import pymysql
import re

#If values are within this distance of 0, set them to 0
TOLERANCE = 1.0e-20

class Rup_Prob_Modify:

    def __init__(self, src_id, rup_id, direction=None, sel_per="0", modify="0", anti_direction=None, anti_sel_per="0"):
        self.src_id = src_id
        self.rup_id = rup_id
        #Part of the fault hypocenters are located in to modify probs
        self.direction=direction
        #What % of hypocenters in this part of the fault to modify probs
        self.select_percentage=sel_per
        #How much to modify the prob by
        if modify.find("%")>-1:
            self.modify = float(modify.split("%")[0])
        elif modify=="all":
            self.modify = "all"
        #The part of the fault to modify so that the overall prob is right for this ruptures
        self.anti_direction=anti_direction
        if self.modify=="all" and self.anti_direction is None:
            if self.direction=="N":
                self.anti_direction = "S"
            elif self.direction=="S":
                self.anti_direction = "N"
        if self.direction==self.anti_direction:
            print("Error: The direction %s was selected to both modify and compensate, which is invalid.  Aborting." % self.direction)
            self.exit(1)
        #What % of hypocenters in this part of the fault to modify probs for
        self.anti_select_percentage=anti_sel_per

    def calculate_new_probs(self, erf_id, rup_var_scenario_id, cur):
        #Figure out the prob for this rupture
        print("Calculating new probabilities for src %d, rup %d." % (self.src_id, self.rup_id))
        query = 'select Prob from Ruptures where ERF_ID=%d and Source_ID=%d and Rupture_ID=%d' % (erf_id, self.src_id, self.rup_id)
        cur.execute(query)
        self.prob = float(cur.fetchone()[0])

        rv_dict = dict()
        query = 'select Rup_Var_ID from Rupture_Variations where ERF_ID=%d and Rup_Var_Scenario_ID=%d and Source_ID=%d and Rupture_ID=%d order by Hypocenter_Lat desc' % (erf_id, rup_var_scenario_id, self.src_id, self.rup_id)
        cur.execute(query)
        res = cur.fetchall()
        num_rvs = len(res)
        rv_default_prob = self.prob/num_rvs
        print("Default prob = %e." % rv_default_prob)
        #What number to select and anti-select
        num_rvs_sel = round(num_rvs*float(self.select_percentage)/100.0)
        if self.modify=="all":
            num_rvs_anti_sel = num_rvs - num_rvs_sel
        else:
            num_rvs_anti_sel = round(num_rvs*float(self.anti_select_percentage)/100.0)
        print("RVs preferred: %d" % num_rvs_sel)
        print("RVs antipreferred: %d" % num_rvs_anti_sel)
        tot_prob = 0.0
        if self.direction=="N":
            #We pick the first num_rvs_sel
            for i in range(0, num_rvs_sel):
                rv_id = int(res[i][0])
                if self.modify=="all":
                    rv_prob = self.prob/num_rvs_sel
                else:
                    rv_prob = rv_default_prob*(1.0 + self.modify/100.0)
                rv_dict[rv_id] = rv_prob
                tot_prob += rv_prob
            #Add the centers unchanged
            for i in range(num_rvs_sel, num_rvs-num_rvs_anti_sel):
                rv_id = int(res[i][0])
                rv_dict[rv_id] = rv_default_prob
                tot_prob += rv_default_prob
        elif self.direction=="S":
            #We pick the last num_rvs_sel
            for i in range(num_rvs-num_rvs_sel, num_rvs):
                rv_id = int(res[i][0])
                if self.modify=="all":
                    rv_prob = self.prob/num_rvs_sel
                else:
                    rv_prob = rv_default_prob*(1.0 + self.modify/100.0)
                rv_dict[rv_id] = rv_prob
                tot_prob += rv_prob
            #Add the centers unchanged
            for i in range(num_rvs_anti_sel, num_rvs-num_rvs_sel):
                rv_id = int(res[i][0])
                rv_dict[rv_id] = rv_default_prob
                tot_prob += rv_default_prob
        else:
            print("Don't recognize direction '%s', aborting." % self.direction)
            sys.exit(2)
        #Next, change other probs to compensate and preserve overall rupture prob
        print(self.prob, tot_prob)
        rv_anti_prob = (self.prob-tot_prob)/num_rvs_anti_sel
        if abs(rv_anti_prob)<TOLERANCE:
            rv_anti_prob = 0.0
        if rv_anti_prob<0:
            print("Error: some probabilities for src %d, rup %d are %e, which is less than zero.  Aborting." % (self.src_id, self.rup_id, rv_anti_prob))
            sys.exit(1)
        if self.anti_direction=="N":
            for i in range(0, num_rvs_anti_sel):
                rv_id = int(res[i][0])
                rv_dict[rv_id] = rv_anti_prob
        elif self.anti_direction=="S":
            for i in range(num_rvs-num_rvs_anti_sel, num_rvs):
                rv_id = int(res[i][0])
                rv_dict[rv_id] = rv_anti_prob
        elif self.anti_direction is None:
            for i in range(num_rvs-num_rvs_anti_sel, num_rvs):
                rv_id = int(res[i][0])
                rv_dict
        else:
            print("Don't recognize direction '%s', aborting." % self.anti_direction)
            sys.exit(2)
            
        self.rup_var_prob = rv_dict



def parse_rupture_range(rup_str):
    pattern = r'(\d+)-(\d+)|(\d+)'
    ruptures = []
    matches = re.findall(pattern, rup_str)
    for hit in matches:
        #dash-indicated range
        if hit[0] and hit[1]:
            ruptures.extend([r for r in range(int(hit[0]), int(hit[1])+1)])
        #single rupture
        elif hit[2]:
            ruptures.append(int(hit[2]))
    return ruptures


def parse_prob_mod_file(prob_mod_filename):
    rup_prob_modify_list = []
    with open(prob_mod_filename, 'r') as fp_in:
        data = fp_in.readlines()
        index = 0
        while index<len(data):
            #source and ruptures
            line = data[index]
            if line[0]=="#":
                index += 1
                continue
            (src, rups) = line.strip().split(":")
            rup_list = parse_rupture_range(rups)
            #what to modify
            line = data[index+1]
            (direction, select_percent) = line.strip().split(":")
            #How much to modify 
            line = data[index+2]
            mod = line.strip()
            if mod=="all":
                #Don't need to specify percentage or anti direction
                for r in rup_list:
                    rpm = Rup_Prob_Modify(int(src), r, direction=direction, sel_per=select_percent, modify="all")
                    rup_prob_modify_list.append(rpm)
                index += 3
                continue
            #How to cover the slack
            line = data[index+3]
            (anti_direction, anti_select_percent) = line.strip().split(":")
            for r in rup_list:
                rpm = Rup_Prob_Modify(int(src), r, direction=direction, sel_per=select_percent, modify=mod, anti_direction=anti_direction, anti_sel_per=anti_select_percent)
                rup_prob_modify_list.append(rpm)            
        
            index += 4
        fp_in.close()
    return rup_prob_modify_list

if len(sys.argv)<6:
    print("Usage: %s <site> <erf id> <rup_var_scenario_id> <prob mod file> <output file>" % (sys.argv[0]))
    sys.exit(1)

site = sys.argv[1]
erf_id = int(sys.argv[2])
rup_var_scenario_id = int(sys.argv[3])
prob_mod_filename = sys.argv[4]
output_file = sys.argv[5]

rpm_list = parse_prob_mod_file(prob_mod_filename)

conn = pymysql.connect(host='moment', user='cybershk_ro', passwd='CyberShake2007', db='CyberShake')
cur = conn.cursor()

for rpm in rpm_list:
    rpm.calculate_new_probs(erf_id, rup_var_scenario_id, cur)

with open(output_file, 'w') as fp_out:
    fp_out.write("Source_ID,Rupture_ID,Rup_Var_ID,Probability\n")
    for rpm in rpm_list:
        keys = sorted(rpm.rup_var_prob.keys())
        for k in keys:
            fp_out.write("%d,%d,%d,%e\n" % (rpm.src_id, rpm.rup_id, k, rpm.rup_var_prob[k]))
    fp_out.flush()
    fp_out.close()


'''
query = "select R.Source_ID, R.Rupture_ID, R.Prob from CyberShake_Sites S, CyberShake_Site_Ruptures SR, Ruptures R where S.CS_Short_Name='%s' and S.CS_Site_ID=SR.CS_Site_ID and SR.ERF_ID=%d and SR.ERF_ID=R.ERF_ID and SR.Source_ID=R.Source_ID and SR.Rupture_ID=R.Rupture_ID order by R.Source_ID asc, R.Rupture_ID asc" % (site, erf_id)

cur.execute(query)
res = cur.fetchall()
with open(output_file, "w") as fp_out:
    fp_out.write("Source_ID,Rupture_ID,Rup_Var_ID,Probability\n")
    for r in res:
        (src_id, rup_id, prob) = r
        query = "select count(*) from Rupture_Variations where Source_ID=%s and Rupture_ID=%s and ERF_ID=%d and Rup_Var_Scenario_ID=10" % (src_id, rup_id, erf_id)
        cur.execute(query)
        rv = cur.fetchone()
        num_rvs = int(rv[0])
        for i in range(0, num_rvs):
            fp_out.write("%d,%d,%d,%e\n" % (src_id, rup_id, i, prob/num_rvs))
    fp_out.flush()
    fp_out.close()
'''    
