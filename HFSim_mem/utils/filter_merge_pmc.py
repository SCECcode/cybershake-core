#!/usr/bin/env python3

import sys
import os

if len(sys.argv)<4:
    print("Usage: %s <merge pmc infile> <CSV with source,rupture pairs to keep> <merge pmc outfile>" % sys.argv[0])
    sys.exit(1)

infile = sys.argv[1]
rup_infile = sys.argv[2]
outfile = sys.argv[3]

rups_to_include = []

with open(rup_infile, "r") as fp_in:
    data = fp_in.readlines()
    for line in data:
        rups_to_include.append(line.split(',')[0:2])
    fp_in.close()

output_data = []

with open(infile, 'r') as fp_in:
    data = fp_in.readlines()
    fp_in.close()
    #merge in file looks like
    #@ 1 scec::HF_Synth:3.0 HF_Synth_128_958_t9 
    #TASK HF_Synth_128_958_t9 ...
    #<snip>
    #EDGE EDGE HF_Synth_128_958_t9 Combine_HF_Synth_128_958
    i = 0
    job_index = 1
    while(i<len(data)):
        if (i%1000==0):
            print("Line %d of %d." % (i, len(data)))
        line = data[i]
        if line[0:2]=="#@":
            pieces = line.strip().split()
            #jobtype = pieces[2]
            if line.find(":HF_Synth:")>-1:
                (src_id, rup_id) = pieces[3].split("_")[2:4]
            elif line.find("Combine_HF_Synth_")>-1:
                (src_id, rup_id) = pieces[3].split("_")[3:5]
            elif line.find("Combine_PGA_")>-1:
                (src_id, rup_id) = pieces[3].split("_")[2:4]
            elif line.find(":LF_Site_Response:")>-1:
                (src_id, rup_id) = pieces[3].split("_")[3:5]
            elif line.find(":MergeIM:")>-1:
                (src_id, rup_id) = pieces[3].split("_")[2:4]
            else:
                print("Jobtype %s not recognized, aborting." % line)
                sys.exit(2)
            for (s, r) in rups_to_include:
                if src_id==s and rup_id==r:
                    output_data.append("#@ %d %s %s\n" % (job_index, pieces[2], pieces[3]))
                    output_data.append(data[i+1])
                    job_index += 1
                    break
            i += 2
        elif line[0:4]=="EDGE":
            #Last two in line are (src id, rup_id)
            (src_id, rup_id) = line.strip().split("_")[-2:]
            for (s, r) in rups_to_include:
                if src_id==s and rup_id==r:
                    output_data.append(line)
                    break
            i += 1
        elif line=="\n":
            #Empty line between tasks and edges
            i += 1
        else:
            print("Not sure what to do with line #%d ('%s'), aborting." % (i, line))
            sys.exit(3)
        
with open(outfile, 'w') as fp_out:
    for line in output_data:
        fp_out.write(line)
    fp_out.flush()
    fp_out.close()


