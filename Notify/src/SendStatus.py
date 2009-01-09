#!/usr/bin/python

import sys
import socket
import os
import pwd
import smtplib


# Global vars
notify_from = ""
notify_to = []
site = ""
workflow = ""
stage_info = []


def doSGT(args):
    stage = args[0]
    stages = ["preCVM", "vMeshGen", "vMeshMerge", "sgtGenXY", "sgtMergeXY"]
    subject = "Status - " + site + " " + workflow + " Workflow" 
    msg = "Workflow is currently on this stage:\r\n\r\n"

    try:
        numstage = len(stages)
        istage = stages.index(stage)
        for i in range(0, numstage):
            if (i <= istage):
                msg = msg + stages[i] + ": Complete\r\n"
            else:
                msg = msg + stages[i] + ": Scheduled\r\n"
    except:
        print "Stage " + stage + " not found"
        return 1
    
    # Send the email
    sendNotification(subject, msg, notify_to)

    return 0


def doPP(args):
    stage = args[0]
    daxnum = int(args[1])
    maxdax = int(args[2])
    stages = ["CheckSgt", "DAX"]
    subject = "Status - " + site + " " + workflow + " Workflow"
    msg = "Workflow is currently on this stage:\r\n\r\n"   
    
    try:
        istage = stages.index(stage)
    except:
        print "Stage " + stage + " not found"
        return 1
    
    if (daxnum > maxdax):
        print "DAX number " + str(daxnum) + " is greater than max DAX " + str(maxdax)
        return 1
    
    msg = msg + stages[0] + ": Complete\r\n"
    if (istage > 0):
        msg = msg + stages[1] + ": Number " + str(daxnum) + " of approx " + str(maxdax) + " completed successfully\r\n"
    else:
        msg = msg + stages[1] + ": Scheduled\r\n"
                        
    # Send the email
    sendNotification(subject, msg, notify_to)
    
    return 0


# Workflow definitions
# Mapping of workflow name -> tuple (number of arguments, handler)
WORKFLOWS = {"SGT":(1, doSGT), \
             "PP":(3, doPP)}


# Send email msg using SMTP
def sendNotification(subject, msg, notify_user):
    to_str = ""
    for n in notify_user:
        if (to_str == ""):
            to_str = n
        else:
            to_str = to_str + "," + n 
    msg = "From: " + notify_from + \
        "\r\nTo: " + to_str + \
        "\r\nSubject: " + subject + \
        "\r\n" + msg + \
        "\r\n---------------------------------------------------\r\nAutomated msg from Workflow Status\r\n"
    
    server = smtplib.SMTP('localhost')
    #server.set_debuglevel(1)
    server.sendmail(notify_from, notify_user, msg)
    server.quit()

    return 0


def init():
    global notify_from
    global notify_to
    global site
    global workflow
    global stage_info
    
    # Get the current user id
    userid = pwd.getpwuid(os.getuid())[0]
    domain = socket.getfqdn()
    notify_from = userid + "@" + domain
    
    # Get number of command-line arguments
    argc = len(sys.argv)
    
    # Parse command line arguments
    if (argc < 5):
        print "Usage: " + sys.argv[0] + " <site> <workflow> <notify_list file> <stage info>"
        print "Example: " + sys.argv[0] + " USC SGT notify.file PreCVM"
        print "Example: " + sys.argv[0] + " USC PP notify.file CheckSgt 12 80"
        return 1
            
    site = sys.argv[1]
    workflow = sys.argv[2]
    notify_file = sys.argv[3]
    stage_info = sys.argv[4:]

    print "Configuration:"
    print "Notify From:\t" + notify_from
    print "Notify File:\t" + notify_file
    print "Site:\t\t" + site
    print "Workflow:\t" + workflow
    print "Stage Info:\t" + str(stage_info)

    # Check that the workflow is valid and the correct number of arguments were supplied
    try:
        wfinfo = WORKFLOWS[workflow]
        numargs = wfinfo[0]
        if (len(stage_info) != numargs):
            print "Workflow " + workflow + " requires " + str(numargs) + " argument(s)"
            return 1
    except:
        print "Unable to find " + workflow
        return 1        
    
    # Load the notify list from the file
    try:
        file = open(notify_file)
        lines = file.read()
        lines = lines.splitlines()
        for l in lines:
            if (len(l) > 0):
                notify_to.append(l)
    except:
        print "Unable to read file " + notify_file
        return 1

    print "Notification List:"
    for n in notify_to:
        print " " + n
    
    return 0


def main():
    if (len(notify_to) == 0):
        print "No users specified in notify file - notifications disabled"
    else:
        # Execute the workflow handler
        retcode = WORKFLOWS[workflow][1](stage_info)
        if (retcode != 0):
            print "Error sending email notification"
            return 1

    return 0


def cleanup():
    return 0


if __name__ == '__main__':
    if (init() != 0):
        sys.exit(1)
    if (main() != 0):
        sys.exit(1)
    cleanup()
    sys.exit(0)
