#!/usr/bin/env python

"""
Module CSEPStatus
"""

import sys, os, inspect, subprocess


#--------------------------------------------------------------------------------
#
# CSEPStatus.
#
# This class is designed to acquire current CSEP system and software status.
#
class CSEPStatus:

	# Static data members
	
	# Name of the file with system status information.
	SystemType = "SystemStatus"

	# Name of the file with CSEP software status information.
	SoftwareType = "SoftwareStatus"

	# Dictionary of commands that used to capture external software version and
	# flag if command output is on stderr (True). If command output is redirected
	# to the stderr, the flag should be set to 'True' so CSEP would not trigger
	# it as a failure
	__allPackages = {"gcc --version" : False, # gcc version
					 "awk --version" : False, # awk version
					 "sed --version" : False, # sed version
					 "java -version" : True,  # Java output info to stderr
					 "echo -n $GLOBUS_LOCATION" : False, # globus location
					 "globus-url-copy -version" : True, # globus-url-copy 
					 "globus-job-run -version" : False, # globus-job-run
					 "globus-version" : False, # globus-version
					 "condor_q -version" : False # condor_q
					 }
	
	__options = None
	
	__classMethodForExternalSoftware = 'externalSoftwareVersions'
	
	
	#--------------------------------------------------------------------
	#
	# Initialization.
	#
	# Input:
	#		options - Command-line options including defaults ones used
	#				  by caller program. Default is None.
	# 
	def __init__ (self, sysfile="sys.info", softfile="soft.info"):
		""" Initialization for CSEPStatus class"""
		
		self.__sysfile = sysfile
		self.__softfile = softfile
		

	#--------------------------------------------------------------------
	#
	# Get names of the files for system status.
	#
	# Input: None.
	#
	# Output: A tuple of data and corresponding metadata filenames. 
	#			
	def systemFilename (self):
		""" Get the name of the system status file."""

		return self.__sysfile
	 
	 
	#--------------------------------------------------------------------
	#
	# Get names of the files for system status.
	#
	# Input: None.
	#
	# Output: A list of data filename and corresponding metadata filename. 
	#			
	def softwareFilename (self):
		""" Get the name of the software status file."""
		
		return self.__softfile   
	
	def writeProp(self, fp, prop, value):
		fp.write(str(prop) + ":\n" + str(value).replace("\n", "\\\n") + "\n\n")
	
	def getCommandOutput(self, command, stderr = False):
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		output = proc.communicate()
		out = output[0]
		err = output[1]
		
		if stderr:
			return err
		else:
			return out

	#--------------------------------------------------------------------
	#
	# Capture status of the system.
	#
	# Input:
	#		  filenames - Names of the file and metadata file to capture status 
	#						  information to. Default is None.
	#
	# Output: None.
	#
	def system (self, filenames = None):
		""" Capture system status to the file."""

		if filenames == None:
		   filenames = self.systemFilename()

		# Unpack the sequence
		datafile = filenames

		# Create data file
		fhandle = open(datafile, "w")
		
		# Store host info
		self.writeProp(fhandle, "os.uname", os.uname())
		
		# Store user info
		command = "id"
		self.writeProp(fhandle, command, self.getCommandOutput(command))
		
		# Store environment variables
		self.writeProp(fhandle, "os.environ", os.environ)

		# Store executable and command-line options
		self.writeProp(fhandle, "sys.argv", sys.argv)
		
		# Store executable and command-line options including the default ones
		if self.__options is not None:
		   self.writeProp(fhandle, 
								  "command-line options (including defaults)", 
								  self.__options)
		
		# Close the file
		fhandle.close()
		

	#--------------------------------------------------------------------
	#
	# Capture status of the software used by the system.
	#
	# Input:
	#		  program_name - Name of the calling program.
	#		  program_version - Version of the calling program.
	#		  filenames - Names of the file and metadata file to capture status 
	#						  information to. Default is None.
	#
	# Output: None.
	#
	def software (self, program_name, program_version, filenames = None):
		""" Capture software status to the file."""
			
		if filenames == None:
		   filenames = self.softwareFilename()

		# Unpack the sequence
		datafile = filenames

		# Create data file
		fhandle = open(datafile, "w")

		# Store version of calling program
		self.writeProp(fhandle, program_name, program_version)
		
		# Store python version
		self.writeProp(fhandle, "python sys.version", sys.version)

		for command, output_on_stderr in CSEPStatus.__allPackages.iteritems():
			self.writeProp(fhandle, 
								   command,
								   self.getCommandOutput(command,
															output_on_stderr))

		# Close the file
		fhandle.close()
		
		
	#--------------------------------------------------------------------
	#
	# Get user name of the running process.
	#
	# Input: None.
	#
	# Output: A username.
	#			
	def userName ():
		""" Get the user name of the running process."""

		name = self.getCommandOutput("whoami")
		
		# Strip newline if any
		return name.replace("\n", "")
	 
	userName = staticmethod(userName) 
	

# Invoke the module
if __name__ == '__main__':
	
	numargs = len(sys.argv)
	
	if numargs == 1:
		status = CSEPStatus()
	elif numargs == 3:
		status = CSEPStatus(sys.argv[1], sys.argv[2])
	else:
		sys.stderr.write("USAGE: " + sys.argv[0] + " [sysfile softfile]\n")
		sys.exit(1)
	
	# System status
	filenames = status.systemFilename()
	status.system(filenames)
	
	# Software status
	filenames = status.softwareFilename()
	status.software(CSEPStatus.SoftwareType, "0.0.1", filenames)	 

# end of main	 
						
