#--------------------------------------------------------------------
#
# Compare two ASCII format files.
#
# Input: 
#        filename1 - File name of the first file to compare.
#        filename2 - File name of the second file to compare.
#        precision - Given precision for value comparison. Default is 1e-13.
# 
# Output: 
#        True - files compared OK, exception is raised otherwise.
#
def compare (filename1, filename2, precision = 1E-13):
     """ Compare two files with given precision."""
     
     # Read line at a time from each file, and compare the values
     fhandle1 = open(filename1)
     fhandle2 = open(filename2)

     try:
         while fhandle1 or fhandle2:
             line1 = fhandle1.readline().strip()
             line2 = fhandle2.readline().strip()

             # Ignore lines that begin with strings
             if len(line1) != 0 and len(line2) != 0 and \
                line1[0].isalpha() == True and line2[0].isalpha() == True:

                #CSEPLogging.getLogger(ModuleName).debug("Skipping string lines: %s and %s" \
                #                                        %(line1, line2))
                continue

             line1_tokens = line1.split()
             line2_tokens = line2.split()
             
             if len(line1_tokens) != len(line2_tokens):
                error_msg = "Inconsistent number of elements in lines '%s' vs. '%s'\n" \
                            %(line1, line2)
#                CSEPLogging.getLogger(ModuleName).error(error_msg)
                return False
             
             if len(line1_tokens) == 0:
                break
             
             line1_values = [ float( token ) for token in line1_tokens ]
             line2_values = [ float( token ) for token in line2_tokens ]          

             #CSEPLogging.getLogger(ModuleName).debug("File #1: %s" %line1_values)
             #CSEPLogging.getLogger(ModuleName).debug("File #2: %s" %line2_values)                
             
             for value1, value2 in zip(line1_values, line2_values):
                 diff = value1 - value2
                 
                 if abs(diff) > precision:
                     # Value difference exceeds accepted tolerance, report the error
                     error_msg = "Difference %s (value %s vs. value %s) exceeds \
allowed tolerance %s. Line (%s) vs. line (%s)" \
                                 %(diff, value1, value2, precision, line1, line2)
#                     CSEPLogging.getLogger(ModuleName).error(error_msg)
                     
                     return False
             
     except StandardError, e:
           error_msg = "Error comparing files '%s' and '%s': (%s)" \
                       %(filename1, filename2, e)
#           CSEPLogging.getLogger(ModuleName).error(error_msg)            
           
           return False
     
     # Return True if files compared OK
     return True
 
