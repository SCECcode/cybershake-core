package util;

import java.io.File;
import java.io.FilenameFilter;


public class BSAFilenameFilter implements FilenameFilter {

	public boolean accept(File dir, String name) {
		if (name.endsWith(".bsa"))
			return true;
		else
			return false;
		
	}


	
	

}
