package util;

import java.io.File;
import java.io.FilenameFilter;

public class ZipFilenameFilter implements FilenameFilter {

	public boolean accept(File dir, String name) {
		if (name.endsWith("PSA.zip"))
			return true;
		else
			return false;
		
	}

}
