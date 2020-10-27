package util;

import java.io.File;
import java.io.FilenameFilter;

public class DurationFilenameFilter implements FilenameFilter {

	public boolean accept(File dir, String name) {
		if (name.endsWith(".dur"))
			return true;
		else
			return false;
	}

}
