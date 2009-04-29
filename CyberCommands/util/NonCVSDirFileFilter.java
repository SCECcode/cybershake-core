package util;

import java.io.File;
import java.io.FileFilter;

public class NonCVSDirFileFilter implements FileFilter {

	public boolean accept(File pathname) {
		if (pathname.isDirectory() && !pathname.equals("CVS")) 
			return true;
		else 
			return false;
	}

}
