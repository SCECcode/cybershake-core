package util;

import java.io.File;
import java.io.FileFilter;

public class DirFileFilter implements FileFilter {

	public boolean accept(File pathname) {
		if (pathname.isDirectory()) 
			return true;
		else 
			return false;
	}

}
