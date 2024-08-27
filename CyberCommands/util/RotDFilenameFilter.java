package util;

import java.io.File;
import java.io.FilenameFilter;

import javax.swing.filechooser.FileFilter;

public class RotDFilenameFilter implements FilenameFilter {

	public boolean accept(File dir, String name) {
		if (name.startsWith("RotD") && name.endsWith(".rotd")) {
			return true;
		}
		return false;
	}


}
