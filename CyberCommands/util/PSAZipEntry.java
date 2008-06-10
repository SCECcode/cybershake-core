package util;

import java.util.zip.ZipEntry;

public class PSAZipEntry extends ZipEntry {

	public PSAZipEntry(ZipEntry e) {
		super(e);
	}

	public PSAZipEntry(String name) {
		super(name);
	}

}
