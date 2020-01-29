import java.io.BufferedReader;

import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import acm.util.ErrorException;

/*
 * File: NameSurferDataBase.java
 * -----------------------------
 * This class keeps track of the complete database of names.
 * The constructor reads in the database from a file, and
 * the only public method makes it possible to look up a
 * name and get back the corresponding NameSurferEntry.
 * Names are matched independent of case, so that "Eric"
 * and "ERIC" are the same names.
 */

public class NameSurferDataBase implements NameSurferConstants {

	/* Constructor: NameSurferDataBase(filename) */
	/**
	 * Creates a new NameSurferDataBase and initializes it using the data in the
	 * specified file. The constructor throws an error exception if the requested
	 * file does not exist or if an error occurs as the file is being read.
	 */
	public NameSurferDataBase(String filename) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String line = "";
			while (true) {
				line = br.readLine();
				if (line == null)
					break;
				entries.put(getName(line), new NameSurferEntry(line));
			}
			br.close();
		} catch (IOException e) {
			throw new ErrorException(e);
		}
	}

	// gets name separated from the given line
	private String getName(String line) {
		String name = line.substring(0, line.indexOf(" "));
		return name;
	}

	/* Method: findEntry(name) */
	/**
	 * Returns the NameSurferEntry associated with this name, if one exists. If the
	 * name does not appear in the database, this method returns null.
	 */
	public NameSurferEntry findEntry(String name) {
		name = convertCase(name);
		if (!entries.containsKey(name)) {
			return null;
		}
		return entries.get(name);
	}

	// converts entered name so that hashmap can perceive it
	private String convertCase(String name) {
		if (name.length() > 0)
			name = name.substring(0, 1).toUpperCase() + name.substring(1).toLowerCase();
		return name;
	}

	private HashMap<String, NameSurferEntry> entries = new HashMap<>();
}