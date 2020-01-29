
/*
 * File: HangmanLexicon.java
 * -------------------------
 * This file contains a stub implementation of the HangmanLexicon
 * class that you will reimplement for Part III of the assignment.
 */

import acm.util.*;
import java.io.*;
import java.util.*;

public class HangmanLexicon {
	private ArrayList<String> arr = new ArrayList<String>();

	// This is the HangmanLexicon constructor
	public HangmanLexicon() {
		try {
			BufferedReader br = new BufferedReader(new FileReader("HangmanLexicon.txt"));
			while (true) {
				String str = br.readLine();
				if (str == null)
					break;
				arr.add(str);
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/** Returns the number of words in the lexicon. */
	public int getWordCount() {
		return arr.size();
	}

	/** Returns the word at the specified index. */
	public String getWord(int index) {
		return arr.get(index);
	}
}