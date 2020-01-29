
import acm.util.*;
import java.io.*;
import java.util.*;

public class HangmanLexiconExtension1 {
	private ArrayList<String> arr3 = new ArrayList<String>();

	// This is the HangmanLexicon constructor
	public HangmanLexiconExtension1() {
		try {
			BufferedReader rd = new BufferedReader(new FileReader("gameList.txt"));
			while (true) {
				String s = rd.readLine();
				if (s == null)
					break;
				arr3.add(s);
			}
			rd.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/** Returns the number of words in the lexicon. */
	public int getWordCounti() {
		return arr3.size();
	}

	/** Returns the word at the specified index. */
	public String getWordi(int index) {
		String s = arr3.get(index);
		return s;
	}

}