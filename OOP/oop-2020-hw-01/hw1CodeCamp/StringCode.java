import java.util.HashSet;
import java.util.Set;

// CS108 HW1 -- String static methods

public class StringCode {

	/**
	 * Given a string, returns the length of the largest run.
	 * A a run is a series of adajcent chars that are the same.
	 * @param str
	 * @return max run length
	 */
	public static int maxRun(String str) {
		if(str.length() == 0) return 0;

		int maxRun = 1;
		char curr = str.charAt(0);
		int currRun = 1;

		for(int i = 1; i < str.length(); i++){
			if(str.charAt(i) == curr) currRun++;
			else {
				if(currRun > maxRun) maxRun = currRun;
				currRun = 1;
				curr = str.charAt(i);
			}
		}

		return Math.max(currRun, maxRun); //last run might be maximal
	}

	
	/**
	 * Given a string, for each digit in the original string,
	 * replaces the digit with that many occurrences of the character
	 * following. So the string "a3tx2z" yields "attttxzzz".
	 * @param str string
	 * @return blown up string
	 */
	public static String blowup(String str) {
		if(str.length() == 0) return "";
		
		StringBuilder res = new StringBuilder();
		int length = str.length();
		if(Character.isDigit(str.charAt(str.length()-1))) length--;

		for(int i = 0; i < length; i++){
			char ch = str.charAt(i);
			if(Character.isDigit(ch))
				//makes string out of next char repeated ch's value times and appends it to main string
				res.append(Character.toString(str.charAt(i + 1)).repeat(ch - '0'));

			else res.append(ch);
		}

		return res.toString();
	}
	
	/**
	 * Given 2 strings, consider all the substrings within them
	 * of length len. Returns true if there are any such substrings
	 * which appear in both strings.
	 * Compute this in linear time using a HashSet. Len will be 1 or more.
	 */
	public static boolean stringIntersect(String a, String b, int len) {
		Set<String> substrings = new HashSet<String>();

		for(int i = 0; i < a.length() - len + 1; i++){
			substrings.add(a.substring(i, i + len));
		}

		for(int i = 0; i < b.length() - len + 1; i++){
			if(substrings.contains(b.substring(i, i + len))) return true;
		}

		return false;
	}
}
