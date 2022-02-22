import java.util.*;

public class Appearances {
	
	/**
	 * Returns the number of elements that appear the same number
	 * of times in both collections. Static method. (see handout).
	 * @return number of same-appearance elements
	 */
	public static <T> int sameCount(Collection<T> a, Collection<T> b) {
		Map<T, Integer> data = new HashMap<T, Integer>();
		for(T elem: a)
			data.merge(elem, 1, Integer::sum); //short syntax for incrementing value by 1

		for(T elem: b)
			data.merge(elem, -1, Integer::sum);

		int result = 0;
		for(T key : data.keySet())
			if(data.get(key) == 0) result++;

		return result;
	}
	
}
