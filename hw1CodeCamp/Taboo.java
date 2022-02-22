
/*
 HW1 Taboo problem class.
 Taboo encapsulates some rules about what objects
 may not follow other objects.
 (See handout).
*/

import java.util.*;

public class Taboo<T> {
	private Map<T, Set <T> > taboos;
	
	/**
	 * Constructs a new Taboo using the given rules (see handout.)
	 * @param rules rules for new Taboo
	 */

	//doesn't make copies of objects in rules list
	public Taboo(List<T> rules) {
		taboos = new HashMap<>();
		for(int i = 0; i < rules.size() - 1; i++){
			T elem = rules.get(i);
			T nextElem = rules.get(i+1);
			if(elem == null || nextElem == null) continue;
			if(!taboos.containsKey(elem)){
				Set<T> s = new HashSet<>();
				taboos.put(elem, s);
			}
			taboos.get(elem).add(nextElem);
		}
	}
	
	/**
	 * Returns the set of elements which should not follow
	 * the given element.
	 * @param elem
	 * @return elements which should not follow the given element
	 */

	//makes copies of rules so that user can't change rules(accidentaly) that were active
	// 													by the time user called this function
	public Set<T> noFollow(T elem) {
		Set<T> copy = taboos.getOrDefault(elem, Collections.emptySet());
		return copy;
	}
	
	/**
	 * Removes elements from the given list that
	 * violate the rules (see handout).
	 * @param list collection to reduce
	 */
	public void reduce(List<T> list) {
		for(int i = 0; i < list.size()-1; i++){
			T elem = list.get(i);
			T nextElem = list.get(i+1);
			if(taboos.getOrDefault(elem, Collections.emptySet()).contains(nextElem)){
				list.remove(i+1);
				i--;
			}
		}
	}
}
