import org.junit.jupiter.api.Test;

import java.lang.reflect.Array;
import java.util.*;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class AppearancesTest{
	// utility -- converts a string to a list with one
	// elem for each char.
	private List<String> stringToList(String s) {
		List<String> list = new ArrayList<String>();
		for (int i=0; i<s.length(); i++) {
			list.add(String.valueOf(s.charAt(i)));
			// note: String.valueOf() converts lots of things to string form
		}
		return list;
	}

	@Test
	public void testSameCount1() {
		List<String> a = stringToList("abbccc");
		List<String> b = stringToList("cccbba");
		assertEquals(3, Appearances.sameCount(a, b));
	}

	@Test
	public void testSameCount2() {
		// basic List<Integer> cases
		List<Integer> a = Arrays.asList(1, 2, 3, 1, 2, 3, 5);
		assertEquals(1, Appearances.sameCount(a, Arrays.asList(1, 9, 9, 1)));
		assertEquals(2, Appearances.sameCount(a, Arrays.asList(1, 3, 3, 1)));
		assertEquals(1, Appearances.sameCount(a, Arrays.asList(1, 3, 3, 1, 1)));
	}
	
	@Test
	public void testSameCount3(){
		Set<Integer> s1 = new HashSet<>();
		for(int i = 0 ; i < 90; i+=2){
			s1.add(3*i*i - 5*i + 13);
		}

		assertEquals(8, Appearances.sameCount(s1, Arrays.asList(22805, 18813, 9141, 2925, 531, 263, 15, 13)));
		assertEquals(3, Appearances.sameCount(s1, Arrays.asList(11235, 12751, 8491, 1321, 896)));
	}

	@Test
	public void testSameCountEmpty(){
		ArrayList<String> arr = new ArrayList<>();
		assertEquals(0, arr.size());
		List<String> ls = new LinkedList<>();
		assertEquals(0, ls.size());

		assertEquals(0, Appearances.sameCount(arr, ls));
		assertEquals(0, Appearances.sameCount(ls, Arrays.asList("asvsd", "sadasa", "sadds", "sdawas", "kbsc")));
	}

	@Test
	public void testSameCountFullMatch(){
		int repeatNum = 4;
		Queue<Integer> q = new PriorityQueue<>();
		Stack<Integer> s = new Stack<>();
		for(int i = 0; i < 300; i+=3){
			int elem = 5*i*i*i - 3*i*i + 14*i + 27;
			for(int j = 0; j < repeatNum; j++){
				q.add(elem);
				s.add(elem);
			}
		}
		assertEquals(q.size(), s.size());

		assertEquals(q.size() / repeatNum, Appearances.sameCount(q, s));
	}

	@Test
	public void testSameCountNoMatch(){
		Queue<Integer> q = new PriorityQueue<>();
		Stack<Integer> s = new Stack<>();
		for(int i = 0; i < 300; i+=3){
			int elem = 5*i*i*i - 3*i*i + 14*i + 27;
			q.add(elem);
			s.add(elem+1);
		}
		assertEquals(q.size(), s.size());

		assertEquals(0, Appearances.sameCount(q, s));
	}
}































