// TabooTest.java
// Taboo class tests -- nothing provided.

import org.junit.jupiter.api.Test;

import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

public class TabooTest {

    //
    // constructor
    //

    @Test
    public void testConstructorSimple() {
        List<Integer> rules1 = Arrays.asList(1, 9, 9, 1);
        Taboo<Integer> t1 = new Taboo<Integer>(rules1);

        List<Character> rules2 = Arrays.asList('a', 'b', 'c', 'd');
        Taboo<Character> t2 = new Taboo<Character>(rules2);

        List<String> rules3 = Arrays.asList("a", "c", "a", "bd");
        Taboo<String> t3 = new Taboo<String>(rules3);

        assertNotNull(t1);
        assertNotNull(t2);
        assertNotNull(t3);

    }

    @Test
    public void testConstructorEmpty() {
        List<Object> rules = Collections.emptyList();
        Taboo<Object> empty = new Taboo<Object>(rules);

        assertNotNull(empty);
    }

    //
    // noFollow
    //

    @Test
    public void testNoFollowEmpty(){
        List<Object> rules = Collections.emptyList();
        Taboo<Object> empty = new Taboo<Object>(rules);

        assertEquals(Collections.emptySet(), empty.noFollow(Collections.emptySet()));
        assertEquals(Collections.emptySet(), empty.noFollow(null));
        assertEquals(Collections.emptySet(), empty.noFollow(""));
        assertEquals(Collections.emptySet(), empty.noFollow("a"));
        assertEquals(Collections.emptySet(), empty.noFollow(1));
        assertEquals(Collections.emptySet(), empty.noFollow('a'));
        assertEquals(Collections.emptySet(), empty.noFollow(Arrays.asList(1, 9, 9, 1)));
        assertEquals(Collections.emptySet(), empty.noFollow(Arrays.asList('a', 'b', 'c', 'd')));
        assertEquals(Collections.emptySet(), empty.noFollow(Arrays.asList("a", "c", "a", "bd")));
    }

    @Test
    public void testNoFollowStringsSimple(){
        List<String> rules = Arrays.asList("a", "c", "a", "bd");
        Taboo<String> t = new Taboo<String>(rules);

        Set<String> forA = new HashSet<String>();
        Collections.addAll(forA, "c", "bd");
        Set<String> forC = new HashSet<String>();
        Collections.addAll(forC, "a");

        assertEquals(forA, t.noFollow("a"));
        assertEquals(Collections.emptySet(), t.noFollow("bd"));
        assertEquals(forC, t.noFollow("c"));
        assertEquals(Collections.emptySet(), t.noFollow("x"));
    }

    @Test
    public void testNoFollowStringsComplex() {
        List<String> rules = Arrays.asList("sadwadsa", "dwadsad c", "xcvssaa", "bwwwssxad",
                "asdawsd", "dadsxv", "adsvbgsd", "awdafdohluikbugbj,", "asn,mzvcxkhi",
                "sakldjhud", "knl.jhluo", "sn,adjkqh2l3iu", "mn,dfol2j3n4m", ".nlkajhr3;ewdm,",
                "sdawkljedsvcmx", "as", "b", "b", "b", "b", "asn,dkvhc", "21r34564", "213");

        Taboo<String> t = new Taboo<String>(rules);

        Map<String, Set<String>> strings = new HashMap<String, Set<String>>();

        for (int i = 0; i < rules.size() - 1; i++){
            if(!strings.containsKey(rules.get(i))){
                Set<String> s = new HashSet<String>();
                strings.put(rules.get(i), s);
            }
            strings.get(rules.get(i)).add(rules.get(i+1));
        }

        for (int i = 0; i < rules.size() - 1; i++) {
            assertEquals(strings.get(rules.get(i)), t.noFollow(rules.get(i)));
        }
    }

    @Test
    public void testNoFollowIntegersSimple(){
        List<Integer> rules = Arrays.asList(1, 9, 9, 1);
        Taboo<Integer> t = new Taboo<Integer>(rules);

        Set<Integer> one = new HashSet<Integer>();
        Collections.addAll(one, 9);
        Set<Integer> nine = new HashSet<Integer>();
        Collections.addAll(nine, 9, 1);

        assertEquals(one, t.noFollow(1));
        assertEquals(Collections.emptySet(), t.noFollow(0));
        assertEquals(nine, t.noFollow(9));
    }

    @Test
    public void testNoFollowIntegersComplex(){
        List<Integer> rules = Arrays.asList(1, 3, 2, 20, 13, 2 , 15, 1, 1, 7);
        Taboo<Integer> t = new Taboo<Integer>(rules);

        List<Set<Integer>> arr = new LinkedList<Set<Integer>>();

        for(int i = 0; i < 21; i++){
            Set<Integer> s = new HashSet<>();
            arr.add(s);
        }

        for(int i = 0; i < rules.size() - 1; i++)
            arr.get(rules.get(i)).add(rules.get(i+1));

        for (Integer rule : rules) {
            assertEquals(arr.get(rule), t.noFollow(rule));
        }
    }

    @Test
    public void testNoFollowArrays(){
        List<List<Integer>> rules = Arrays.asList(Arrays.asList(1, 3, 4, 8),
                                                Arrays.asList(1, 3, 2, 9));

        Taboo<List<Integer>> t = new Taboo<List<Integer>>(rules);

        Set<List<Integer>> expected = new HashSet<>();
        expected.add(Arrays.asList(1, 3, 2, 9));

        assertEquals(expected, t.noFollow(Arrays.asList(1, 3, 4, 8)));
    }

    //
    // reduce
    //

    @Test
    public void testReduceAllEmpty(){
        List<Object> rules = Collections.emptyList();

        Taboo<Object> empty1 = new Taboo<Object>(rules);
        Taboo<String> empty2 = new Taboo<String>(Arrays.asList(""));

        List<Object> before1 = Collections.emptyList();
        List<String> before2 = Arrays.asList("");

        List<Object> after1 = Collections.emptyList();
        List<String> after2 = Arrays.asList("");

        empty1.reduce(before1);
        empty2.reduce(before2);

        assertEquals(after1, before1);
        assertEquals(after2, before2);
    }

    @Test
    public void testReduceEmptyRules(){
        List<Object> rules = Collections.emptyList();
        Taboo<Object> empty = new Taboo<Object>(rules);


        List<Object> before = Arrays.asList("", "a", "b", "c", "asda");
        List<Object> after = Arrays.asList("", "a", "b", "c", "asda");

        empty.reduce(before);

        assertEquals(after, before);
    }

    @Test
    public void testReduceEmptyString(){
        List<Object> rules = Arrays.asList("a", "c", "a", "bd");
        Taboo<Object> t = new Taboo<Object>(rules);

        List<Object> before = Collections.emptyList();
        List<Object> after = Collections.emptyList();

        t.reduce(before);

        assertEquals(after, before);
    }

    @Test
    public void testReduceLengthOne(){
        List<String> rules = Arrays.asList("a", "c", "a", "bd");
        Taboo<String> t = new Taboo<String>(rules);

        List<String> before = Arrays.asList("a");
        List<String> after = Arrays.asList("a");

        t.reduce(before);

        assertEquals(after, before);
    }

    @Test
    public void testReduceIntegers(){
        List<Integer> rules = Arrays.asList(1, 3, 2, 20, 13, 2 , 15, 1, 1, 7);
        Taboo<Integer> t = new Taboo<Integer>(rules);

        List<Integer> before = new ArrayList<>();
        Collections.addAll(before, 1, 0, 14, 13, null, 15, 2, 1, 1, 7);
        List<Integer> after = Arrays.asList(1, 0, 14, 13, null, 15, 2, 1);

        t.reduce(before);

        assertEquals(after, before);
    }

    @Test
    public void testReduceStrings(){
        List<String> rules = Arrays.asList("a", "c", "a", "b");
        Taboo<String> t = new Taboo<String>(rules);

        List<String> before = new ArrayList<>();
        Collections.addAll(before, "a", "c", "b", "x", "c", "a");
        List<String> after = Arrays.asList("a", "x", "c");

        t.reduce(before);

        assertEquals(after, before);
    }

    @Test
    public void testReduceArrays(){
        List<List<Integer>> rules = Arrays.asList(Arrays.asList(1, 3, 4, 8),
                Arrays.asList(1, 3, 2, 9));
        Taboo<List<Integer>> t = new Taboo<List<Integer>>(rules);

        List<List<Integer>> before = new ArrayList<>();
        Collections.addAll(before, Arrays.asList(1, 3, 2, 9), Arrays.asList(1, 3, 4, 8),
                Arrays.asList(1, 3, 2, 9) );
        List<List<Integer>> after = Arrays.asList(Arrays.asList(1, 3, 2, 9), Arrays.asList(1, 3, 4, 8));

        t.reduce(before);

        assertEquals(after, before);

    }

    @Test
    public void testReduceNoReduce(){
        List<Integer> rules = Arrays.asList(1, 3, 2, 20, 13, 2 , 15, 1, 1, 7);
        Taboo<Integer> t = new Taboo<Integer>(rules);

        List<Integer> before = new ArrayList<>();
        Collections.addAll(before, 1, 0, 14, 13, null, 15, 2, 1);
        List<Integer> after = Arrays.asList(1, 0, 14, 13, null, 15, 2, 1);

        t.reduce(before);

        assertEquals(after, before);
    }

    @Test
    public void testReduceNullTaboo(){
        List<Integer> rules = Arrays.asList(1, 3, 2 ,null, 20, 13, null, 2 , 15, 1, 1, null, 7);
        Taboo<Integer> t = new Taboo<Integer>(rules);

        List<Integer> before = new ArrayList<>();
        Collections.addAll(before, 1, 0, 14, 13, null, 15, 2, 1, 1, 7);
        List<Integer> after = Arrays.asList(1, 0, 14, 13, null, 15, 2, 1, 7);

        t.reduce(before);

        assertEquals(after, before);
    }

    @Test
    void testReduceRepeated() {
        List<String> rules = Arrays.asList("a", "a", "a", "a");
        Taboo<String> t = new Taboo<String>(rules);

        List<String> before = new ArrayList<>();
        Collections.addAll(before, "a", "a", "a", "a", "a", "a");
        List<String> after = Arrays.asList("a");

        t.reduce(before);

        assertEquals(after, before);
    }
}
