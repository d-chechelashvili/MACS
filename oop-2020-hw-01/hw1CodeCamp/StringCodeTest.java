// StringCodeTest
// Some test code is provided for the early HW1 problems,
// and much is left for you to add.

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class StringCodeTest {
	//
	// blowup
	//

	@Test
	public void testBlowup1() {
		// basic cases
		assertEquals("xxaaaabb", StringCode.blowup("xx3abb"));
		assertEquals("xxxZZZZ", StringCode.blowup("2x3Z"));
	}

	@Test
	public void testBlowup2() {
		// things with digits
		
		// digit at end
		assertEquals("axxx", StringCode.blowup("a2x3"));
		
		// digits next to each other
		assertEquals("a33111", StringCode.blowup("a231"));
		
		// try a 0
		assertEquals("aabb", StringCode.blowup("aa0bb"));
	}

	@Test
	public void testBlowup3() {
		// weird chars, empty string
		assertEquals("AB&&,- ab", StringCode.blowup("AB&&,- ab"));
		assertEquals("", StringCode.blowup(""));
		
		// string with only digits
		assertEquals("", StringCode.blowup("2"));
		assertEquals("33", StringCode.blowup("23"));
	}

	@Test
	public void blowupLengthOne(){
		assertEquals("", StringCode.blowup("1"));
		assertEquals("a", StringCode.blowup("a"));
		assertEquals("5", StringCode.blowup("15"));
	}

	@Test
	public void blowupContainsZero(){
		assertEquals("", StringCode.blowup("0"));
		assertEquals("a", StringCode.blowup("0a"));
		assertEquals("", StringCode.blowup("00"));
		assertEquals("", StringCode.blowup("000"));
		assertEquals("5", StringCode.blowup("015"));
		assertEquals("00", StringCode.blowup("0020"));
	}
	
	
	//
	// maxRun
	//

	@Test
	public void testRun1() {
		assertEquals(2, StringCode.maxRun("hoopla"));
		assertEquals(3, StringCode.maxRun("hoopllla"));
	}

	@Test
	public void testRun2() {
		assertEquals(3, StringCode.maxRun("abbcccddbbbxx"));
		assertEquals(0, StringCode.maxRun(""));
		assertEquals(3, StringCode.maxRun("hhhooppoo"));
	}

	@Test
	public void testRun3() {
		// "evolve" technique -- make a series of test cases
		// where each is change from the one above.
		assertEquals(1, StringCode.maxRun("123"));
		assertEquals(2, StringCode.maxRun("1223"));
		assertEquals(2, StringCode.maxRun("112233"));
		assertEquals(3, StringCode.maxRun("1112233"));
	}

	@Test
	public void maxRunLengthOne(){
		assertEquals(1, StringCode.maxRun("a"));
		assertEquals(1, StringCode.maxRun("1"));
		assertEquals(1, StringCode.maxRun(" "));
	}

	@Test
	public void maxRunBigString(){
		assertEquals(80, StringCode.maxRun(("adsadcxjkkuiq2143jk253ntrwm,enncxkikadsjncx" +
				"asmndajsfauisekjdsffffffffffffffffffffffsadASEOIDCXZM.Asldjiawjsdk,z" +
				"sakjddddddddddddddddddddaweeeeeeeeewafdsmnvxc,,,,mmmmmmmmmmmmmmmmmmawidasjncmx" +
				"kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk")));
	}

	@Test
	public void maxRunFullBigRun(){
		assertEquals(245, StringCode.maxRun("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" +
				"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" +
				"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" +
				"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"));

	}

	//
	//stringIntersect
	//

	@Test
	public void intersectLengthZero(){
		assertFalse(StringCode.stringIntersect("", "a", 1));
		assertFalse(StringCode.stringIntersect("", "dbc", 2));
		assertFalse(StringCode.stringIntersect("gcc", "", 3));
		assertFalse(StringCode.stringIntersect("", "", 1));
	}

	@Test
	public void intersectLengthOne(){
		assertFalse(StringCode.stringIntersect("a", "b", 1));
		assertFalse(StringCode.stringIntersect("a", "z", 1));
		assertTrue(StringCode.stringIntersect("a", "a", 1));
		assertTrue(StringCode.stringIntersect("g", "g", 1));
		assertFalse(StringCode.stringIntersect("h", "q", 2));
		assertFalse(StringCode.stringIntersect("h", "q", 10));
	}

	@Test
	public void intersectNormalStrings(){
		assertTrue(StringCode.stringIntersect("Tom Marvolo Riddle", "I Am Lord Voldemort", 2));
		assertTrue(StringCode.stringIntersect("Tiste Andii", "Tiste Edur", 6));
		assertFalse(StringCode.stringIntersect("Kvothe", "Ambrose", 2));
		assertFalse(StringCode.stringIntersect("Kvothe", "Ambrose", 7));
		assertFalse(StringCode.stringIntersect("Logen Ninefingers", "Bloody-Nine", 10));
		assertFalse(StringCode.stringIntersect("Geralt of Rivia", "Butcher of Blaviken", 7));
		assertFalse(StringCode.stringIntersect("Bloody Rose", "Golden Gabe", 5));
	}

	@Test
	public void intersectNumbers(){
		assertTrue(StringCode.stringIntersect("123333", "333", 3));
		assertFalse(StringCode.stringIntersect("29.03.2020", "......", 2));
		assertTrue(StringCode.stringIntersect("29.03.2020", "2020", 4));
		assertFalse(StringCode.stringIntersect("29.03.2020", "28.03.2020", 11));
	}

	@Test
	public void intersectSame(){
		assertTrue(StringCode.stringIntersect("Whiskeyjack", "Whiskeyjack", 4));
		assertTrue(StringCode.stringIntersect("Whiskeyjack", "Whiskeyjack", 11));
		assertFalse(StringCode.stringIntersect("ITACHI", "Itachi", 4));
		assertFalse(StringCode.stringIntersect("ITACHI", "Itachi", 4));
		assertFalse(StringCode.stringIntersect("Itachi", "Itachi", 7));
	}
}
