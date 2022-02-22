
// Test cases for CharGrid -- a few basic tests are provided.

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

public class CharGridTest{

	//
	// constructor
	//

	@Test
	public void testConstructorSimple(){
		char[][] grid = new char[][] {
				{'a', 'y', ' '},
				{'x', 'a', 'z'},
		};

		CharGrid cg = new CharGrid(grid);

		assertNotNull(cg);
	}

	@Test
	public void testConstructorEmpty(){
		char[][] grid = new char[][]{

		};

		CharGrid cg = new CharGrid(grid);
		assertNotNull(cg);
	}

	//
	// charArea
	//

	@Test
	public void testCharArea1() {
		char[][] grid = new char[][] {
				{'a', 'y', ' '},
				{'x', 'a', 'z'},
			};
		
		
		CharGrid cg = new CharGrid(grid);
				
		assertEquals(4, cg.charArea('a'));
		assertEquals(1, cg.charArea('z'));
		assertEquals(0, cg.charArea('t'));
	}

	@Test
	public void testCharArea2() {
		char[][] grid = new char[][] {
				{'c', 'a', ' '},
				{'b', ' ', 'b'},
				{' ', ' ', 'a'}
			};
		
		CharGrid cg = new CharGrid(grid);
		
		assertEquals(6, cg.charArea('a'));
		assertEquals(3, cg.charArea('b'));
		assertEquals(1, cg.charArea('c'));
		assertEquals(9, cg.charArea(' '));
		assertEquals(0, cg.charArea('ჭ'));
	}

	@Test
	public void charAreaOneRow(){
		char[][] grid = new char[][]{
				{'r' , ' ', 'o', 'w', 'r'}
		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(1, cg.charArea('o'));
		assertEquals(5, cg.charArea('r'));
		assertEquals(1, cg.charArea(' '));
		assertEquals(0, cg.charArea('a'));
	}

	@Test
	public void charAreaOneColumn(){
		char[][] grid = new char[][]{
				{'c'},
				{'o'},
				{' '},
				{'l'},
				{'l'},
				{'l'}
		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(1, cg.charArea('o'));
		assertEquals(3, cg.charArea('l'));
		assertEquals(1, cg.charArea(' '));
		assertEquals(0, cg.charArea('a'));
	}

	@Test
	public void charAreaOneXOne(){
		char[][] grid = new char[][]{
				{'1'}
		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(1, cg.charArea('1'));
		assertEquals(0, cg.charArea('0'));
	}

	@Test
	public void charAreaEmptyGrid(){
		char[][] grid = new char[][]{

		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(0, cg.charArea('1'));
		assertEquals(0, cg.charArea('0'));
	}

	@Test
	public void charAreaJaggedArrays1(){
		char[][] grid = new  char[][]{
				{'1'},
				{'b', 'o', 'c'},
				{'d', 'a' , 'o'},
				{'b', 'b'}
		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(1, cg.charArea('1'));
		assertEquals(4, cg.charArea('o'));
		assertEquals(6, cg.charArea('b'));
		assertEquals(0, cg.charArea('ჭ'));
	}

	@Test
	public void charAreajaggedArrays2(){
		char[][] grid = new  char[][]{
				{'1', 't', 'k', 'g', 'b'},
				{'b', 'o', 'c'},
				{'d', 'a' , 'o'},
				{'b', 'b'}
		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(1, cg.charArea('g'));
		assertEquals(4, cg.charArea('o'));
		assertEquals(20, cg.charArea('b'));
		assertEquals(0, cg.charArea('ჭ'));
	}


	//
	// countPlus
	//

	@Test
	public void countPlusNormal(){
		char[][] grid = new char[][]{
				{' ', ' ', 'p', ' ', ' ', ' ', ' ', ' ', ' '},
				{' ', ' ', 'p', ' ', ' ', ' ', ' ', 'x', ' '},
				{'p', 'p', 'p', 'p', 'p', ' ', 'x', 'x', 'x'},
				{' ', ' ', 'p', ' ', ' ', 'y', ' ', 'x', ' '},
				{' ', ' ', 'p', ' ', 'y', 'y', 'y', ' ', ' '},
				{'z', 'z', 'z', 'z', 'z', 'y', 'z', 'z', 'z'},
				{' ', ' ', 'x', 'x', ' ', 'y', ' ', ' ', ' '},

		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(2, cg.countPlus());
	}

	@Test
	public void countPlusJagged(){
		char[][] grid = new char[][]{
				{' ', ' ', 'p'},
				{' ', ' ', 'p', ' ', ' ', ' ', ' ', 'x'},
				{'p', 'p', 'p', 'p', 'p', ' ', 'x', 'x', 'x'},
				{' ', ' ', 'p', ' ', ' ', 'y', ' ', 'x'},
				{' ', ' ', 'p', ' ', 'y', 'y', 'y'},
				{'z', 'z', 'z', 'z', 'z', 'y', 'z', 'z', 'z'},
				{' ', ' ', 'x', 'x', ' '},

		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(3, cg.countPlus());
	}

	@Test
	public void countPlusPlus(){
		char[][] grid = new char[][]{
				{' ', '+'},
				{'+', '+', '+'},
				{' ', '+'}
		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(1, cg.countPlus());
	}

	@Test
	public void countPlusEmpty(){
		char[][] grid = new char[][]{

		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(0, cg.countPlus());
	}

	@Test
	public void countPlusOneColumn(){
		char[][] grid = new char[][]{
				{'a'},
				{'a'},
				{'a'},
				{'a'},
				{'a'},
				{'a'}
		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(0, cg.countPlus());
	}

	@Test
	public void countPlusOneRow(){
		char[][] grid = new char[][]{
				{'a' , 'a', 'a', 'a', 'a'}
		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(0, cg.countPlus());
	}

	@Test
	public void countPlusOneXOne(){
		char[][] grid = new char[][]{
				{'+'}
		};

		CharGrid cg = new CharGrid(grid);

		assertEquals(0, cg.countPlus());
	}
}
