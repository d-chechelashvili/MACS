import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

public class TetrisGridTest {

	//
	//constructor
	//
	@Test
	public void testConstructorSimple(){
		boolean[][] grid = {
				{false, true, false},
				{false, true, true}
		};

		TetrisGrid tetris = new TetrisGrid(grid);
		assertNotNull(tetris);
	}

	@Test
	public void testConstructorChange(){
		//constructor description says that it doesn't make a copy
		//so tests if constructor has a bug which changes passed grid

		boolean[][] grid = {
				{false, true, false},
				{false, true, true}
		};

		boolean[][] gridCopy = grid.clone();

		TetrisGrid tetris = new TetrisGrid(grid);

		assertTrue(Arrays.deepEquals(grid, gridCopy));
	}

	//
	//getGrid
	//
	@Test
	public void testGetterNormal(){
		boolean[][] grid = {
				{true, false, false},
				{false, true, true}
		};

		TetrisGrid tetris = new TetrisGrid(grid);

		assertTrue(Arrays.deepEquals(grid, tetris.getGrid()));
	}

	@Test
	public void testGetterOneXOne(){
		boolean[][] grid1 = {
				{true}
		};

		boolean[][] grid2 = {
				{false}
		};

		TetrisGrid tetris1 = new TetrisGrid(grid1);
		TetrisGrid tetris2 = new TetrisGrid(grid2);

		assertTrue(Arrays.deepEquals(grid1, tetris1.getGrid()));
		assertTrue(Arrays.deepEquals(grid2, tetris2.getGrid()));
	}



	//
	//clearRows
	//

	// Provided simple clearRows() test
	// width 2, height 3 grid

	@Test
	public void testClearSimple() {
		boolean[][] before =
		{
			{true, true, false},
			{false, true, true}
		};
		
		boolean[][] after =
		{	
			{true, false, false},
			{false, true, false}
		};
		
		TetrisGrid tetris = new TetrisGrid(before);
		tetris.clearRows();

		assertTrue(Arrays.deepEquals(after, tetris.getGrid()));
	}

	@Test
	public void testClearComplex() {
		boolean[][] before = {
				{true, true, false, true, true, false},
				{false, true, false, true, true, true},
				{false, true, true, true, true, false}
		};

		boolean[][] after ={
				{true,  false,  false,  false,  false, false},
				{false,  false, true, false,  false, false},
				{false, true, false, false,  false,  false}
		};


		TetrisGrid tetris = new TetrisGrid(before);
		tetris.clearRows();

		assertTrue(Arrays.deepEquals(after, tetris.getGrid()));
	}

	@Test
	public void testClearOneXOne(){
		boolean[][] before1 = {
				{true}
		};

		boolean[][] after1 = {
				{false}
		};

		boolean[][] before2 = {
				{false}
		};

		boolean[][] after2 = {
				{false}
		};

		TetrisGrid tetris1 = new TetrisGrid(before1);
		tetris1.clearRows();

		TetrisGrid tetris2 = new TetrisGrid(before2);
		tetris2.clearRows();

		assertTrue(Arrays.deepEquals(after1, tetris1.getGrid()));
		assertTrue(Arrays.deepEquals(after2, tetris2.getGrid()));
	}

	@Test
	public void testClearAlreadyClear(){
		boolean[][] before = {
				{true, true, false, true, false, false},
				{false, false, false, false, true, true},
				{false, true, true, true, true, false}
		};

		boolean[][] after = {
				{true, true, false, true, false, false},
				{false, false, false, false, true, true},
				{false, true, true, true, true, false}
		};

		TetrisGrid tetris = new TetrisGrid(before);
		tetris.clearRows();

		assertTrue(Arrays.deepEquals(after, tetris.getGrid()));
	}

	@Test
	public void testClearAllTrue(){
		boolean[][] before = {
				{true, true, true, true, true, true, true},
				{true, true, true, true, true, true, true},
				{true, true, true, true, true, true, true},
				{true, true, true, true, true, true, true},
				{true, true, true, true, true, true, true}
		};

		boolean[][] after = {
				{false, false, false, false, false, false, false},
				{false, false, false, false, false, false, false},
				{false, false, false, false, false, false, false},
				{false, false, false, false, false, false, false},
				{false, false, false, false, false, false, false}
		};

		TetrisGrid tetris = new TetrisGrid(before);
		tetris.clearRows();

		assertTrue(Arrays.deepEquals(after, tetris.getGrid()));
	}
	
}
