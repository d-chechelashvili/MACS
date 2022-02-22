import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;


public class BoardTest{

	private static Piece[] pieces;

	@BeforeAll
	static void init(){
		pieces = Piece.getPieces();
	}

	//
	// constructor
	//

	@Test
	public void testConstructorSimple(){
		Board b = new Board(3, 6);
		Board b1 = new Board(0, 0);
		Board b2 = new Board(1, 0);
		Board b3 = new Board(0, 1);
	}

	//
	//getters
	//

	@Test
	public void testGetHeight(){
		Board b = new Board(3, 6);
		assertEquals(6 , b.getHeight());

		Board b1 = new Board(0, 3);
		assertEquals(3 , b1.getHeight());

		Board b2 = new Board(3, 0);
		assertEquals(0 , b2.getHeight());

		Board b3 = new Board(0, 0);
		assertEquals(0 , b3.getHeight());
	}

	@Test
	public void testGetWidth(){
		Board b = new Board(6, 3);
		assertEquals(6 , b.getWidth());

		Board b1 = new Board(3, 0);
		assertEquals(3 , b1.getWidth());

		Board b2 = new Board(0, 3);
		assertEquals(0 , b2.getWidth());

		Board b3 = new Board(0, 0);
		assertEquals(0 , b3.getWidth());
	}


	//
	// getGrid
	//

	@Test
	public void testGetGridFalse(){
		Board b = new Board(7,9);

		for(int i = 0 ; i < b.getWidth(); i++){
			for (int j = 0; j < b.getHeight(); j++)
				assertFalse(b.getGrid(i , j));
		}
	}


	@Test
	public void testGetGridInvalid(){
		Board b = new Board(7,9);

		assertTrue(b.getGrid(-1, 5));
		assertTrue(b.getGrid(1, -5));
		assertTrue(b.getGrid(-2, -3));
		assertTrue(b.getGrid(8, 3));
		assertTrue(b.getGrid(7, 3));
		assertTrue(b.getGrid(3, 10));
		assertTrue(b.getGrid(3, 9));
		assertTrue(b.getGrid(10, 10));
		assertTrue(b.getGrid(7, 9));
	}

	@Test
	public void testGetGridTrue(){
		Board b = new Board(8, 8);


		b.place(pieces[Piece.STICK], 0, 0 );
		b.commit();
		for(int i = 0; i < pieces[Piece.STICK].getHeight(); i++)
			assertTrue(b.getGrid(0, i));

		String hundredPercentCoverage = b.toString();
	}

	//
	// place
	//

	@Test
	public void testPlaceOk(){
		Board b = new Board(7, 9);

		int result = b.place(pieces[Piece.PYRAMID],3,5);
		b.commit();
		assertEquals(Board.PLACE_OK, result);
		result = b.place(pieces[Piece.STICK].fastRotation(), 0, 0);
		b.commit();
		assertEquals(Board.PLACE_OK, result);
		result = b.place(pieces[Piece.STICK].fastRotation().fastRotation(), 6, 0);
		b.commit();
		assertEquals(Board.PLACE_OK, result);
	}

	@Test
	public void testPlaceFillBoard(){
		Board b = new Board(3, 4);

		getFilledBoard3X4(b);

		for (int i = 0; i < b.getWidth(); i++) {
			for (int j = 0; j < b.getHeight(); j++)
				assertTrue(b.getGrid(i,j));
		}

	}

	private void getFilledBoard3X4(Board b){
		int result = b.place(pieces[Piece.PYRAMID], 0, 0);
		b.commit();
		assertEquals(Board.PLACE_ROW_FILLED, result);

		Piece fourthPyramidRot = pieces[Piece.PYRAMID].fastRotation().
				fastRotation().fastRotation();
		result = b.place(fourthPyramidRot, 0, 1);
		b.commit();
		assertEquals(Board.PLACE_OK, result);

		Piece thirdL1Rot = pieces[Piece.L1].fastRotation().fastRotation();
		result = b.place(thirdL1Rot, 1, 1);
		b.commit();
		assertEquals(Board.PLACE_ROW_FILLED, result);
	}

	@Test
	public void testPlaceBad(){
		Board b = new Board(7, 9);

		int result = b.place(pieces[Piece.STICK], 0, 0 );
		b.commit();
		assertEquals(Board.PLACE_OK, result);

		result = b.place(pieces[Piece.STICK], 0, 0 );
		b.commit();
		assertEquals(Board.PLACE_BAD, result);

		result = b.place(pieces[Piece.L1], 0, 0);
		b.commit();
		assertEquals(Board.PLACE_BAD, result);

		result = b.place(pieces[Piece.PYRAMID], 0, 2);
		b.commit();
		assertEquals(Board.PLACE_BAD, result);
	}

	@Test
	public void testPlaceOutBounds(){
		Board b = new Board(7, 9);

		int result = b.place(pieces[Piece.STICK], -1, 0 );
		b.commit();
		assertEquals(Board.PLACE_OUT_BOUNDS, result);

		result = b.place(pieces[Piece.SQUARE], 7, 9 );
		b.commit();
		assertEquals(Board.PLACE_OUT_BOUNDS, result);

		result = b.place(pieces[Piece.L2].fastRotation().fastRotation().fastRotation(),
				5, 0);
		b.commit();
		assertEquals(Board.PLACE_OUT_BOUNDS, result);

		result = b.place(pieces[Piece.PYRAMID], 6, 1);
		b.commit();
		assertEquals(Board.PLACE_OUT_BOUNDS, result);
	}

	@Test
	public void testPlaceFillRow(){
		Board b = new Board(7, 12);

		int result = b.place(pieces[Piece.STICK].fastRotation(), 0 , 0);
		b.commit();
		assertEquals(Board.PLACE_OK, result);

		result = b.place(pieces[Piece.PYRAMID], 4, 0);
		b.commit();
		assertEquals(Board.PLACE_ROW_FILLED, result);

		result = b.place(pieces[Piece.STICK].fastRotation(), 0 , 1);
		b.commit();
		assertEquals(Board.PLACE_OK, result);

		result = b.place(pieces[Piece.STICK], 4, 1);
		b.commit();
		assertEquals(Board.PLACE_OK, result);

		result = b.place(pieces[Piece.STICK], 6, 1);
		b.commit();
		assertEquals(Board.PLACE_ROW_FILLED, result);
	}

	@Test
	public void testPlaceNotCommited(){
		Board b = new Board(7, 12);

		int result = b.place(pieces[Piece.STICK].fastRotation(), 0 , 0);
		assertEquals(Board.PLACE_OK, result);

		assertThrows(RuntimeException.class,
				()-> b.place(pieces[Piece.STICK].fastRotation(), 0 , 0));
	}


	//
	// getColumnHeight
	//

	@Test
	public void testGetColumnHeightEmpty(){
		Board b = new Board(3, 0);
		for (int i = 0; i < b.getWidth(); i++)
			assertEquals(0, b.getColumnHeight(i));


		Board b1 = new Board(3, 3);
		for (int i = 0; i < b1.getWidth(); i++)
			assertEquals(0, b1.getColumnHeight(i));

	}

	@Test
	public void testGetColumnHeightNormal(){
		Board b = new Board( 7, 9);

		b.place(pieces[Piece.SQUARE], 0, 0);
		b.commit();
		assertEquals(2, b.getColumnHeight(0));
		assertEquals(2, b.getColumnHeight(1));

		b.place(pieces[Piece.STICK], 3, 1);
		b.commit();
		assertEquals(5, b.getColumnHeight(3));

		b.place(pieces[Piece.PYRAMID].fastRotation(), 5, 3);
		b.commit();
		assertEquals(5, b.getColumnHeight(5));
		assertEquals(6, b.getColumnHeight(6));


		b.place(pieces[Piece.S1], 0, 6);
		b.commit();
		assertEquals(7, b.getColumnHeight(0));
		assertEquals(8, b.getColumnHeight(1));
		assertEquals(8, b.getColumnHeight(2));
	}



	@Test
	public void testGetColumnHeightFull(){
		Board b = new Board(3, 4);

		getFilledBoard3X4(b);

		for (int i = 0; i < b.getWidth(); i++)
			assertEquals(4, b.getColumnHeight(i));
	}


	//
	// getColumnHeight
	//

	@Test
	public void testGetRowWidthEmpty(){
		Board b = new Board(0, 3);
		for (int i = 0; i < b.getHeight(); i++)
			assertEquals(0, b.getRowWidth(i));


		Board b1 = new Board(3, 3);
		for (int i = 0; i < b1.getHeight(); i++)
			assertEquals(0, b1.getRowWidth(i));

	}

	@Test
	public void testGetRowWidthNormal(){
		Board b = new Board( 7, 9);

		b.place(pieces[Piece.SQUARE], 0, 0);
		b.commit();
		assertEquals(2, b.getRowWidth(0));
		assertEquals(2, b.getRowWidth(1));

		b.place(pieces[Piece.STICK], 3, 1);
		b.commit();
		assertEquals(3, b.getRowWidth(1));
		assertEquals(1, b.getRowWidth(2));
		assertEquals(1, b.getRowWidth(3));
		assertEquals(1, b.getRowWidth(4));

		b.place(pieces[Piece.PYRAMID].fastRotation(), 5, 3);
		b.commit();
		assertEquals(2, b.getRowWidth(3));
		assertEquals(3, b.getRowWidth(4));
		assertEquals(1, b.getRowWidth(5));


		b.place(pieces[Piece.S1], 0, 6);
		b.commit();
		assertEquals(2, b.getRowWidth(6));
		assertEquals(2, b.getRowWidth(7));
	}



	@Test
	public void testGetRowWidthFull(){
		Board b = new Board(3, 4);

		getFilledBoard3X4(b);

		for (int i = 0; i < b.getHeight(); i++)
			assertEquals(3, b.getRowWidth(i));
	}

	//
	// maxHeight
	//

	@Test
	public void testMaxHeightEmpty(){
		Board b = new Board(0, 0);
		assertEquals(0, b.getMaxHeight());

		Board b1 = new Board(0, 3);
		assertEquals(0, b1.getMaxHeight());

		Board b2 = new Board(3, 0);
		assertEquals(0, b2.getMaxHeight());
	}

	@Test
	public void testMaxHeightNormal(){
		Board b = new Board(15, 17);
		assertEquals(0, b.getMaxHeight());

		b.place(pieces[Piece.STICK], 0, 0);
		b.commit();
		assertEquals(4, b.getMaxHeight());

		b.place(pieces[Piece.STICK], 0, 4);
		b.commit();
		assertEquals(8, b.getMaxHeight());

		for (int i = 0; i < 4; i++){
			b.place(pieces[Piece.STICK], 1, 4 * i);
			b.commit();
		}


		assertEquals(16, b.getMaxHeight());

	}

	@Test
	public void testMaxHeightFull(){
		Board b = new Board(3, 4);

		getFilledBoard3X4(b);

		assertEquals(4, b.getMaxHeight());
	}

	//
	// dropHeight
	//

	@Test
	public void dropHeightEmpty(){
		Board b = new Board(15,15);
		assertEquals(0, b.dropHeight(pieces[Piece.STICK], 0));

		Board b1 = new Board(2, 0);
		assertEquals(0, b1.dropHeight(pieces[Piece.STICK], 0));
	}

	@Test
	public void dropHeightNormal(){
		Board b = new Board(15, 15);

		b.place(pieces[Piece.STICK], 3, 4);
		b.commit();
		b.place(pieces[Piece.STICK], 4, 4);
		b.commit();

		assertEquals(8, b.dropHeight(pieces[Piece.S1].fastRotation(), 3));
		assertEquals(6,
				b.dropHeight(pieces[Piece.L2].fastRotation().fastRotation(), 2) );

		b.place(pieces[Piece.PYRAMID], 10, 7);
		b.commit();
		assertEquals(9,
				b.dropHeight(pieces[Piece.PYRAMID].fastRotation().fastRotation(), 10));

		b.place(pieces[Piece.L2], 7, 0);
		b.commit();
		assertEquals(1,
				b.dropHeight(pieces[Piece.L2].fastRotation().fastRotation(),7));

		b.place(pieces[Piece.STICK].fastRotation(), 10, 12);
		b.commit();
		assertEquals(13, b.dropHeight(pieces[Piece.STICK], 10));
	}

	@Test
	public void dropHeightFull(){
		Board b = new Board(3, 4);

		getFilledBoard3X4(b);

		for(int i = 0; i < b.getWidth(); i++)
			assertEquals(b.getHeight(), b.dropHeight(pieces[Piece.STICK], 0));
	}

	//
	// clearRows
	//

	@Test
	public void testClearRowsEmpty(){
		Board oldBoard = new Board(15, 15);
		Board newBoard = new Board(15, 15);

		assertEquals(0, newBoard.clearRows());
		newBoard.commit();

		checkBoardSimilarity(oldBoard, newBoard);
	}


	@Test
	public void testClearRowsSimple(){
		Board b = new Board(3, 4);

		b.place(pieces[Piece.STICK], 0, 0);
		b.commit();

		assertEquals(0, b.clearRows());
	}


	@Test
	public void testClearRowsFull(){
		Board emptyBoard = new Board(3, 4);
		Board fullBoard = new Board(3, 4);

		getFilledBoard3X4(fullBoard);

		assertEquals(4, fullBoard.clearRows());
		fullBoard.commit();

		checkBoardSimilarity(emptyBoard, fullBoard);

	}


	@Test
	public void testClearRowsNormal(){
		Board b = new Board(5, 6);

		b.place(pieces[Piece.STICK].fastRotation(), 0, 0);
		b.commit();
		b.place(pieces[Piece.STICK].fastRotation(), 0, 1);
		b.commit();
		b.place(pieces[Piece.STICK].fastRotation(), 0, 2);
		b.commit();

		assertEquals(0, b.clearRows());
		b.commit();

		b.place(pieces[Piece.STICK], 4, 0);
		b.commit();

		assertEquals(3, b.clearRows());
		b.commit();

		b.place(pieces[Piece.STICK].fastRotation(), 0, 0);
		b.commit();

		assertEquals(1, b.clearRows());
	}

	@Test
	public void testClearRowsOnlyTopFilled(){
		Board testBoard = new Board(4, 4);
		Board emptyBoard = new Board(4, 4 );

		testBoard.place(pieces[Piece.STICK].fastRotation(), 0, 3);
		testBoard.commit();

		testBoard.clearRows();
		testBoard.commit();

		checkBoardSimilarity(emptyBoard, testBoard);
	}


	@Test
	public void testClearRowsWithGap(){
		Board testBoard = new Board(4, 4);
		Board answerBoard = new Board(4, 4 );

		testBoard.place(pieces[Piece.PYRAMID], 0, 0);
		testBoard.commit();
		answerBoard.place(pieces[Piece.PYRAMID], 0, 0);
		answerBoard.commit();

		testBoard.place(pieces[Piece.STICK].fastRotation(), 0, 3);
		testBoard.commit();

		testBoard.clearRows();
		testBoard.commit();

		checkBoardSimilarity(answerBoard, testBoard);
	}


	//
	// undo
	//

	@Test
	public void testUndoCommited(){
		Board b = new Board(5, 6);
		b.place(pieces[Piece.STICK], 0, 0);
		b.commit();
		b.undo();
		b.clearRows();
		b.commit();
		b.undo();
		assertEquals(4, b.getColumnHeight(0));
		for(int i = 0; i < b.getColumnHeight(0); i++)
			assertEquals(1, b.getRowWidth(i));
		assertEquals(4, b.getMaxHeight());
	}

	@Test
	public void testDoubleUndo(){
		Board testBoard = new Board(5, 6);
		Board answerBoard = new Board(5, 6);

		testBoard.place(pieces[Piece.PYRAMID], 0, 0);
		testBoard.undo();
		testBoard.undo();

		checkBoardSimilarity(answerBoard, testBoard);

		testBoard.place(pieces[Piece.PYRAMID], 0, 0);
		answerBoard.place(pieces[Piece.PYRAMID], 0, 0);
		testBoard.commit();
		answerBoard.commit();
		testBoard.place(pieces[Piece.PYRAMID], 0, 2);
		testBoard.undo();
		testBoard.undo();

		checkBoardSimilarity(answerBoard, testBoard);
	}


	@Test
	public void testUndoPlace(){
		Board newBoard = new Board(5, 6);
		Board oldBoard = new Board(5, 6);

		newBoard.place(pieces[Piece.SQUARE], 0, 0);
		newBoard.undo();
		checkBoardSimilarity(newBoard, oldBoard);

		newBoard.place(pieces[Piece.PYRAMID], 3, 1);
		newBoard.commit();
		oldBoard.place(pieces[Piece.PYRAMID], 3, 1);
		oldBoard.commit();
		checkBoardSimilarity(newBoard, oldBoard);
		newBoard.undo();
		checkBoardSimilarity(newBoard, oldBoard);

		newBoard.place(pieces[Piece.S2], 2, 3);
		newBoard.undo();
		checkBoardSimilarity(newBoard, oldBoard);
	}

	private void checkBoardSimilarity(Board b1, Board b2) {
		assertEquals(b1.getHeight(), b2.getHeight());
		assertEquals(b1.getWidth(), b2.getWidth());
		for(int i = 0; i < b1.getWidth(); i++){
			assertEquals(b1.getColumnHeight(i), b2.getColumnHeight(i));
			for (int j = 0; j < b1.getHeight(); j++) {
				assertEquals(b1.getRowWidth(j), b2.getRowWidth(j));
				assertEquals(b1.getGrid(i, j), b2.getGrid(i, j));
			}
		}
		assertEquals(b1.getMaxHeight(), b2.getMaxHeight());
	}

	@Test
	public void testUndoClearRowsEmpty(){
		Board newBoard = new Board(5, 6);
		Board oldBoard = new Board(5, 6);

		newBoard.clearRows();
		checkBoardSimilarity(newBoard, oldBoard);

		newBoard.undo();
		checkBoardSimilarity(newBoard, oldBoard);

	}

	@Test
	public void testUndoClearRowsFull(){
		Board testBoard = new Board(3, 4);
		Board fullBoard = new Board(3, 4);
		Board emptyBoard = new Board(3, 4);

		getFilledBoard3X4(fullBoard);
		getFilledBoard3X4(testBoard);

		testBoard.clearRows();
		checkBoardSimilarity(emptyBoard, testBoard);
		testBoard.undo();
		checkBoardSimilarity(fullBoard, testBoard);
	}

	@Test
	public void testUndoClearRowsNormal(){
		Board testBoard = new Board(5, 6);
		Board backupBoard = new Board(5, 6);

		assertEquals(Board.PLACE_OK, testBoard.place(pieces[Piece.STICK].fastRotation(), 0, 1));
		testBoard.commit();
		testBoard.place(pieces[Piece.STICK].fastRotation(), 0, 2);
		testBoard.commit();
		testBoard.place(pieces[Piece.STICK].fastRotation(), 0, 3);
		testBoard.commit();
		assertEquals(Board.PLACE_ROW_FILLED, testBoard.place(pieces[Piece.STICK], 4, 1));
		testBoard.commit();

		backupBoard.place(pieces[Piece.STICK].fastRotation(), 0, 1);
		backupBoard.commit();
		backupBoard.place(pieces[Piece.STICK].fastRotation(), 0, 2);
		backupBoard.commit();
		backupBoard.place(pieces[Piece.STICK].fastRotation(), 0, 3);
		backupBoard.commit();
		backupBoard.place(pieces[Piece.STICK], 4, 1);
		backupBoard.commit();


		assertEquals(3, testBoard.clearRows());
		testBoard.undo();
		checkBoardSimilarity(backupBoard, testBoard);

		assertEquals(3, testBoard.clearRows());
		testBoard.commit();
		assertEquals(3, backupBoard.clearRows());
		backupBoard.commit();

		testBoard.place(pieces[Piece.STICK].fastRotation(), 0, 1);
		testBoard.commit();
		backupBoard.place(pieces[Piece.STICK].fastRotation(), 0, 1);
		backupBoard.commit();

		assertEquals(1, testBoard.clearRows());
		testBoard.undo();
		checkBoardSimilarity(backupBoard, testBoard);
	}


	@Test
	public void testUndoBoth(){
		Board testBoard = new Board(4, 8);
		Board backupBoard = new Board(4, 8);

		testBoard.place(pieces[Piece.PYRAMID], 0, 0);
		testBoard.commit();
		backupBoard.place(pieces[Piece.PYRAMID], 0, 0);
		backupBoard.commit();

		assertEquals(Board.PLACE_ROW_FILLED, testBoard.place(pieces[Piece.STICK], 3, 0));
		assertEquals(1, testBoard.clearRows());

		testBoard.undo();
		checkBoardSimilarity(testBoard, backupBoard);
	}

	//
	// sanityCheck
	//

	@Test
	public void testSanityCheck(){
		Board b = new Board(3, 4);
		b.place(pieces[Piece.PYRAMID], 0, 0);
		b.commit();

		Piece fourthPyramidRot = pieces[Piece.PYRAMID].fastRotation().
				fastRotation().fastRotation();
		b.place(fourthPyramidRot, 0, 1);
		b.commit();

		int[] heights = new int[]{4, 3, 1};
		int[] widths = new int[]{1, 2, 2, 3};
		int maxHeigth = 4;

		assertDoesNotThrow(()-> b.checkHeights(heights));
		assertDoesNotThrow(()-> b.checkWidths(widths));
		assertDoesNotThrow(()-> b.checkMaxHeight(maxHeigth));

		int[] wrongHeights = new int[]{4, 3, 2};
		int[] wrongWidths = new int[]{1, 2, 0, 3};
		int wrongMaxHeight = 3;

		assertThrows(RuntimeException.class, ()-> b.checkHeights(wrongHeights));
		assertThrows(RuntimeException.class, ()-> b.checkWidths(wrongWidths));
		assertThrows(RuntimeException.class, ()-> b.checkMaxHeight(wrongMaxHeight));

		b.setDebugMode(false);
		Piece thirdL1Rot = pieces[Piece.L1].fastRotation().fastRotation();
		b.place(thirdL1Rot, 1, 1);
		b.commit();
		b.setDebugMode(true);
	}
}
