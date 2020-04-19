import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

import static org.junit.jupiter.api.Assertions.*;

/*
  Unit test for Piece class -- starter shell.
 */
public class PieceTest{

	// Contains all possible 19 pieces of Tetris in this order:
	// 0 - 1 all possible rotations of "The Stick"
	// 2 - 5 all possible rotations of "The L1"
	// 6 - 9 all possible rotations of "The L2"
	// 10 - 11 all possible rotations of "The S1"
	// 12 - 13 all possible rotations of "The S2"
	// 14 "The Square"
	// 15 - 18 all possible rotations of "The Pyramid"
	private static Piece[] pieces;


	@BeforeAll
	private static void init(){
		pieces = new Piece[19];
		pieces[STICK_0] = new Piece("0 0  0 1  0 2  0 3");
		pieces[STICK_1] = new Piece("0 0  1 0  2 0  3 0");
		pieces[L1_0] = new Piece("0 0  0 1  0 2  1 0");
		pieces[L1_1] = new Piece("0 0  1 0  2 0  2 1");
		pieces[L1_2] = new Piece("1 0  1 1  1 2  0 2");
		pieces[L1_3] = new Piece("0 0  0 1  1 1  2 1");
		pieces[L2_0] = new Piece("0 0  1 0  1 1  1 2");
		pieces[L2_1] = new Piece("2 1  1 1  0 1  2 0");
		pieces[L2_2] = new Piece("0 0  0 1  0 2  1 2");
		pieces[L2_3] = new Piece("2 0  1 0  0 0  0 1");
		pieces[S1_0] = new Piece("0 0  1 0  1 1  2 1");
		pieces[S1_1] = new Piece("0 1  0 2  1 0  1 1");
		pieces[S2_0] = new Piece("0 1  1 1  1 0  2 0");
		pieces[S2_1] = new Piece("0 0  0 1  1 1  1 2");
		pieces[SQUARE_0] = new Piece("0 0  0 1  1 0  1 1");
		pieces[P_0] = new Piece("0 0  1 0  1 1  2 0");
		pieces[P_1] = new Piece("0 1  1 0  1 1  1 2");
		pieces[P_2] = new Piece("0 1  1 0  1 1  2 1");
		pieces[P_3] = new Piece("0 0  0 1  1 1  0 2");
	}


	//
	// Constructor, getters
	//

	@Test
	public void checkConstructorImmutability(){
		TPoint[] points = new TPoint[4];
		points[0] = new TPoint(0,0);
		points[1] = new TPoint(1,0);
		points[2] = new TPoint(1,1);
		points[3] = new TPoint(3,2);

		TPoint[] pointsCopy = points.clone();
		Piece p = new Piece(points);
		points[2] = new TPoint(1225, 223);

		assertArrayEquals(pointsCopy, p.getBody());
	}

	@Test
	public void testWrongConstructor(){
		assertThrows(RuntimeException.class, () -> {
			Piece WrongPiece = new Piece("a b 1 0  1 1  2 1");
		});

		Piece p = new Piece("");
	}

	//little bit too much, but at least its exhaustive
	@Test
	public void testHeight(){
		assertEquals(4, pieces[STICK_0].getHeight());
		assertEquals(1, pieces[STICK_1].getHeight());
		assertEquals(3, pieces[L1_0].getHeight());
		assertEquals(2, pieces[L1_1].getHeight());
		assertEquals(3, pieces[L1_2].getHeight());
		assertEquals(2, pieces[L1_3].getHeight());
		assertEquals(3, pieces[L2_0].getHeight());
		assertEquals(2, pieces[L2_1].getHeight());
		assertEquals(3, pieces[L2_2].getHeight());
		assertEquals(2, pieces[L2_3].getHeight());
		assertEquals(2, pieces[S1_0].getHeight());
		assertEquals(3, pieces[S1_1].getHeight());
		assertEquals(2, pieces[S2_0].getHeight());
		assertEquals(3, pieces[S2_1].getHeight());
		assertEquals(2, pieces[SQUARE_0].getHeight());
		assertEquals(2, pieces[P_0].getHeight());
		assertEquals(3, pieces[P_1].getHeight());
		assertEquals(2, pieces[P_2].getHeight());
		assertEquals(3, pieces[P_3].getHeight());
	}

	@Test
	public void testWidth(){
		assertEquals(1, pieces[STICK_0].getWidth());
		assertEquals(4, pieces[STICK_1].getWidth());
		assertEquals(2, pieces[L1_0].getWidth());
		assertEquals(3, pieces[L1_1].getWidth());
		assertEquals(2, pieces[L1_2].getWidth());
		assertEquals(3, pieces[L1_3].getWidth());
		assertEquals(2, pieces[L2_0].getWidth());
		assertEquals(3, pieces[L2_1].getWidth());
		assertEquals(2, pieces[L2_2].getWidth());
		assertEquals(3, pieces[L2_3].getWidth());
		assertEquals(3, pieces[S1_0].getWidth());
		assertEquals(2, pieces[S1_1].getWidth());
		assertEquals(3, pieces[S2_0].getWidth());
		assertEquals(2, pieces[S2_1].getWidth());
		assertEquals(2, pieces[SQUARE_0].getWidth());
		assertEquals(3, pieces[P_0].getWidth());
		assertEquals(2, pieces[P_1].getWidth());
		assertEquals(3, pieces[P_2].getWidth());
		assertEquals(2, pieces[P_3].getWidth());
	}



	@Test
	public void testSkirt(){
		assertArrayEquals(new int[]{0}, pieces[STICK_0].getSkirt());
		assertArrayEquals(new int[]{0, 0, 0, 0}, pieces[STICK_1].getSkirt());
		assertArrayEquals(new int[]{0, 0}, pieces[L1_0].getSkirt());
		assertArrayEquals(new int[]{0, 0, 0}, pieces[L1_1].getSkirt());
		assertArrayEquals(new int[]{2, 0}, pieces[L1_2].getSkirt());
		assertArrayEquals(new int[]{0, 1, 1}, pieces[L1_3].getSkirt());
		assertArrayEquals(new int[]{0, 0}, pieces[L2_0].getSkirt());
		assertArrayEquals(new int[]{1, 1, 0}, pieces[L2_1].getSkirt());
		assertArrayEquals(new int[]{0, 2}, pieces[L2_2].getSkirt());
		assertArrayEquals(new int[]{0, 0, 0}, pieces[L2_3].getSkirt());
		assertArrayEquals(new int[]{0, 0, 1}, pieces[S1_0].getSkirt());
		assertArrayEquals(new int[]{1, 0}, pieces[S1_1].getSkirt());
		assertArrayEquals(new int[]{1, 0, 0}, pieces[S2_0].getSkirt());
		assertArrayEquals(new int[]{0, 1}, pieces[S2_1].getSkirt());
		assertArrayEquals(new int[]{0, 0}, pieces[SQUARE_0].getSkirt());
		assertArrayEquals(new int[]{0, 0, 0}, pieces[P_0].getSkirt());
		assertArrayEquals(new int[]{1, 0}, pieces[P_1].getSkirt());
		assertArrayEquals(new int[]{1, 0, 1}, pieces[P_2].getSkirt());
		assertArrayEquals(new int[]{0, 1}, pieces[P_3].getSkirt());
	}


	//
	// equals
	//

	@Test
	public void testEqualityNull(){
		TPoint[] points = new TPoint[4];
		points[0] = new TPoint(0,0);
		points[1] = new TPoint(1,0);
		points[2] = new TPoint(1,1);
		points[3] = new TPoint(3,2);
		Piece p = new Piece(points);

		assertNotEquals(null, p);
		assertNotEquals(null, pieces[P_3]);
		assertNotEquals(new Piece(""), p);
	}

	@Test
	public void testEqualityItself(){
		assertEquals(pieces[STICK_0], pieces[STICK_0]);
		assertEquals(pieces[STICK_1], pieces[STICK_1]);
		assertEquals(pieces[L1_0], pieces[L1_0]);
		assertEquals(pieces[L1_1], pieces[L1_1]);
		assertEquals(pieces[L1_2], pieces[L1_2]);
		assertEquals(pieces[L1_3], pieces[L1_3]);
		assertEquals(pieces[L2_0], pieces[L2_0]);
		assertEquals(pieces[L2_1], pieces[L2_1]);
		assertEquals(pieces[L2_2], pieces[L2_2]);
		assertEquals(pieces[L2_3], pieces[L2_3]);
		assertEquals(pieces[S1_0], pieces[S1_0]);
		assertEquals(pieces[S1_1], pieces[S1_1]);
		assertEquals(pieces[S2_0], pieces[S2_0]);
		assertEquals(pieces[S2_1], pieces[S2_1]);
		assertEquals(pieces[SQUARE_0], pieces[SQUARE_0]);
		assertEquals(pieces[P_0], pieces[P_0]);
		assertEquals(pieces[P_1], pieces[P_1]);
		assertEquals(pieces[P_2], pieces[P_2]);
		assertEquals(pieces[P_3], pieces[P_3]);
	}

	@Test
	public void testEqualityNormal(){
		TPoint[] points = new TPoint[4];
		points[0] = new TPoint(0,0);
		points[1] = new TPoint(1,0);
		points[2] = new TPoint(1,1);
		points[3] = new TPoint(3,2);
		Piece p = new Piece(points);
		TPoint[] pointsCopy = points.clone();
		Piece pCopy = new Piece(pointsCopy);

		assertEquals(p, pCopy);
		assertNotEquals(p, pieces[STICK_1]);
		assertNotEquals(pCopy, pieces[L2_1]);
		List<String> ls = new ArrayList<>();
		assertNotEquals(p, ls);

		points[3] = new TPoint(2, 0);
		p = new Piece(points);
		assertEquals(p, pieces[P_0]);
	}

	//
	// computeNextRotation
	//

	@Test
	public void testComputeNextRotation(){
		assertEquals(pieces[STICK_1], pieces[STICK_0].computeNextRotation());
		assertEquals(pieces[STICK_0], pieces[STICK_1].computeNextRotation());
		assertEquals(pieces[L1_1], pieces[L1_0].computeNextRotation());
		assertEquals(pieces[L1_2], pieces[L1_1].computeNextRotation());
		assertEquals(pieces[L1_3], pieces[L1_2].computeNextRotation());
		assertEquals(pieces[L1_0], pieces[L1_3].computeNextRotation());
		assertEquals(pieces[L2_1], pieces[L2_0].computeNextRotation());
		assertEquals(pieces[L2_2], pieces[L2_1].computeNextRotation());
		assertEquals(pieces[L2_3], pieces[L2_2].computeNextRotation());
		assertEquals(pieces[L2_0], pieces[L2_3].computeNextRotation());
		assertEquals(pieces[S1_1], pieces[S1_0].computeNextRotation());
		assertEquals(pieces[S1_0], pieces[S1_1].computeNextRotation());
		assertEquals(pieces[S2_1], pieces[S2_0].computeNextRotation());
		assertEquals(pieces[S2_0], pieces[S2_1].computeNextRotation());
		assertEquals(pieces[SQUARE_0], pieces[SQUARE_0].computeNextRotation());
		assertEquals(pieces[P_1], pieces[P_0].computeNextRotation());
		assertEquals(pieces[P_2], pieces[P_1].computeNextRotation());
		assertEquals(pieces[P_3], pieces[P_2].computeNextRotation());
		assertEquals(pieces[P_0], pieces[P_3].computeNextRotation());
	}

	//
	// getPieces
	//

	@Test
	public void testGetPiecesOriginal(){
		Piece[] figures = Piece.getPieces();
		assertEquals(pieces[STICK_0], figures[Piece.STICK]);
		assertEquals(pieces[L1_0], figures[Piece.L1]);
		assertEquals(pieces[L2_0], figures[Piece.L2]);
		assertEquals(pieces[S1_0], figures[Piece.S1]);
		assertEquals(pieces[S2_0], figures[Piece.S2]);
		assertEquals(pieces[SQUARE_0], figures[Piece.SQUARE]);
		assertEquals(pieces[P_0], figures[Piece.PYRAMID]);
	}


	//
	// fastRotation
	//

	@Test
	public void testFastRotations(){
		Piece[] figures = Piece.getPieces();
		testRotation(figures[Piece.STICK], STICK_0);
		testRotation(figures[Piece.L1], L1_0);
		testRotation(figures[Piece.L2], L2_0);
		testRotation(figures[Piece.S1], S1_0);
		testRotation(figures[Piece.S2], S2_0);
		testRotation(figures[Piece.SQUARE], SQUARE_0);
		testRotation(figures[Piece.PYRAMID], P_0);

	}

	private void testRotation(Piece root, int index){
		Piece lastPiece = root;
		int rotationNum = getRotationNum(index);
		for(int i = 0; i < 1000; i++){
			assertEquals(pieces[index + i%rotationNum], lastPiece);
			lastPiece = lastPiece.fastRotation();
		}
	}

	private int getRotationNum(int index) {
		if(index == STICK_0 || index == S1_0 || index == S2_0) return 2;
		if (index == SQUARE_0) return 1;
		return 4;
	}


	//enums don't work in java as in c/c++ ...
	private static final int STICK_0 = 0;
	private static final int STICK_1 = 1;
	private static final int L1_0 = 2;
	private static final int L1_1 = 3;
	private static final int L1_2 = 4;
	private static final int L1_3 = 5;
	private static final int L2_0 = 6;
	private static final int L2_1 = 7;
	private static final int L2_2 = 8;
	private static final int L2_3 = 9;
	private static final int S1_0 = 10;
	private static final int S1_1 = 11;
	private static final int S2_0 = 12;
	private static final int S2_1 = 13;
	private static final int SQUARE_0 = 14;
	private static final int P_0= 15;
	private static final int P_1 = 16;
	private static final int P_2 = 17;
	private static final int P_3 = 18;
}