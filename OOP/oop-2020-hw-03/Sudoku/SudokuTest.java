import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import static org.junit.jupiter.api.Assertions.*;

class SudokuTest {
    @Test
    void sudoku() {
        Sudoku s1 = new Sudoku(Sudoku.easyGrid);
        Sudoku s2 = new Sudoku("1 6 4 0 0 0 0 0 2\n"
                + "2 0 0 4 0 3 9 1 0\n"
                + "0 0 5 0 8 0 4 0 7\n"
                + "0 9 0 0 0 6 5 0 0\n"
                + "5 0 0 1 0 2 0 0 8\n"
                + "0 0 8 9 0 0 0 3 0\n"
                + "8 0 9 0 4 0 2 0 0\n"
                + "0 7 3 5 0 9 0 0 1\n"
                + "4 0 0 0 0 0 6 7 9");

        assertThrows(RuntimeException.class, () -> {
            Sudoku s3 = new Sudoku("1 6 4 0 0 0 0 0 2\n"
                    + "2 0 0 4 0 3 9 1 0\n"
                    + "0 0 5 0 8 0 4 0 7\n"
                    + "0 9 0 0 0 6 5 0 0\n"
                    + "5 0 0 1 0 2 0 0 8\n"
                    + "0 0 8 9 0 0 0 3 0\n"
                    + "8 0 9 0 4 0 2 0 0\n"
                    + "0 7 3 5 0 9 0 0 1\n"
                    + "4 0 0 0 0 0 6 7 9 234r54567576");
        });
    }

    @Test
    void main() {
        String[] arr = new String[1];
        Sudoku.main(arr);
    }

    @Test
    void solve() {
        Sudoku s = new Sudoku(Sudoku.hardGrid);
        assertEquals(1, s.solve());
        Sudoku changedHardGrid = new Sudoku(
                "3 0 0 0 0 0 0 8 0\n"
                        + "0 0 1 0 9 3 0 0 0\n"
                        + "0 4 0 7 8 0 0 0 3\n"
                        + "0 9 3 8 0 0 0 1 2\n"
                        + "0 0 0 0 4 0 0 0 0\n"
                        + "5 2 0 0 0 6 7 9 0\n"
                        + "6 0 0 0 2 1 0 4 0\n"
                        + "0 0 0 5 3 0 9 0 0\n"
                        + "0 3 0 0 0 0 0 5 1");
        assertEquals(6, changedHardGrid.solve());
        Sudoku empty = new Sudoku(Sudoku.emptyGrid);
        assertEquals(100, empty.solve());

        Sudoku unsolvable = new Sudoku(
                "3 0 0 0 0 0 0 8 0\n"
                        + "0 0 1 0 9 3 0 0 0\n"
                        + "0 4 0 7 8 0 0 0 3\n"
                        + "0 9 3 8 0 0 0 1 2\n"
                        + "0 0 0 0 4 0 0 0 0\n"
                        + "5 2 0 0 0 6 7 9 0\n"
                        + "6 0 0 0 2 1 0 4 7\n"
                        + "0 0 0 5 3 0 9 6 1\n"
                        + "0 3 0 0 0 0 8 5 0");

        assertEquals(0, unsolvable.solve());

        Sudoku alreadySolved = new Sudoku("3 7 5 1 6 2 4 8 9\n"
                + "8 6 1 4 9 3 5 2 7\n"
                + "2 4 9 7 8 5 1 6 3\n"
                + "4 9 3 8 5 7 6 1 2\n"
                + "7 1 6 2 4 9 8 3 5\n"
                + "5 2 8 3 1 6 7 9 4\n"
                + "6 5 7 9 2 1 3 4 8\n"
                + "1 8 2 5 3 4 9 7 6\n"
                + "9 3 0 6 7 8 2 5 1");
        assertEquals(1, alreadySolved.solve());
    }

    //
    // Spot
    //

    //constructor, setter, getter
    @Test
    public void spot() {
        Sudoku sud = new Sudoku(Sudoku.easyGrid);
        Sudoku.Spot s1 = sud.new Spot(0, 0);
        Sudoku.Spot s2 = sud.new Spot(0, 0);
        assertEquals(s1.getValue(), s2.getValue());
        s1.setValue(5);
        assertEquals(s1.getValue(), 5);
        assertEquals(s1.getValue(), s2.getValue());
        assertEquals(0, s1.compareTo(s2));
    }

    @Test
    public void getPossibleNumbers() {
        Sudoku sud = new Sudoku(Sudoku.hardGrid);
        Sudoku.Spot s1 = sud.new Spot(0, 2);
        assertEquals(Sudoku.EMPTY_SPOT_VALUE, s1.getValue());
        ArrayList<Integer> expected = new ArrayList<>();
        Collections.addAll(expected, 2, 5, 6, 9);
        assertEquals(expected, s1.getPossibleNumbers());
    }
}