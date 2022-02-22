import java.util.*;

/*
 * Encapsulates a Sudoku grid to be solved.
 * CS108 Stanford.
 */
public class Sudoku {

    private int[][] grid;
    private String firstSolution;
    private long timeElapsed;
    // Provided grid data for main/testing
    // The instance variable strategy is up to you.

    // Provided easy 1 6 grid
    // (can paste this text into the GUI too)
    public static final int[][] easyGrid = Sudoku.stringsToGrid(
            "1 6 4 0 0 0 0 0 2",
            "2 0 0 4 0 3 9 1 0",
            "0 0 5 0 8 0 4 0 7",
            "0 9 0 0 0 6 5 0 0",
            "5 0 0 1 0 2 0 0 8",
            "0 0 8 9 0 0 0 3 0",
            "8 0 9 0 4 0 2 0 0",
            "0 7 3 5 0 9 0 0 1",
            "4 0 0 0 0 0 6 7 9");

    public static final int[][] emptyGrid = Sudoku.stringsToGrid(
            "0 0 0 0 0 0 0 0 0",
            "0 0 0 0 0 0 0 0 0",
            "0 0 0 0 0 0 0 0 0",
            "0 0 0 0 0 0 0 0 0",
            "0 0 0 0 0 0 0 0 0",
            "0 0 0 0 0 0 0 0 0",
            "0 0 0 0 0 0 0 0 0",
            "0 0 0 0 0 0 0 0 0",
            "0 0 0 0 0 0 0 0 0");


    // Provided medium 5 3 grid
    public static final int[][] mediumGrid = Sudoku.stringsToGrid(
            "530070000",
            "600195000",
            "098000060",
            "800060003",
            "400803001",
            "700020006",
            "060000280",
            "000419005",
            "000080079");

    // Provided hard 3 7 grid
    // 1 solution this way, 6 solutions if the 7 is changed to 0
    public static final int[][] hardGrid = Sudoku.stringsToGrid(
            "3 7 0 0 0 0 0 8 0",
            "0 0 1 0 9 3 0 0 0",
            "0 4 0 7 8 0 0 0 3",
            "0 9 3 8 0 0 0 1 2",
            "0 0 0 0 4 0 0 0 0",
            "5 2 0 0 0 6 7 9 0",
            "6 0 0 0 2 1 0 4 0",
            "0 0 0 5 3 0 9 0 0",
            "0 3 0 0 0 0 0 5 1");


    public static final int SIZE = 9;  // size of the whole 9x9 puzzle
    public static final int PART = 3;  // size of each 3x3 part
    public static final int MAX_SOLUTIONS = 100;
    public static final int EMPTY_SPOT_VALUE = 0;

    // Provided various static utility methods to
    // convert data formats to int[][] grid.

    /**
     * Returns a 2-d grid parsed from strings, one string per row.
     * The "..." is a Java 5 feature that essentially
     * makes "rows" a String[] array.
     * (provided utility)
     *
     * @param rows array of row strings
     * @return grid
     */
    public static int[][] stringsToGrid(String... rows) {
        int[][] result = new int[rows.length][];
        for (int row = 0; row < rows.length; row++) {
            result[row] = stringToInts(rows[row]);
        }
        return result;
    }


    /**
     * Given a single string containing 81 numbers, returns a 9x9 grid.
     * Skips all the non-numbers in the text.
     * (provided utility)
     *
     * @param text string of 81 numbers
     * @return grid
     */
    public static int[][] textToGrid(String text) {
        int[] nums = stringToInts(text);
        if (nums.length != SIZE * SIZE) {
            throw new RuntimeException("Needed 81 numbers, but got:" + nums.length);
        }

        int[][] result = new int[SIZE][SIZE];
        int count = 0;
        for (int row = 0; row < SIZE; row++) {
            for (int col = 0; col < SIZE; col++) {
                result[row][col] = nums[count];
                count++;
            }
        }
        return result;
    }


    /**
     * Given a string containing digits, like "1 23 4",
     * returns an int[] of those digits {1 2 3 4}.
     * (provided utility)
     *
     * @param string string containing ints
     * @return array of ints
     */
    public static int[] stringToInts(String string) {
        int[] a = new int[string.length()];
        int found = 0;
        for (int i = 0; i < string.length(); i++) {
            if (Character.isDigit(string.charAt(i))) {
                a[found] = Integer.parseInt(string.substring(i, i + 1));
                found++;
            }
        }
        int[] result = new int[found];
        System.arraycopy(a, 0, result, 0, found);
        return result;
    }


    // Provided -- the deliverable main().
    // You can edit to do easier cases, but turn in
    // solving hardGrid.
    public static void main(String[] args) {
        Sudoku sudoku;
        sudoku = new Sudoku(hardGrid);
        System.out.println(sudoku); // print the raw problem
        int count = sudoku.solve();
        System.out.println("solutions:" + count);
        System.out.println("elapsed:" + sudoku.getElapsed() + "ms");
        System.out.println(sudoku.getSolutionText());
    }


    /**
     * Sets up based on the given ints.
     * doesn't make a copy
     */
    public Sudoku(int[][] ints) {
        grid = ints;
        firstSolution = "";
    }

    public Sudoku(String text) {
        this(textToGrid(text));
    }

    /**
     * Solves the puzzle, invoking the underlying recursive search.
     */
    public int solve() {
        long initialTime = System.currentTimeMillis();
        List<Spot> emptySpots = new ArrayList<>();
        getEmptySpots(emptySpots);
        boolean[] mutableBoolean = {true};
        int solutionCount = getSolutionCount(emptySpots, 0, mutableBoolean);
        timeElapsed = System.currentTimeMillis() - initialTime;
        return solutionCount;
    }

    private void getEmptySpots(List<Spot> emptySpots) {
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                if (grid[i][j] == EMPTY_SPOT_VALUE) {
                    Spot emptySpot = new Spot(i, j);
                    emptySpot.getPossibleNumbers();
                    emptySpots.add(emptySpot);
                }
            }
        }
        Collections.sort(emptySpots);
    }

    private int getSolutionCount(List<Spot> emptySpots, int currIndex, boolean[] firstSolution) {
        if (currIndex == emptySpots.size()) {
            if (firstSolution[0]) {
                firstSolution[0] = false;
                this.firstSolution = toString();
            }
            return 1;
        }

        Spot spot = emptySpots.get(currIndex);
        ArrayList<Integer> possibleNums = spot.getPossibleNumbers();

        int solutionCount = 0;
        for (int possibleNum : possibleNums) {
            spot.setValue(possibleNum);
            solutionCount += getSolutionCount(emptySpots, currIndex + 1, firstSolution);
            spot.setValue(EMPTY_SPOT_VALUE);
            if (solutionCount >= MAX_SOLUTIONS) return MAX_SOLUTIONS;
        }
        return solutionCount;
    }

    public String getSolutionText() {
        return firstSolution;
    }

    public long getElapsed() {
        return timeElapsed;
    }

    @Override
    public String toString() {
        StringBuilder board = new StringBuilder();
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                board.append(grid[i][j]);
                if (j != SIZE - 1) board.append(" ");
            }
            board.append("\n");
        }
        return board.toString();
    }


    public class Spot implements Comparable<Spot> {
        private int row;
        private int col;
        private ArrayList<Integer> possibleNums;

        public Spot(int row, int col) {
            this.row = row;
            this.col = col;
            possibleNums = new ArrayList<>();
        }

        public int getValue() {
            return grid[row][col];
        }

        public void setValue(int value) {
            grid[row][col] = value;
        }

        //returns pointer to internal arrayList, doesn't make a copy
        public ArrayList<Integer> getPossibleNumbers() {
            possibleNums.clear();
            int[] arr = new int[10]; //all integers from 0 to 9 used for determining what is left
            getUsedNumbers(arr);
            for (int i = 1; i < arr.length; i++) {
                if (arr[i] == 0) possibleNums.add(i);
            }
            return possibleNums;
        }

        private void getUsedNumbers(int[] arr) {
            for (int i = 0; i < SIZE; i++) {
                arr[grid[row][i]]++;
                arr[grid[i][col]]++;
            }

            //for quadratic iteration
            for (int i = row - row % PART; i < row - row % PART + PART; i++) {
                for (int j = col - col % PART; j < col - col % PART + PART; j++) {
                    arr[grid[i][j]]++;
                }
            }
        }

        @Override
        public int compareTo(Spot other) {
            return Integer.compare(this.possibleNums.size(), other.possibleNums.size());
        }
    }

}
