// Board.java

/**
 CS108 Tetris Board.
 Represents a Tetris board -- essentially a 2-d grid
 of booleans. Supports tetris pieces and row clearing.
 Has an "undo" feature that allows clients to add and remove pieces efficiently.
 Does not do any drawing or have any idea of pixels. Instead,
 just represents the abstract 2-d board.
*/
public class Board	{
	// Some ivars are stubbed out for you:
	private int width;
	private int height;
	private boolean[][] grid;
	private boolean DEBUG = true;
	boolean committed;
	private int[] widths;
	private int[] heights;
	private int maxHeight;

	private boolean[][] gridBackup;
	private int[] widthsBackup;
	private int[] heightsBackup;
	private int maxHeightBackup;
	private boolean dataBackedUp;


	// Here a few trivial methods are provided:
	
	/**
	 Creates an empty board of the given width and height
	 measured in blocks.
	*/
	public Board(int width, int height) {
		this.width = width;
		this.height = height;
		grid = new boolean[width][height];
		widths = new int[height];
		heights = new int[width];
		maxHeight = 0;
		widthsBackup = new int[height];
		heightsBackup = new int[width];
		maxHeightBackup = 0;
		gridBackup = new boolean[width][height];
		dataBackedUp = false;
		committed = true;
	}
	
	
	/**
	 Returns the width of the board in blocks.
	*/
	public int getWidth() {
		return width;
	}
	
	
	/**
	 Returns the height of the board in blocks.
	*/
	public int getHeight() {
		return height;
	}
	
	
	/**
	 Returns the max column height present in the board.
	 For an empty board this is 0.
	*/
	public int getMaxHeight() {
		return maxHeight;
	}
	
	
	/**
	 Checks the board for internal consistency -- used
	 for debugging.
	*/
	public void sanityCheck() {
		if (!DEBUG) return;
		checkHeights(heights);
		checkWidths(widths);
		checkMaxHeight(maxHeight);
	}

	protected void checkHeights(int[] heights) {
		for(int i = 0; i < width; i++){
			int colHeight = 0;
			for (int j = 0; j < height; j++) {
				if(grid[i][j]){
					colHeight = height - j;
					break;
				}
			}
			if(colHeight != heights[i])
				throw new RuntimeException
						("Column height is wrong, expected " + colHeight + " got " + heights[i]);
		}
	}

	protected void checkWidths(int[] widths) {
		for(int i = 0; i < height; i++){
			int rowWidth = 0;
			for (int j = 0; j < width; j++) {
				if(grid[j][i]) rowWidth++;
			}
			if(rowWidth != widths[i])
				throw new RuntimeException("Row width is wrong, expected " + rowWidth + " got " + widths[i]);
		}
	}

	protected void checkMaxHeight(int maxHeight) {
		int maxH = 0;
		for (int height : heights)
			maxH = Math.max(maxH, height);

		if(maxH != maxHeight)
			throw new RuntimeException("Max height is wrong, expected " + maxH + " got " + maxHeight);
	}

	protected void setDebugMode(boolean mode){
		DEBUG = mode;
	}
	/**
	 Given a piece and an x, returns the y
	 value where the piece would come to rest
	 if it were dropped straight down at that x.
	 
	 <p>
	 Implementation: use the skirt and the col heights
	 to compute this fast -- O(skirt length).
	*/
	public int dropHeight(Piece piece, int x) {
		int[] skirt = piece.getSkirt();
		int initialY = height;
		int minDifference = height;

		for (int i = 0; i < skirt.length; i++)
			minDifference = Math.min(minDifference,
					initialY + skirt[i] - getColumnHeight(x + i));

		return initialY - minDifference;
	}

	/**
	 Returns the height of the given column --
	 i.e. the y value of the highest block + 1.
	 The height is 0 if the column contains no blocks.
	*/
	public int getColumnHeight(int x) {
		return heights[x];
	}
	
	
	/**
	 Returns the number of filled blocks in
	 the given row.
	*/
	public int getRowWidth(int y) {
		 return widths[height - 1 - y];
	}
	
	
	/**
	 Returns true if the given block is filled in the board.
	 Blocks outside of the valid width/height area
	 always return true.
	*/
	public boolean getGrid(int x, int y) {
		return checkCell(x, height - 1 - y) != PLACE_OK;
	}
	
	
	public static final int PLACE_OK = 0;
	public static final int PLACE_ROW_FILLED = 1;
	public static final int PLACE_OUT_BOUNDS = 2;
	public static final int PLACE_BAD = 3;
	
	/**
	 Attempts to add the body of a piece to the board.
	 Copies the piece blocks into the board grid.
	 Returns PLACE_OK for a regular placement, or PLACE_ROW_FILLED
	 for a regular placement that causes at least one row to be filled.
	 
	 <p>Error cases:
	 A placement may fail in two ways. First, if part of the piece may falls out
	 of bounds of the board, PLACE_OUT_BOUNDS is returned.
	 Or the placement may collide with existing blocks in the grid
	 in which case PLACE_BAD is returned.
	 In both error cases, the board may be left in an invalid
	 state. The client can use undo(), to recover the valid, pre-place state.
	*/
	public int place(Piece piece, int x, int y) {
		// flag !committed problem
		if (!committed) throw new RuntimeException("place commit problem");

		createBackup();
		boolean[] results = new boolean[PLACE_BAD + 1];

		TPoint[] arr = piece.getBody();
		for (TPoint tPoint : arr) {
			int result = fillCell(x + tPoint.x, height - 1 - y - tPoint.y);
			results[result] = true;
		}

		if(results[PLACE_OUT_BOUNDS]) return PLACE_OUT_BOUNDS;
		if(results[PLACE_BAD]) return PLACE_BAD;
		sanityCheck();
		if(results[PLACE_ROW_FILLED]) return PLACE_ROW_FILLED;
		return PLACE_OK;
	}

	private void createBackup() {
		committed = false;
		for (int i = 0; i < grid.length; i++)
			System.arraycopy(grid[i], 0, gridBackup[i], 0, grid[i].length);

		System.arraycopy(widths, 0, widthsBackup, 0, widths.length);
		System.arraycopy(heights, 0, heightsBackup, 0, heights.length);
		maxHeightBackup = maxHeight;
		dataBackedUp = true;
	}

	private int fillCell(int x, int y) {
		int result = checkCell(x, y);

		if(result == PLACE_OK) {
			grid[x][y] = true;
			widths[y]++;
			heights[x] = Math.max(heights[x], height - y);
			maxHeight = Math.max(heights[x], maxHeight);
			if(widths[y] == width) result = PLACE_ROW_FILLED;
		}

		return result;
	}

	private int checkCell(int x, int y) {
		if(x >= width || x < 0 || y >= height || y < 0) return PLACE_OUT_BOUNDS;

		if(grid[x][y])  return PLACE_BAD;

		return PLACE_OK;
	}


	/**
	 Deletes rows that are filled all the way across, moving
	 things above down. Returns the number of rows cleared.
	*/
	public int clearRows() {
		int rowsCleared = 0;
		committed = false;
		if(!dataBackedUp) createBackup();
		for(int i = height - 1; i > height - 1 - maxHeight; i--){
			slideDown(rowsCleared, i);
			if(widths[i] == width) rowsCleared++;
		}
		makeNewRows(rowsCleared);
		sanityCheck();
		return rowsCleared;
	}

	private void slideDown(int rowsCleared, int i) {
		for(int j = 0; j < width; j++) {
			grid[j][i + rowsCleared] = grid[j][i];
			widths[i + rowsCleared] = widths[i];
		}
	}

	private void makeNewRows(int rowsCleared) {
		for (int i = rowsCleared; i > 0; i--) {
			widths[height - 1 - maxHeight + i] = 0;
			for (int j = 0; j < width; j++)
				grid[j][height - 1 - maxHeight + i] = false;
		}
		updateHeights(rowsCleared);
	}

	private void updateHeights(int rowsCleared) {
		int maxH = 0;
		for(int i = 0; i < width; i++){
			int colHeight = 0;
			for (int j = height - heights[i] + rowsCleared; j < height; j++) {
				if(grid[i][j]){
					colHeight = height - j;
					break;
				}
			}
			heights[i] = colHeight;
			maxH = Math.max(maxH, colHeight);
		}
		maxHeight = maxH;
	}

	/**
	 Reverts the board to its state before up to one place
	 and one clearRows();
	 If the conditions for undo() are not met, such as
	 calling undo() twice in a row, then the second undo() does nothing.
	 See the overview docs.
	*/
	public void undo() {
		if(!dataBackedUp) return;

		boolean[][] tempGrid = grid;
		grid = gridBackup;
		gridBackup = tempGrid;

		int[] tempWidths = widths;
		widths = widthsBackup;
		widthsBackup = tempWidths;

		int[] tempHeights = heights;
		heights = heightsBackup;
		heightsBackup = tempHeights;

		int tempMaxHeight = maxHeight;
		maxHeight = maxHeightBackup;
		maxHeightBackup = tempMaxHeight;
		sanityCheck();
		commit();
	}
	
	
	/**
	 Puts the board in the committed state.
	*/
	public void commit() {
		committed = true;
		dataBackedUp = false;
	}


	
	/*
	 Renders the board state as a big String, suitable for printing.
	 This is the sort of print-obj-state utility that can help see complex
	 state change over time.
	 (provided debugging utility) 
	 */
	public String toString() {
		StringBuilder buff = new StringBuilder();
		for (int y = height-1; y>=0; y--) {
			buff.append('|');
			for (int x=0; x<width; x++) {
				if (getGrid(x,y)) buff.append('+');
				else buff.append(' ');
			}
			buff.append("|\n");
		}
		for (int x=0; x<width+2; x++) buff.append('-');
		return(buff.toString());
	}
}
