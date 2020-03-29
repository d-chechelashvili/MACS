//
// TetrisGrid encapsulates a tetris board and has
// a clearRows() capability.

import java.util.Arrays;

public class TetrisGrid {
	private boolean[][] grid;
	/**
	 * Constructs a new instance with the given grid.
	 * Does not make a copy.
	 * @param grid
	 */
	public TetrisGrid(boolean[][] grid) {
		this.grid = grid;
	}
	
	
	/**
	 * Does row-clearing on the grid (see handout).
	 */
	public void clearRows() {
		int rowsCleared = 0;
		for(int i = 0 ; i < this.grid[0].length; i++){
			slideDown(rowsCleared, i); //slideleft
			if(rowFull(i)) rowsCleared++;
		}
		makeNewRows(rowsCleared);
	}

	private void slideDown(int rowsCleared, int i) {
		for(int j = 0; j < this.grid.length; j++){
			this.grid[j][i - rowsCleared] = this.grid[j][i];
		}
	}

	private boolean rowFull(int i) {
		for(int j = 0; j < this.grid.length; j++){
			if(this.grid[j][i] == false) return false;
		}
		return true;
	}

	private void makeNewRows(int rowsCleared) {
		for(int i = 1 ; i <= rowsCleared; i++){
			for(int j = 0; j < this.grid.length; j++){
				this.grid[j][this.grid[0].length - i] = false;
			}
		}
	}

	/**
	 * Returns the internal 2d grid array.
	 * @return 2d grid array
	 */
	boolean[][] getGrid() {
		//makes a copy to prevent any error from passing pointer to user
		return this.grid.clone();
	}
}
