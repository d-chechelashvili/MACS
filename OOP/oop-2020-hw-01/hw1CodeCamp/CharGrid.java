// HW1 2-d array Problems
// CharGrid encapsulates a 2-d grid of chars and supports
// a few operations on the grid.

public class CharGrid {
	private char[][] grid;

	/**
	 * Constructs a new CharGrid with the given grid.
	 * Does not make a copy.
	 * @param grid
	 */
	public CharGrid(char[][] grid) {
		this.grid = grid;
	}
	
	/**
	 * Returns the area for the given char in the grid. (see handout).
	 * @param ch char to look for
	 * @return area for given char
	 */
	public int charArea(char ch) {
		//captures uttermost points in all 4 directions and calculates area of rectangle made out of them
		if(this.grid.length < 1) return 0;

		int rightmost = -this.grid[0].length, leftmost = 1 - rightmost,
				lowest = -this.grid.length, highest = 1 - lowest;

		for(int i = 0; i < this.grid.length; i++){
			for(int j = 0; j < this.grid[i].length; j++){
				if(this.grid[i][j] == ch){
					if(j < leftmost) leftmost = j;
					if(j > rightmost) rightmost = j;
					if(i < highest) highest = i;
					if(i > lowest) lowest = i;
				}
			}
		}

		return (Math.abs(rightmost) - leftmost + 1) * (Math.abs(lowest) - highest + 1);
	}
	
	/**
	 * Returns the count of '+' figures in the grid (see handout).
	 * @return number of + in grid
	 */
	public int countPlus() {
		int plusNum = 0;

		for(int i = 0; i < this.grid.length; i++){
			for(int j = 0; j < this.grid[i].length; j++){
				if(isCentre(i, j, this.grid[i][j])) plusNum++;
			}
		}

		return plusNum;
	}

	private boolean isCentre(int i, int j, char ch) {
		int topArm = getArmLength(i, j, -1, 0, ch);
		int bottomArm = getArmLength(i, j, 1, 0, ch);
		int leftArm = getArmLength(i, j, 0, -1, ch);
		int rightArm = getArmLength(i, j, 0, 1, ch);

		if(topArm != bottomArm || topArm != leftArm || topArm != rightArm || topArm < 2) return false;
		return true;
	}

	private int getArmLength(int i, int j, int deltaI, int deltaJ, char ch) {
		if(!inbounds(i, j) || this.grid[i][j] != ch) return 0;
		return getArmLength(i + deltaI, j + deltaJ, deltaI, deltaJ, ch) + 1;
	}

	private boolean inbounds(int i, int j) {
		if( i < 0 || j < 0 || i >= this.grid.length || j >= this.grid[i].length)
			return false;
		return true;
	}


}
