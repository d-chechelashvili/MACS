
/*
 * File: Pyramid.java
 * Name: 
 * Section Leader: 
 * ------------------
 * This file is the starter file for the Pyramid problem.
 * It includes definitions of the constants that match the
 * sample run in the assignment, but you should make sure
 * that changing these values causes the generated display
 * to change accordingly.
 */

import acm.graphics.*;
import acm.program.*;
import java.awt.*;

public class Pyramid extends GraphicsProgram {

	/** Width of each brick in pixels */
	private static final int BRICK_WIDTH = 30;

	/** Width of each brick in pixels */
	private static final int BRICK_HEIGHT = 12;

	/** Number of bricks in the base of the pyramid */
	private static final int BRICKS_IN_BASE = 14;

	/*
	 * this program must draw pyramid which is symmetrical to the vertical line at
	 * the middle of the screen. pyramid must have 14 tiers of rows each row has
	 * blocks equal to its number numbering starts from bottom(e.g. 14th tier has 14
	 * blocks, 7th tier has 7 blocks e.t.c)
	 */
	public void run() {
		// x0&y0 are coordinates of the starting point of leftmost block in the first
		// row
		int x0 = getWidth() / 2 - 7 * BRICK_WIDTH;
		int y0 = getHeight() - BRICK_HEIGHT;
		for (int i = BRICKS_IN_BASE; i > 0; i--) {
			drawTier(x0, y0, i);
			x0 += BRICK_WIDTH / 2;
			y0 -= BRICK_HEIGHT;
		}
	}

	private void drawTier(int x0, int y0, int i) {
		while (i > 0) {
			GRect block = new GRect(x0, y0, BRICK_WIDTH, BRICK_HEIGHT);
			add(block);
			x0 += BRICK_WIDTH;
			i--;
		}
	}
}
