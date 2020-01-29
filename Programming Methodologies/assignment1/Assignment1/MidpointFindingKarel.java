
/*
 * File: MidpointFindingKarel.java
 * -------------------------------
 * When you finish writing it, the MidpointFindingKarel class should
 * leave a beeper on the corner closest to the center of 1st Street
 * (or either of the two central corners if 1st Street has an even
 * number of corners).  Karel can put down additional beepers as it
 * looks for the midpoint, but must pick them up again before it
 * stops.  The world may be of any size, but you are allowed to
 * assume that it is at least as tall as it is wide.
 */

import stanford.karel.*;

public class MidpointFindingKarel extends SuperKarel {

	public void run() {
		findTheMiddlePoint();
		putBeeper();

	}

	private void findTheMiddlePoint() {
		findTheMiddleColumn();
		moveToFirstRow();

	}

	private void findTheMiddleColumn() {

		/*
		 * because map is square if karel moves twice upwards for each movement forward
		 it will stop at the middle point of top row if map has odd number of cells,
		 if map has even number it will stop at one of the middle cells on top row.
		 */
		// precondition - (1,1) >
		// postcondition - (Ymax,Xmid) facing up

		turnLeft(); // to prevent fail at the first move.
		while (frontIsClear()) {
			turnRight();
			move();
			moveUp();

		}

	}

	private void moveUp() {
		turnLeft();
		move();
		if (frontIsClear()) {
			move();

		}

	}

	private void moveToFirstRow() {
		turnAround();
		while (frontIsClear()) {
			move();
		}

	}

}
