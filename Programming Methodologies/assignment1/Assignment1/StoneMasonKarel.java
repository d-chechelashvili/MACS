
/*
 * File: StoneMasonKarel.java
 * --------------------------
 *	Karel must repair every pillar( which consist of beepers) in this world. 
 *	pillars are distanced from each other by 4 columns  
 *	Karel mustn't put more than one beeper on one cell.
 *	on the last column of the world there is a pillar.
 */

import stanford.karel.*;

public class StoneMasonKarel extends SuperKarel {
	public void run() {
		repairAllPillars();
		getToFloor();

	}

	private void repairAllPillars() {

		// lets repair pillars one by one until none is left

		while (frontIsClear()) {
			repairPillar();
			// repairPillar means repairing pillar + moving in position to repair new one
		}
		finishOneColumn();
	}

	private void repairPillar() {

		finishOneColumn();
		moveToStartNew();

	}

	private void finishOneColumn() {
		// pre - (x,1) >
		// post - (x,y) facing up (needs turnRight for >)

		// to prevent OBOB on first row
		if (beepersPresent()) {
			turnLeft();
		} else {
			putBeeper();
			turnLeft();
		}

		while (frontIsClear()) {

			if (beepersPresent()) {
				move();
			} else {
				putBeeper();
				move();
			}

		}
		if (noBeepersPresent()) {
			putBeeper();
		}

	}

	private void moveToStartNew() {

		getToFloor();
		moveToNextPillar();

	}

	private void getToFloor() {
		// precondition - (x,y) facing up (needs turnRight for >)
		// postcondition - (x,1) >

		turnAround();
		while (frontIsClear()) {
			move();
		}
		turnLeft();
	}

	private void moveToNextPillar() {
		// pre - (x,1) >
		// post - (x+4,1) >
		move();
		move();
		move();
		move();
	}

}
