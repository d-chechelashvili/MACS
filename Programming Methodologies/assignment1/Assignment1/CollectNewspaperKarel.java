/*
 * File: CollectNewspaperKarel.java
 * --------------------------------
 * At present, the CollectNewspaperKarel subclass does nothing.
 * Your job in the assignment is to add the necessary code to
 * instruct Karel to walk to the door of its house, pick up the
 * newspaper (represented by a beeper, of course), and then return
 * to its initial position in the upper left corner of the house.
 */

import stanford.karel.*;
public class CollectNewspaperKarel extends SuperKarel {

//Karel is in the house , we want to pick the newspaper which is in front of the door and return to starting point.
	
	public void run() {


		moveToTheDoor();
		pickTheNewspaper();
		returnBack();

	}

	private void moveToTheDoor() {

		moveToTheFrontWall();
		moveToTheEntrance();

	}

	private void moveToTheFrontWall() {
		// precondition - top left corner of the house >
		// postcondition - top right corner of the house >

		while (frontIsClear()) {
			move();
		}
	}

	private void moveToTheEntrance() {
		// pre - see "post" above
		// post - at the entrance >

		turnRight();
		while (leftIsBlocked()) {
			move();

		}
		turnLeft();

	}

	private void pickTheNewspaper() {
		move();
		pickBeeper();

	}

	private void returnBack() {
		moveToTheBackWall();
		moveToTheStartingPoint();
	}

	private void moveToTheBackWall() {
		// pre - in front of the entrance >
		// post - facing the back wall <
		turnAround();
		while (frontIsClear()) {
			move();
		}

	}

	private void moveToTheStartingPoint() {
		// pre - see "post" above
		// post - top left corner >
		turnRight();
		while (frontIsClear()) {
			move();
		}
		turnRight();

	}
}