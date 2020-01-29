
/*
 * File: FindRange.java
 * Name: 
 * Section Leader: 
 * --------------------
 * This file is the starter file for the FindRange problem.
 */

import acm.program.*;

public class FindRange extends ConsoleProgram {
	private static final int trigger = 36;

	public void run() {
		println("This program finds the largest and smallest numbers");
		int smallest = 0; // we need some base value for both integers
		int largest = 0;
		for (int i = 0; i >= 0; i++) {
			int n = readInt("? ");
			// if the first number entered is trigger this program should display Bingo
			if (n == trigger && i == 0) {
				println("Bingo");
				break;
			}
			// whenever the trigger is entered program should display
			// largest and smallest numbers entered and it should end
			if (n == trigger) {
				println("smallest: " + smallest);
				println("largest: " + largest);
				break;
			}
			if (i == 0) {
				smallest = n;
				largest = n;
			}
			smallest = Math.min(n, smallest);
			largest = Math.max(n, largest);
		}
	}
}
