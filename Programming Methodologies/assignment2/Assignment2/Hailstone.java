
/*
 * File: Hailstone.java
 * Name: 
 * Section Leader: 
 * --------------------
 * This file is the starter file for the Hailstone problem.
 */

import acm.program.*;

public class Hailstone extends ConsoleProgram {
	private static final int destination = 1;

	public void run() {
		int tryNum = 0; // tryNum stands for number of tries
		int n = readInt("Enter a number: ");
		// if entered number is less than 1 program shouldn't start computing anything
		if (n >= 1) {
			while (n != destination) {
				if (n % 2 == 0) {
					n = wheneven(n);
				} else {
					n = whenodd(n);
				}
				tryNum++;
			}
			println("The process took " + tryNum + " to reach " + destination);
		} else
			println("Still not funny");
	}

	private int wheneven(int n) {
		println(n + " is even, so I take half: " + n / 2);
		n /= 2;
		return n;
	}

	private int whenodd(int n) {
		println(n + " is odd, so I make 3n+1: " + (3 * n + 1));
		n = 3 * n + 1;
		return n;
	}
}
