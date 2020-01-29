/*
 * File: PythagoreanTheorem.java
 * Name: 
 * Section Leader: 
 * -----------------------------
 * This file is the starter file for the PythagoreanTheorem problem.
 */

import acm.program.*;

public class PythagoreanTheorem extends ConsoleProgram {
	public void run() {
		println("Enter values to compute Pythagorean theorem");
		int a = readInt("a: ");
		int b = readInt("b: ");
		//because sides of a right triangle must be more than zero
		if (a >= 0 && b >= 0) {
			println("c: " + Hypotenuse(a, b));
		} else {
			println("Not Funny");
		}
	}
	//calculates the hypotenuse for a right triangle with a and b as other two sides
	private double Hypotenuse(int a, int b) {
		return Math.sqrt(Math.pow(a, 2) + Math.pow(b, 2));
	}
}
