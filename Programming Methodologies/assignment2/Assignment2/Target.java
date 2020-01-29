
/*
 * File: Target.java
 * Name: 
 * Section Leader: 
 * -----------------
 * This file is the starter file for the Target problem.
 */

import acm.graphics.*;
import acm.program.*;
import java.awt.*;
/*this program should draw an archery target which consists of 3 circles
 * its centre is located at the centre of the screen */

public class Target extends GraphicsProgram {
	// with each layer of circles difference between radiuses is 0.89cm
	private static final double Rdifference = 0.89;

	public void run() {
		double R = 2.54; // radius of the biggest circle
		for (int i = 0; i < 3; i++) {
			double Rpix = ConvertToPixel(R); // converts R which is given in cms to pixels
			double x = getWidth() / 2 - Rpix;
			double y = getHeight() / 2 - Rpix;
			double d = 2 * Rpix; // d stands for diameter
			GOval Circle = new GOval(x, y, d, d);
			Circle.setFilled(true);
			if (i % 2 == 0) {
				Circle.setColor(Color.red);
			} else {
				Circle.setColor(Color.white);
			}
			add(Circle);
			R -= Rdifference;
		}
	}

	// returns the value of R converted to pixels.
	private double ConvertToPixel(double R) {
		return (R * 72) / 2.54;
	}
}