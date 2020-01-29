
/*
 * File: ProgramHierarchy.java
 * Name: 
 * Section Leader: 
 * ---------------------------
 * This file is the starter file for the ProgramHierarchy problem.
 */

import acm.graphics.*;
import acm.program.*;
import java.awt.*;

/*
 * This program draws program hierarchy diagram The Diagram is aligned in the
 * centre so that centre of the middle line (which connects Program rectangle to
 * ConsoleProgram rectangle) is also the centre of the screen.
 */
public class ProgramHierarchy extends GraphicsProgram {
	private static final int w = 200; // width of each rectangle
	private static final int h = 70; // height of each rectangle
	private static final int distance = 33; // distance between rectangles on second tier

	public void run() {
		int hgW = getWidth() / 2; // stands for half getWidth
		int hgH = getHeight() / 2; // stands for half getHeight
		drawMother(hgW, hgH);
		drawChildren(hgW, hgH);
		connectRects(hgW, hgH);
		addTexts(hgW, hgH);
	}

	// this method draws mother rectangle
	private void drawMother(int hgW, int hgH) {
		int x = hgW - w / 2;
		int y = hgH - (h + h / 2);
		GRect mother = new GRect(x, y, w, h);
		add(mother);
	}

	// this method draws children rectangle, distance between each of them is 33
	private void drawChildren(int hgW, int hgH) {
		int x = hgW - (w / 2 + distance + w);
		int y = hgH + h / 2;
		for (int i = 0; i < 3; i++) {
			GRect child = new GRect(x, y, w, h);
			x += w + distance;
			add(child);
		}
	}

	// this method draws lines that connect mother rectangle to children
	private void connectRects(int hgW, int hgH) {
		int x0 = hgW;
		int y0 = hgH - h / 2;
		int x = hgW - (w / 2 + distance + w / 2);
		int y = hgH + h / 2;
		for (int i = 0; i < 3; i++) {
			GLine connector = new GLine(x0, y0, x, y);
			add(connector);
			x += (w / 2 + distance + w / 2);
		}
	}

	private void addTexts(int hgW, int hgH) {
		addMotherText(hgW, hgH);
		addChildrenTexts(hgW, hgH);
	}

	private void addMotherText(int hgW, int hgH) {
		GLabel motherText = new GLabel("Program");
		motherText.setFont("black-20");
		double x = hgW - motherText.getWidth() / 2;
		double y = hgH - (h / 2 + h / 2 - motherText.getAscent() / 2);
		motherText.setLocation(x, y);
		add(motherText);
	}

	private void addChildrenTexts(int hgW, int hgH) {
		addChildTextOne(hgW, hgH);
		addChildTextTwo(hgW, hgH);
		addChildTextThree(hgW, hgH);
	}

	private void addChildTextOne(int hgW, int hgH) {
		GLabel childText1 = new GLabel("GraphicsProgram");
		childText1.setFont("black-20");
		double x1 = hgW - (w / 2 + distance + w / 2 + childText1.getWidth() / 2);
		double y1 = hgH + (h / 2 + h / 2 + childText1.getAscent() / 2);
		childText1.setLocation(x1, y1);
		add(childText1);
	}

	private void addChildTextTwo(int hgW, int hgH) {
		GLabel childText2 = new GLabel("ConsoleProgram");
		childText2.setFont("black-20");
		double x2 = hgW - childText2.getWidth() / 2;
		double y2 = hgH + (h / 2 + h / 2 + childText2.getAscent() / 2);
		childText2.setLocation(x2, y2);
		add(childText2);
	}

	private void addChildTextThree(int hgW, int hgH) {
		GLabel childText3 = new GLabel("DialogProgram");
		childText3.setFont("black-20");
		double x3 = hgW + (w / 2 + distance + w / 2 - childText3.getWidth() / 2);
		double y3 = hgH + (h / 2 + h / 2 + childText3.getAscent() / 2);
		childText3.setLocation(x3, y3);
		add(childText3);
	}
}
