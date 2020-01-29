
/*
 * File: HangmanCanvas.java
 * ------------------------
 * This file keeps track of the Hangman display.
 */

import acm.graphics.*;

public class HangmanCanvas extends GCanvas {
	/* Constants for the simple version of the picture (in pixels) */
	private static final int SCAFFOLD_HEIGHT = 360;
	private static final int BEAM_LENGTH = 144;
	private static final int ROPE_LENGTH = 18;
	private static final int HEAD_RADIUS = 36;
	private static final int BODY_LENGTH = 144;
	private static final int ARM_OFFSET_FROM_HEAD = 28;
	private static final int UPPER_ARM_LENGTH = 72;
	private static final int LOWER_ARM_LENGTH = 44;
	private static final int HIP_WIDTH = 36;
	private static final int LEG_LENGTH = 108;
	private static final int FOOT_LENGTH = 28;

	private GLine scaffold;
	private GLine beam;
	private GLine rope;
	private GLabel enteredWord = new GLabel("");
	private GLabel incGuessList = new GLabel("");
	private String incLetters = "";
	private int incGuessNum;
	private GOval head;
	private GLine body;

	/** Resets the display so that only the scaffold appears */
	public void reset() {
		removeAll();
		resetCounters();
		drawScaffold();
	}

	private void resetCounters() {
		incLetters = "";
		incGuessNum = 0;
	}

	private void drawScaffold() {
		drawGallows();
		drawBeam();
		drawRope();
	}

	private void drawGallows() {
		double scaffoldX = getWidth() / 2 - BEAM_LENGTH;
		double scaffoldY = getHeight() / 2 + 150;
		scaffold = new GLine(scaffoldX, scaffoldY, scaffoldX, scaffoldY - SCAFFOLD_HEIGHT);
		add(scaffold);
	}

	private void drawBeam() {
		double beamX = scaffold.getEndPoint().getX();
		double beamY = scaffold.getEndPoint().getY();
		beam = new GLine(beamX, beamY, beamX + BEAM_LENGTH, beamY);
		add(beam);
	}

	private void drawRope() {
		double ropeX = beam.getEndPoint().getX();
		double ropeY = beam.getEndPoint().getY();
		rope = new GLine(ropeX, ropeY, ropeX, ropeY + ROPE_LENGTH);
		add(rope);
	}

	/**
	 * Updates the word on the screen to correspond to the current state of the
	 * game. The argument string shows what letters have been guessed so far;
	 * unguessed letters are indicated by hyphens.
	 */
	public void displayWord(String word) {
		enteredWord.setLabel(word);
		enteredWord.setLocation(scaffold.getStartPoint().getX(), scaffold.getStartPoint().getY() + 40);
		enteredWord.setFont("LONDON-35");
		add(enteredWord);
	}

	/**
	 * Updates the display to correspond to an incorrect guess by the user. Calling
	 * this method causes the next body part to appear on the scaffold and adds the
	 * letter to the list of incorrect guesses that appears at the bottom of the
	 * window.
	 */
	public void noteIncorrectGuess(char letter) {
		incGuessNum++;
		updateIncGuessList(letter);
		drawBodyParts();
	}

	private void updateIncGuessList(char letter) {
		incLetters += letter;
		incGuessList.setLabel(incLetters);
		incGuessList.setLocation(enteredWord.getX(), enteredWord.getY() + 25);
		incGuessList.setFont("SanSerif-20");
		add(incGuessList);
	}

	private void drawBodyParts() {
		if (incGuessNum == 1) {
			drawHead();
		} else if (incGuessNum == 2) {
			drawBody();
		} else if (incGuessNum == 3 || incGuessNum == 4) {
			drawArm();
		} else if (incGuessNum == 5 || incGuessNum == 6) {
			drawLeg();
		} else if (incGuessNum == 7 || incGuessNum == 8) {
			drawFoot();
		}
	}

	private void drawHead() {
		double headX = rope.getEndPoint().getX() - HEAD_RADIUS;
		double headY = rope.getEndPoint().getY();
		double diam = 2 * HEAD_RADIUS;
		head = new GOval(headX, headY, diam, diam);
		add(head);
	}

	private void drawBody() {
		double bodyX = head.getX() + HEAD_RADIUS;
		double bodyY = head.getY() + 2 * HEAD_RADIUS;
		body = new GLine(bodyX, bodyY, bodyX, bodyY + BODY_LENGTH);
		add(body);
	}

	private void drawArm() {
		double armX0 = body.getX();
		double armY0 = head.getY() + 2 * HEAD_RADIUS + ARM_OFFSET_FROM_HEAD;
		if (incGuessNum == 3) {
			drawLeftArm(armX0, armY0);
		} else if (incGuessNum == 4) {
			drawRightArm(armX0, armY0);
		}
	}

	private void drawLeftArm(double X0, double Y0) {
		GLine arm = new GLine(X0, Y0, X0 - UPPER_ARM_LENGTH, Y0);
		double handX = arm.getEndPoint().getX();
		GLine hand = new GLine(handX, Y0, handX, Y0 + LOWER_ARM_LENGTH);
		add(arm);
		add(hand);
	}

	private void drawRightArm(double X0, double Y0) {
		GLine arm = new GLine(X0, Y0, X0 + UPPER_ARM_LENGTH, Y0);
		double handX = arm.getEndPoint().getX();
		GLine hand = new GLine(handX, Y0, handX, Y0 + LOWER_ARM_LENGTH);
		add(arm);
		add(hand);
	}

	private void drawLeg() {
		double hipX0 = body.getEndPoint().getX();
		double hipY0 = body.getEndPoint().getY();
		if (incGuessNum == 5) {
			drawLeftLeg(hipX0, hipY0);
		} else if (incGuessNum == 6) {
			drawRightLeg(hipX0, hipY0);
		}
	}

	private void drawLeftLeg(double X0, double Y0) {
		GLine hip = new GLine(X0, Y0, X0 - HIP_WIDTH / 2, Y0);
		double legX0 = hip.getEndPoint().getX();
		GLine leg = new GLine(legX0, Y0, legX0, Y0 + LEG_LENGTH);
		add(hip);
		add(leg);
	}

	private void drawRightLeg(double X0, double Y0) {
		GLine hip = new GLine(X0, Y0, X0 + HIP_WIDTH / 2, Y0);
		double legX0 = hip.getEndPoint().getX();
		GLine leg = new GLine(legX0, Y0, legX0, Y0 + LEG_LENGTH);
		add(hip);
		add(leg);
	}

	private void drawFoot() {
		double footY0 = body.getEndPoint().getY() + LEG_LENGTH;
		if (incGuessNum == 7) {
			drawLeftFoot(footY0);
		} else if (incGuessNum == 8) {
			drawRightFoot(footY0);
		}
	}

	private void drawLeftFoot(double Y0) {
		double X0 = body.getX() - HIP_WIDTH / 2;
		GLine foot = new GLine(X0, Y0, X0 - FOOT_LENGTH, Y0);
		add(foot);
	}

	private void drawRightFoot(double Y0) {
		double X0 = body.getX() + HIP_WIDTH / 2;
		GLine foot = new GLine(X0, Y0, X0 + FOOT_LENGTH, Y0);
		add(foot);
	}
}