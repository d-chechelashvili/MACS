
/*
 * File: HangmanCanvasExtension.java
 * ------------------------
 * This file keeps track of the Hangman display.
 */

import java.applet.*;
import java.awt.*;
import java.awt.*;
import java.awt.event.*;

import acm.graphics.*;
import acm.util.MediaTools;

public class HangmanCanvasExtension extends GCanvas {

	private GLabel enteredWord = new GLabel("");
	private GLabel incGuessList = new GLabel("");
	private String incLetters = "";
	private int incGuessNum;
	private GImage losingGif;
	private GImage winningGif;
	private GRect firstRect;
	private GRect secondRect;
	private GRect thirdRect;
	private GRect fourthRect;
	private GRect fifthRect;
	private GRect sixthRect;
	private GRect seventhRect;
	private GRect eighthRect;
	private double rWidth = 1200 / 6;
	private double rHeight = 1400 / 6;

	private void initialize(int gameCount, int i) {
		if (gameCount == 1) {
			GImage im = new GImage("images/hangedMan.jpg");
			add(im);
		} else {
			if (i == 0) {
				GImage im = new GImage("images/im2.jpg");
				add(im);
			} else if (i == 1) {
				GImage im = new GImage("images/im3.jpg");
				add(im);
			} else if (i == 2) {
				GImage im = new GImage("images/im4.jpg");
				add(im);
			} else if (i == 3) {
				GImage im = new GImage("images/im5.jpg");
				add(im);
			} else if (i == 4) {
				GImage im = new GImage("images/im6.jpg");
				add(im);
			} else if (i == 5) {
				GImage im = new GImage("images/im7.jpg");
				add(im);
			} else if (i == 6) {
				GImage im = new GImage("images/im8.png");
				add(im);
			} else if (i == 7) {
				GImage im = new GImage("images/im9.jpg");
				add(im);
			} else if (i == 8) {
				GImage im = new GImage("images/im10.png");
				add(im);
			} else if (i == 9) {
				GImage im = new GImage("images/im11.png");
				add(im);
			} else if (i == 10) {
				GImage im = new GImage("images/im12.jpg");
				add(im);
			} else if (i == 11) {
				GImage im = new GImage("images/im13.jpg");
				add(im);
			} else if (i == 12) {
				GImage im = new GImage("images/im14.png");
				add(im);
			} else if (i == 13) {
				GImage im = new GImage("images/im15.jpg");
				add(im);
			} else if (i == 14) {
				GImage im = new GImage("images/im16.png");
				add(im);
			} else if (i == 15) {
				GImage im = new GImage("images/im17.jpg");
				add(im);
			} else if (i == 16) {
				GImage im = new GImage("images/im18.jpg");
				add(im);
			} else if (i == 17) {
				GImage im = new GImage("images/im19.jpg");
				add(im);
			} else if (i == 18) {
				GImage im = new GImage("images/im20.jpg");
				add(im);
			} else if (i == 19) {
				GImage im = new GImage("images/im21.png");
				add(im);
			} else if (i == 20) {
				GImage im = new GImage("images/im22.jpg");
				add(im);
			} else if (i == 21) {
				GImage im = new GImage("images/im23.jpg");
				add(im);
			} else if (i == 22) {
				GImage im = new GImage("images/im24.jpg");
				add(im);
			} else if (i == 23) {
				GImage im = new GImage("images/im25.jpg");
				add(im);
			} else if (i == 24) {
				GImage im = new GImage("images/im26.png");
				add(im);
			} else if (i == 25) {
				GImage im = new GImage("images/im27.jpg");
				add(im);
			} else if (i == 26) {
				GImage im = new GImage("images/im28.jpg");
				add(im);
			} else if (i == 27) {
				GImage im = new GImage("images/im29.png");
				add(im);
			} else if (i == 28) {
				GImage im = new GImage("images/im30.jpg");
				add(im);
			}
		}
		coverScreen();
	}

	private void coverScreen() {
		fillFirstColumn();
		fillSecondColumn();
		fillThirdColumn();
	}

	private void fillFirstColumn() {
		firstRect = new GRect(rWidth, rHeight);
		secondRect = new GRect(0, rHeight, rWidth, rHeight);
		thirdRect = new GRect(0, 2 * rHeight, rWidth, rHeight);
		firstRect.setFilled(true);
		secondRect.setFilled(true);
		thirdRect.setFilled(true);
		firstRect.setColor(Color.white);
		secondRect.setColor(Color.white);
		thirdRect.setColor(Color.white);
		add(firstRect);
		add(secondRect);
		add(thirdRect);
	}

	private void fillSecondColumn() {
		fourthRect = new GRect(rWidth, 0, rWidth, rHeight);
		fifthRect = new GRect(rWidth, 2 * rHeight, rWidth, rHeight);
		fourthRect.setFilled(true);
		fifthRect.setFilled(true);
		fourthRect.setColor(Color.white);
		fifthRect.setColor(Color.white);
		add(fourthRect);
		add(fifthRect);
	}

	private void fillThirdColumn() {
		sixthRect = new GRect(2 * rWidth, 0, rWidth, rHeight);
		seventhRect = new GRect(2 * rWidth, rHeight, rWidth, rHeight);
		eighthRect = new GRect(2 * rWidth, 2 * rHeight, rWidth, rHeight);
		sixthRect.setFilled(true);
		seventhRect.setFilled(true);
		eighthRect.setFilled(true);
		sixthRect.setColor(Color.white);
		seventhRect.setColor(Color.white);
		eighthRect.setColor(Color.white);
		add(sixthRect);
		add(seventhRect);
		add(eighthRect);
	}

	/** Resets the display so that only the scaffold appears */
	public void reset(int gameCount, int index) {
		removeAll();
		resetCounters();
		initialize(gameCount, index);
	}

	private void resetCounters() {
		incLetters = "";
		incGuessNum = 0;
	}

	/**
	 * Updates the word on the screen to correspond to the current state of the
	 * game. The argument string shows what letters have been guessed so far;
	 * unguessed letters are indicated by hyphens.
	 */
	public void displayWord(String word) {
		enteredWord.setLabel(word);
		enteredWord.setLocation(getWidth() / 2 + 20, 500);
		enteredWord.setFont("LONDON-30-ITALIC");
		enteredWord.setColor(Color.RED);
		if (enteredWord.getX() + enteredWord.getWidth() > getWidth()) {
			enteredWord.setLocation(getWidth() - enteredWord.getWidth(), enteredWord.getY());
		}
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
		drawImageParts();
	}

	private void updateIncGuessList(char letter) {
		incLetters += letter;
		incGuessList.setLabel(incLetters);
		incGuessList.setLocation(enteredWord.getX(), enteredWord.getY() + 25);
		incGuessList.setFont("SanSerif-25-BOLD");
		incGuessList.setColor(Color.CYAN);
		add(incGuessList);
	}

	private void drawImageParts() {
		if (incGuessNum == 1) {
			moveAway(firstRect);
		} else if (incGuessNum == 2) {
			moveAway(secondRect);
		} else if (incGuessNum == 3) {
			moveAway(thirdRect);
		} else if (incGuessNum == 4) {
			moveAway(fourthRect);
		} else if (incGuessNum == 5) {
			moveAway(fifthRect);
		} else if (incGuessNum == 6) {
			moveAway(sixthRect);
		} else if (incGuessNum == 7) {
			moveAway(seventhRect);
		} else if (incGuessNum == 8) {
			moveAway(eighthRect);
		}
	}

	private void moveAway(GRect rect) {
		if (incGuessNum == 1 || incGuessNum == 2 || incGuessNum == 3) {
			while (rect.getX() > 0 - rect.getWidth()) {
				rect.pause(15);
				rect.move(-3.0, 0);
			}
		} else if (incGuessNum == 4) {
			while (rect.getY() > 0 - rect.getHeight()) {
				rect.pause(15);
				rect.move(0, -3.0);
			}
		} else if (incGuessNum == 5) {
			while (rect.getY() < getHeight() + rect.getHeight()) {
				rect.pause(15);
				rect.move(0, 3.0);
			}
		} else if (incGuessNum == 6 || incGuessNum == 7 || incGuessNum == 8) {
			while (rect.getX() < getWidth() + rect.getWidth()) {
				rect.pause(15);
				rect.move(3.0, 0);
			}
		}
	}

	public void playLosingGif(int n) {
		AudioClip losingClip = MediaTools.loadAudioClip("gifs/losing/losingAu.au");
		if (n == 0) {
			losingGif = new GImage("gifs/losing/losingGif0.gif");
		} else if (n == 1) {
			losingGif = new GImage("gifs/losing/losingGif1.gif");
		} else if (n == 2) {
			losingGif = new GImage("gifs/losing/losingGif2.gif");
		} else if (n == 3) {
			losingGif = new GImage("gifs/losing/losingGif3.gif");
		} else if (n == 4) {
			losingGif = new GImage("gifs/losing/losingGif4.gif");
		} else if (n == 5) {
			losingGif = new GImage("gifs/losing/losingGif5.gif");
		} else if (n == 6) {
			losingGif = new GImage("gifs/losing/losingGif6.gif");
		}
		losingClip.play();
		add(losingGif, getWidth() / 2 - losingGif.getWidth() / 2, getHeight() / 2 - losingGif.getHeight() / 2);
	}

	public void playWinningGif(int n) {
		AudioClip winningClip = MediaTools.loadAudioClip("gifs/winning/winningAu.au");
		if (n == 0) {
			winningGif = new GImage("gifs/winning/winningGif0.gif");
		} else if (n == 1) {
			winningGif = new GImage("gifs/winning/winningGif1.gif");
		} else if (n == 2) {
			winningGif = new GImage("gifs/winning/winningGif2.gif");
		} else if (n == 3) {
			winningGif = new GImage("gifs/winning/winningGif3.gif");
		} else if (n == 4) {
			winningGif = new GImage("gifs/winning/winningGif4.gif");
		} else if (n == 5) {
			winningGif = new GImage("gifs/winning/winningGif5.gif");
		} else if (n == 6) {
			winningGif = new GImage("gifs/winning/winningGif6.gif");
		}
		winningClip.play();
		add(winningGif, getWidth() / 2 - winningGif.getWidth() / 2, getHeight() / 2 - winningGif.getHeight() / 2);
	}
}