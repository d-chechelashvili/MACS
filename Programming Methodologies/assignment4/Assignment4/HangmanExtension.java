
/*
 * File: HangmanExtension.java
 * ------------------
 * This program will eventually play the Hangman game from
 * Assignment #4.
 */

import acm.graphics.*;
import acm.program.*;
import acm.util.*;

import java.awt.*;

public class HangmanExtension extends ConsoleProgram {
	private HangmanCanvasExtension canvasExtension;
	private RandomGenerator rand = RandomGenerator.getInstance();

	private static final int APP_WIDTH = 1072;
	private static final int APP_HEIGHT = 625;

	private int guesses = 8;

	private char letter;
	private int guessCount;
	private int gameCount = 1;
	HangmanLexicon lexicon = new HangmanLexicon();
	HangmanLexiconExtension1 lexicon1 = new HangmanLexiconExtension1();
	private int index = rand.nextInt(lexicon.getWordCount());
	private int index1;
	private boolean result;
	private String hiddenWord = chooseWord();
	private String enteredWord = startingWord();
	private boolean space;

	public void init() {
		setSize(APP_WIDTH, APP_HEIGHT);
		canvasExtension = new HangmanCanvasExtension();
		add(canvasExtension);
		greetPlayer();
	}

	private String chooseWord() {
		String s = "";
		if (gameCount == 1) {
			s = lexicon.getWord(index);
			println(s);
		} else {
			s = lexicon1.getWordi(index1);
		}
		
		return s;
	}

	private String startingWord() {
		String s = "";
		for (int i = 0; i < hiddenWord.length(); i++) {
			if (Character.isLetter(hiddenWord.charAt(i))) {
				s += "-";
			} else {
				s += " ";
				space = true;
			}
		}
		return s;
	}

	public void run() {
		if (gameCount == 1) {
			canvasExtension.reset(gameCount, index);
		} else {
			canvasExtension.reset(gameCount, index1);
		}
		playGame();
		checkResult();
		newGame();
	}

	private void greetPlayer() {
		println("Welcome to Hangman !");
	}

	private void playGame() {
		while (guesses > 0 && guessCount < hiddenWord.length()) {
			checkWordState();
			printGuessCount();
			checkGuess();
		}
	}

	private void checkWordState() {
		println("The word now looks like this: " + enteredWord);
		canvasExtension.displayWord(enteredWord);
	}

	private void printGuessCount() {
		println("You have " + guesses + " guesses left");
	}

	private void checkGuess() {
		String s = readLine("Your Guess: ");
		if (s.length() == 1) {
			letter = s.charAt(0);
			if (((letter >= 'a' && letter <= 'z') || (letter >= 'A' && letter <= 'Z'))) {
				convertCase();
				if (!duplicate()) {
					checkMatch();
				} else {
					println("That letter is already revealed");
				}
			} else {
				println("That guess is illegal.Try again");
			}
		} else {
			println("Enter a letter please");
		}
	}

	private void convertCase() {
		if (Character.isLowerCase(letter)) {
			letter = Character.toUpperCase(letter);
		}
	}

	private boolean duplicate() {
		for (int i = 0; i < enteredWord.length(); i++) {
			if (enteredWord.charAt(i) == letter) {
				return true;
			}
		}
		return false;
	}

	private void checkMatch() {
		if (space) {
			for (int i = 0; i < hiddenWord.length(); i++) {
				if (hiddenWord.charAt(i) == ' ') {
					guessCount++;
				}
			}
			space = false;
		}
		if (letterMatch()) {
			println("That guess is correct");
		} else {
			println("There are no " + letter + "'s in the word");
			canvasExtension.noteIncorrectGuess(letter);
			guesses--;
		}
	}

	private boolean letterMatch() {
		boolean match = false;
		for (int i = 0; i < hiddenWord.length(); i++) {
			if (hiddenWord.charAt(i) == letter) {
				enteredWord = enteredWord.substring(0, i) + letter + enteredWord.substring(i + 1);
				match = true;
				guessCount++;
			}
		}
		return match;
	}

	private void checkResult() {
		if (guesses == 0) {
			result = false;
			println("You're completely hung.");
			println("The word was: " + hiddenWord);
			println("You lose.");
			pause(300);
			canvasExtension.playLosingGif(rand.nextInt(7));
		} else {
			result = true;
			canvasExtension.displayWord(enteredWord);
			println("You guessed the word	: " + hiddenWord);
			println("You win.");
			pause(300);
			canvasExtension.playWinningGif(rand.nextInt(7));
		}
	}

	private void newGame() {
		if (result) {
			String trigger = readLine(
					"Enter 0 if you want to play again \nEnter 1 if you want to play second round \nEnter any other symbol to stop playing: ");
			if (trigger.length() == 1) {
				if (trigger.charAt(0) == '0') {
					reset();
					run();
				} else if (trigger.charAt(0) == '1') {
					println("Welcome to the round of video games! ");
					secondRound();
				} else {
					println("Goodbye!");
				}
			} else {
				println("Goodbye!");
			}
		} else {
			String trigger = readLine("Enter any symbol if you want to play again: ");
			if (trigger.length() > 0) {
				reset();
				run();
			} else {
				println("Goodbye!");
			}
		}
	}

	private void secondRound() {
		index1 =rand.nextInt(lexicon1.getWordCounti());
		gameCount++;
		hiddenWord = chooseWord();
		enteredWord = startingWord();
		guesses = 8;
		guessCount = 0;
		run();
	}

	private void reset() {
		index = rand.nextInt(lexicon.getWordCount());
		hiddenWord = chooseWord();
		enteredWord = startingWord();
		guesses = 8;
		guessCount = 0;
		gameCount = 1;
		result = false;
	}
}