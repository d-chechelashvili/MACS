
/*
 * File: Hangman.java
 * ------------------
 * This program will eventually play the Hangman game from
 * Assignment #4.
 */

import acm.graphics.*;
import acm.program.*;
import acm.util.*;

import java.awt.*;

public class Hangman extends ConsoleProgram {
	private HangmanCanvas canvas;
	HangmanLexicon lexicon = new HangmanLexicon();
	private RandomGenerator rand = RandomGenerator.getInstance();

	private String hiddenWord = chooseWord();
	private int guesses = 8;
	private String enteredWord = startingWord();
	private char letter;
	private int guessCount;

	public void init() {
		canvas = new HangmanCanvas();
		add(canvas);
		greetPlayer();
	}

	// chooses a word from lexicon
	private String chooseWord() {
		int i = rand.nextInt(lexicon.getWordCount());
		String s = lexicon.getWord(i);
		return s;
	}

	// ciphers chosen word
	private String startingWord() {
		String s = "";
		for (int i = 0; i < hiddenWord.length(); i++) {
			s += "-";
		}
		return s;
	}

	public void run() {
		canvas.reset();
		playGame();
		checkResult();
		newGame();
	}

	private void greetPlayer() {
		println("Welcome to Hangman !");
	}

	// controls the gameplay
	private void playGame() {
		while (guesses > 0 && guessCount < hiddenWord.length()) {
			checkWordState();
			printGuessCount();
			checkGuess();
		}
	}

	// shows user state of the ciphered word
	private void checkWordState() {
		println("The word now looks like this: " + enteredWord);
		canvas.displayWord(enteredWord);
	}

	// tells user how many guesses are left
	private void printGuessCount() {
		println("You have " + guesses + " guesses left");
	}

	// checks if entered letter is legal and then if it matches one of hiddenWord's
	// letters
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

	// converts entered letter to uppercase if its lowercase
	private void convertCase() {
		if (Character.isLowerCase(letter)) {
			letter = Character.toUpperCase(letter);
		}
	}

	// checks if entered letter is new and user hasn't revealed it already
	private boolean duplicate() {
		for (int i = 0; i < enteredWord.length(); i++) {
			if (enteredWord.charAt(i) == letter) {
				return true;
			}
		}
		return false;
	}

	// checks if entered letter is in hidden word
	private void checkMatch() {
		if (letterMatch()) {
			println("That guess is correct");
		} else {
			println("There are no " + letter + "'s in the word");
			canvas.noteIncorrectGuess(letter);
			guesses--;
		}
	}

	// checks if entered letter is part of hidden word , if yes reveals it in
	// ciphered word
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

	// displays the hidden word and tells user whether he/she won or lost
	private void checkResult() {
		if (guesses == 0) {
			println("You're completely hung.");
			println("The word was: " + hiddenWord);
			println("You lose.");
		} else {
			canvas.displayWord(enteredWord);
			println("You guessed the word	: " + hiddenWord);
			println("You win.");
		}
	}

	// starts a new game if user wants to
	private void newGame() {
		String trigger = readLine("Enter any symbol if you want to play again: ");
		if (trigger.length() > 0) {
			reset();
			run();
		} else {
			println("Goodbye!");
		}
	}

	// prepares for new game. chooses word from lexicon,ciphers it and gives user 8
	// lives
	private void reset() {
		hiddenWord = chooseWord();
		enteredWord = startingWord();
		guesses = 8;
		guessCount = 0;
	}
}