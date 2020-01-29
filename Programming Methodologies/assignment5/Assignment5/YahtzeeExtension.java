import acm.io.*;
import acm.program.*;
import acm.util.*;

public class YahtzeeExtension extends GraphicsProgram implements YahtzeeConstants {

	public static void main(String[] args) {
		new YahtzeeExtension().start(args);
	}

	// registers the players
	public void run() {
		dialog.println("WELCOME TO GRAND DUEL!");
		dialog.println("There are 2 players in grand duel, each controls 2 scorecards");
		dialog.println(
				"Scores in second scorecard are doubled(except yahtzee)\n every time yahtzee is rolled it's score is increased by 50");
		nPlayers = 2 * PLAYERS_EXTENSION;
		initialize(dialog);
	}

	// sets up everything so that game can start
	private void initialize(IODialog dialog) {
		playerNames = new String[nPlayers];
		for (int i = 1; i <= nPlayers / 2; i++) {
			playerNames[i - 1] = dialog.readLine("Enter name for player " + i);
		}
		playerNames[2] = playerNames[1];
		playerNames[3] = playerNames[1];
		playerNames[1] = playerNames[0];
		display = new YahtzeeDisplay(getGCanvas(), playerNames);
		playGame();
	}

	// loop that provides gameplay and end result
	private void playGame() {
		while (roundCount > 0) {
			play();
			roundCount--;
		}
		endGame();
		determineWinner();
	}

	// controls the gameplay of each round
	private void play() {
		for (int player = 1; player <= nPlayers / 2; player++) {
			rollDice(player);
			for (int i = 0; i < N_REROLL; i++) {
				rollAgain();
			}
			chooseCategory(player);
		}
	}

	private void rollDice(int player) {
		display.printMessage(player + "'s turn!Click \"Roll Dice\" button to roll the dice.");
		display.waitForPlayerToClickRoll(player + 1);
		callDisplayDice(0);
	}

	// rolls the dice and displays each roll
	private void callDisplayDice(int rollAgain) {
		int roll;
		if (rollAgain == 0) {
			for (int i = 0; i < N_DICE; i++) {
				roll = rgen.nextInt(1, 6);
				dice[i] = roll;
			}

		} else if (rollAgain == 1) {
			for (int i = 0; i < N_DICE; i++) {
				if (display.isDieSelected(i)) {
					roll = rgen.nextInt(1, 6);
					dice[i] = roll;
				}
			}
		}
		try {
			display.displayDice(dice);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	// method that gives player opportunity to roll any number of dices again
	private void rollAgain() {
		display.printMessage("Select the dice you wish to re-roll and click \"Roll Again\".");
		display.waitForPlayerToSelectDice();
		callDisplayDice(1);
	}

	// this method lets player choose category in which his/her roll belongs, only
	// when chosen category's cell is empty
	private void chooseCategory(int player) {
		display.printMessage("Select a category for this roll.");
		int category = display.waitForPlayerToSelectCategory();
		int scorecard = dialog.readInt("which scorecard?");
		if ((scorecard == 1 || scorecard == 2) && (!filled[player][scorecard][category])) {
			checkCategory(player, category, scorecard);
			if (player == 1) {
				display.updateScorecard(category, player + scorecard - 1, score);
			} else {
				display.updateScorecard(category, player + scorecard, score);
			}
			updatetotals(category, player, scorecard);
		} else {
			chooseCategory(player);
		}
	}

	// checks if player chose upper part of the table or lower
	private void checkCategory(int player, int category, int scorecard) {
		if (category < UPPER_SCORE || category == CHANCE) {
			checkUpper(category, scorecard);
		} else {
			checkLower(category, scorecard);
		}
		filled[player][scorecard][category] = true;
	}

	// writes score in category's cell corresponding to chosen category and rolled
	// dice
	private void checkUpper(int category, int scorecard) {
		int n = 0;
		switch (category) {
		case ONES:
			n = 1;
			calculateScore(n, scorecard);
			break;
		case TWOS:
			n = 2;
			calculateScore(n, scorecard);
			break;
		case THREES:
			n = 3;
			calculateScore(n, scorecard);
			break;
		case FOURS:
			n = 4;
			calculateScore(n, scorecard);
			break;
		case FIVES:
			n = 5;
			calculateScore(n, scorecard);
			break;
		case SIXES:
			n = 6;
			calculateScore(n, scorecard);
			break;
		case CHANCE:
			sumAll(scorecard);
			break;
		}
	}

	private void calculateScore(int n, int scorecard) {
		if (scorecard == 2) {
			for (int i = 0; i < N_DICE; i++) {
				if (dice[i] == n) {
					score += 2 * n;
				}
			}
		} else {
			for (int i = 0; i < N_DICE; i++) {
				if (dice[i] == n) {
					score += n;
				}
			}
		}
	}

	// sums numbers on all dice
	private void sumAll(int scorecard) {
		if (scorecard == 2) {
			for (int i = 0; i < N_DICE; i++) {
				score += 2 * dice[i];
			}
		} else {
			for (int i = 0; i < N_DICE; i++) {
				score += dice[i];
			}
		}
	}

	// writes score in category's cell corresponding to chosen category and rolled
	// dice
	// if rolled dice doesn't fit in the chosen category writes 0 in that category's
	// cell
	private void checkLower(int category, int scorecard) {
		int[] numCount = new int[N_SIDES + 1];
		countEachNum(numCount);
		switch (category) {
		case THREE_OF_A_KIND:
			if (toak(numCount))
				sumAll(scorecard);
			break;
		case FOUR_OF_A_KIND:
			if (foak(numCount))
				sumAll(scorecard);
			break;
		case FULL_HOUSE:
			if (FH(numCount))
				if (scorecard == 2) {
					score = 2 * FULL_HOUSE_SCORE;
				} else {
					score = FULL_HOUSE_SCORE;
				}
			break;
		case SMALL_STRAIGHT:
			if (sStraight(numCount))
				if (scorecard == 2) {
					score = 2 * SMALL_STRAIGHT_SCORE;
				} else {
					score = SMALL_STRAIGHT_SCORE;
				}
			break;
		case LARGE_STRAIGHT:
			if (lStraight(numCount))
				if (scorecard == 2) {
					score = 2 * LARGE_STRAIGHT_SCORE;
				} else {
					score = LARGE_STRAIGHT_SCORE;
				}
			break;
		case YAHTZEE:
			if (yahtzee(numCount)) {
				score = yahtzeeNum * YAHTZEE_SCORE;
				yahtzeeNum++;
			}
			break;
		}
	}

	// counts how many die show each number
	private void countEachNum(int[] numCount) {
		for (int i = 0; i < dice.length; i++) {
			numCount[dice[i]]++;
		}
	}

	// checks if rolled dices have meet necessary requirements for three of a kind
	private boolean toak(int[] numCount) {
		for (int i = 0; i < numCount.length; i++) {
			if (numCount[i] >= 3)
				return true;
		}
		return false;
	}

	// checks if rolled dices have meet necessary requirements for four of a kind
	private boolean foak(int[] numCount) {
		for (int i = 0; i < numCount.length; i++) {
			if (numCount[i] >= 4)
				return true;
		}
		return false;
	}

	//// checks if rolled dices have meet necessary requirements for full house
	private boolean FH(int[] numCount) {
		boolean toak = false;
		boolean twoak = false;
		for (int i = 0; i < numCount.length; i++) {
			if (numCount[i] == 3) {
				toak = true;
			}
			if (numCount[i] == 2) {
				twoak = true;
			}
		}
		if (twoak && toak)
			return true;
		return false;
	}

	// checks if rolled dices have meet necessary requirements for small straight
	private boolean sStraight(int[] numCount) {
		int diffRoll = calculateDiffRoll(numCount);
		if (diffRoll >= 4 && numCount[3] > 0 && numCount[4] > 0) {
			if (privateCase(numCount))
				return false;

			return true;
		}
		return false;
	}

	// calculates number of different numbers rolled
	private int calculateDiffRoll(int[] numCount) {
		int diffRoll = 0;
		for (int i = 0; i < numCount.length; i++) {
			if (numCount[i] > 0)
				diffRoll++;
		}
		return diffRoll;
	}

	// there are only 3 cases that also satisfy if statements given in sStraight
	// method so this method eliminates them
	private boolean privateCase(int[] numCount) {
		if (numCount[1] > 0 && numCount[5] > 0 && numCount[2] == 0 && numCount[6] == 0)
			return true;
		if (numCount[1] > 0 && numCount[6] > 0 && numCount[2] == 0 && numCount[5] == 0)
			return true;
		if (numCount[2] > 0 && numCount[6] > 0 && numCount[1] == 0 && numCount[5] == 0)
			return true;

		return false;
	}

	// checks if rolled dices have meet necessary requirements for large straight
	private boolean lStraight(int[] numCount) {
		int diffRoll = calculateDiffRoll(numCount);
		if (diffRoll >= 5 && numCount[2] > 0 && numCount[3] > 0 && numCount[4] > 0 && numCount[5] > 0)
			return true;
		return false;
	}

	// checks if rolled dices have meet necessary requirements for yahtzee
	private boolean yahtzee(int[] numCount) {
		for (int i = 0; i < numCount.length; i++) {
			if (numCount[i] == N_DICE)
				return true;
		}
		return false;
	}

	// updates total , upper score and lower score for given player
	private void updatetotals(int category, int player, int scorecard) {
		if (player == 2) {
			total[player + scorecard] += score;
			if (category <= SIXES) {
				upperScore[player + scorecard] += score;
			} else {
				lowerScore[player + scorecard] += score;
			}
			score = 0;
			display.updateScorecard(TOTAL, player + scorecard, total[player + scorecard]);
		} else {
			total[player + scorecard - 1] += score;
			if (category <= SIXES) {
				upperScore[player + scorecard - 1] += score;
			} else {
				lowerScore[player + scorecard - 1] += score;
			}
			score = 0;
			display.updateScorecard(TOTAL, player + scorecard - 1, total[player + scorecard - 1]);
		}
	}

	// endgame gamovitano bolos (cilkis garet)rom yvelas ertad aupdatebdes

	// shows lower score and upper score at the end of given player's turns also
	// calculates if upper score bonus should be applied
	private void endGame() {
		for (int i = 1; i <= nPlayers; i++) {
			display.updateScorecard(UPPER_SCORE, i, upperScore[i]);
			display.updateScorecard(LOWER_SCORE, i, lowerScore[i]);
			if (upperScore[i] >= 63) {
				display.updateScorecard(UPPER_BONUS, i, UPPER_BONUS_SCORE);
				total[i] += UPPER_BONUS_SCORE;
			} else
				display.updateScorecard(UPPER_BONUS, i, 0);
			display.updateScorecard(TOTAL, i, total[i]);
		}
	}

	private void determineWinner() {
		int winner = 1;
		total[1] = total[1] + total[2];
		total[2] = total[3] + total[4];
		int maxTotal = Math.max(total[1], total[2]);
		if (maxTotal == total[2]) {
			winner = 2;
		}
		display.printMessage(
				"Congratulations, " + playerNames[winner] + ", you're the winner with a total score of " + maxTotal);
	}

	/* Private instance variables */
	private int nPlayers;
	private String[] playerNames;
	private YahtzeeDisplay display;
	private RandomGenerator rgen = new RandomGenerator();
	private int roundCount = PLAYERS_EXTENSION * N_SCORING_CATEGORIES;
	private int[] dice = new int[N_DICE];
	private int score;
	private int[] upperScore = new int[MAX_PLAYERS + 1];
	private int[] lowerScore = new int[MAX_PLAYERS + 1];
	private boolean[][][] filled = new boolean[MAX_PLAYERS + 1][3][CHANCE + 1];
	private int[] total = new int[MAX_PLAYERS + 1];
	private IODialog dialog = getDialog();
	private static final int PLAYERS_EXTENSION = 2;
	private int yahtzeeNum = 1;
}