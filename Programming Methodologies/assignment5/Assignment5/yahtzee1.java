import acm.io.*;
import acm.program.*;
import acm.util.*;
import java.io.*;
import java.util.*;

public class yahtzee1 extends GraphicsProgram implements YahtzeeConstants {

	public static void main(String[] args) {
		new Yahtzee().start(args);
	}

	int[] currentDice;
	int[][] scoreboard;

	public void run() {
		IODialog dialog = getDialog();
		nPlayers = dialog.readInt("Enter number of players");
		playerNames = new String[nPlayers];
		currentDice = new int[N_DICE];
		scoreboard = new int[N_CATEGORIES][nPlayers];
		for (int i = 1; i <= nPlayers; i++) {
			playerNames[i - 1] = dialog.readLine("Enter name for player " + i);
		}
		display = new YahtzeeDisplay(getGCanvas(), playerNames);
		createScoreBoard(scoreboard);
		playGame();
		getTheWinner();
	}

	private void createScoreBoard(int[][] scoreboard2) {
		for (int i = 0; i < N_CATEGORIES; i++) {
			for (int j = 0; j < nPlayers; j++) {
				scoreboard[i][j] = -1;
			}
		}

	}

	private void playGame() {
		for (int i = 12; i < N_SCORING_CATEGORIES; i++) {
			for (int n = 1; n <= nPlayers; n++) {
				rollTheDice(n);
				for (int j = 0; j < 2; j++) {
					rollTheDiceSecondTime();
				}

				chooseTheCategory(n);
			}
		}
		getTheResults();

	}

	private void rollTheDice(int n) {
		// [ , , , , , ]
		display.printMessage(
				playerNames[n - 1] + "'s turn , Click " + '"' + " Roll Dice " + '"' + " button to roll the dice");
		display.waitForPlayerToClickRoll(n);
		for (int k = 0; k < N_DICE; k++) {
			currentDice[k] = rgen.nextInt(1, 6);
		}
		display.displayDice(currentDice);

	}

	private void rollTheDiceSecondTime() {
		display.printMessage("select the dice you wish to re-roll and click " + '"' + " Roll Again " + '"');
		display.waitForPlayerToSelectDice();
		for (int m = 0; m < N_DICE; m++) {
			if (display.isDieSelected(m)) {
				currentDice[m] = rgen.nextInt(1, 6);
			}
		}
		display.displayDice(currentDice);
	}

	private void chooseTheCategory(int player) {

		int score = 0;
		display.printMessage("select category for this roll");

		while (true) {
			int category = display.waitForPlayerToSelectCategory();
			if (categoryIsValid(player, category)) {
				calculateTheScore(player, category);
				break;
			}
		}
	}

	private void calculateTheScore(int player, int category) {
		int score = 0;
		if (category == ONES) {
			score = 0;
			for (int i = 0; i < N_DICE; i++) {
				if (currentDice[i] == 1)
					score++;

			}
			scoreboard[ONES - 1][player - 1] = score;

		} else if (category == TWOS) {
			score = 0;
			for (int i = 0; i < N_DICE; i++) {
				if (currentDice[i] == 2)
					score = score + 2;

			}
			scoreboard[TWOS - 1][player - 1] = score;

		} else if (category == THREES) {
			score = 0;
			for (int i = 0; i < N_DICE; i++) {
				if (currentDice[i] == 3)
					score = score + 3;
			}
			scoreboard[THREES - 1][player - 1] = score;

		} else if (category == FOURS) {
			score = 0;
			for (int i = 0; i < N_DICE; i++) {
				if (currentDice[i] == 4)
					score = score + 4;

			}
			scoreboard[FOURS - 1][player - 1] = score;

		} else if (category == FIVES) {
			score = 0;
			for (int i = 0; i < N_DICE; i++) {
				if (currentDice[i] == 5)
					score = score + 5;

			}
			scoreboard[FIVES - 1][player - 1] = score;

		} else if (category == SIXES) {
			score = 0;
			for (int i = 0; i < N_DICE; i++) {
				if (currentDice[i] == 6)
					score = score + 6;

			}
			scoreboard[SIXES - 1][player - 1] = score;

		} else if (category == THREE_OF_A_KIND) {
			score = 0;
			for (int i = 0; i < N_DICE; i++) {
				int sameDiceCounter = 0;
				for (int j = 0; j < N_DICE; j++) {
					if (currentDice[i] == currentDice[j])
						sameDiceCounter++;
				}

				if (sameDiceCounter >= 3) {
					score = sumOfAllDices(i);

				}
			}
		} else if (category == FOUR_OF_A_KIND) {
			score = 0;
			for (int i = 0; i < N_DICE; i++) {
				int sameDiceCounter = 0;
				for (int j = 0; j < N_DICE; j++) {
					if (currentDice[i] == currentDice[j])
						sameDiceCounter++;
				}
				if (sameDiceCounter >= 4) {
					score = sumOfAllDices(i);

				}
			}
		} else if (category == FULL_HOUSE) {
			if (checkTheCategory(category)) {
				score = 25;

			} else {
				score = 0;
			}
		}

		else if (category == SMALL_STRAIGHT) {
			if (checkTheCategory(category)) {
				score = 30;

			} else {
				score = 0;
			}
		} else if (category == LARGE_STRAIGHT) {
			if (checkTheCategory(category)) {
				score = 40;
			} else {
				score = 0;
			}
		} else if (category == YAHTZEE) {
			if (checkTheCategory(category)) {
				score = 50;
			} else {
				score = 0;
			}

		} else if (category == CHANCE) {
			for (int dice = 0; dice < N_DICE; dice++) {
				score = score + currentDice[dice];
			}
		}
		display.updateScorecard(category, player, score);
		scoreboard[category - 1][player - 1] = score;
	}

	private boolean checkTheCategory(int category) {
		ArrayList<Integer> ones = new ArrayList<Integer>();
		for (int i = 0; i < N_DICE; i++) {
			if (currentDice[i] == 1)
				ones.add(1);
		}
		ArrayList<Integer> twos = new ArrayList<Integer>();
		for (int i = 0; i < N_DICE; i++) {
			if (currentDice[i] == 2)
				twos.add(1);
		}
		ArrayList<Integer> threes = new ArrayList<Integer>();
		for (int i = 0; i < N_DICE; i++) {
			if (currentDice[i] == 3)
				threes.add(1);
		}
		ArrayList<Integer> fours = new ArrayList<Integer>();
		for (int i = 0; i < N_DICE; i++) {
			if (currentDice[i] == 4)
				fours.add(1);
		}
		ArrayList<Integer> fives = new ArrayList<Integer>();
		for (int i = 0; i < N_DICE; i++) {
			if (currentDice[i] == 5)
				fives.add(1);
		}
		ArrayList<Integer> sixes = new ArrayList<Integer>();
		for (int i = 0; i < N_DICE; i++) {
			if (currentDice[i] == 6)
				sixes.add(1);
		}

		if (category == FULL_HOUSE) {
			if (ones.size() == 3 || twos.size() == 3 || threes.size() == 3 || fours.size() == 3 || fives.size() == 3
					|| sixes.size() == 3) {
				if (ones.size() == 2 || twos.size() == 2 || threes.size() == 2 || fours.size() == 2 || fives.size() == 2
						|| sixes.size() == 2) {
					return true;
				}
			}

		} else if (category == SMALL_STRAIGHT) {
			if (ones.size() >= 1 && twos.size() >= 1 && threes.size() >= 1 && fours.size() >= 1) {
				return true;
			} else if (twos.size() >= 1 && threes.size() >= 1 && fours.size() >= 1 && fives.size() >= 1) {
				return true;
			} else if (threes.size() >= 1 && fours.size() >= 1 && fives.size() >= 1 && sixes.size() >= 1) {
				return true;
			}

		} else if (category == LARGE_STRAIGHT) {
			if (ones.size() >= 1 && twos.size() >= 1 && threes.size() >= 1 && fours.size() >= 1 && fives.size() >= 1) {
				return true;
			} else if (twos.size() >= 1 && threes.size() >= 1 && fours.size() >= 1 && fives.size() >= 1
					&& sixes.size() >= 1) {
				return true;
			}
		}

		else if (category == YAHTZEE) {
			if (ones.size() == 5 || twos.size() == 5 || threes.size() == 5 || fours.size() == 5 || fives.size() == 5
					|| sixes.size() == 5) {
				return true;
			}
		} else if (category == CHANCE) {
			return true;
		}
		return false;
	}

	private boolean categoryIsValid(int player, int category) {
		if (scoreboard[category - 1][player - 1] != -1 || category == UPPER_SCORE || category == UPPER_BONUS
				|| category == LOWER_SCORE || category == TOTAL)
			return false;
		return true;
	}

	private int sumOfAllDices(int i) {
		return currentDice[0] + currentDice[1] + currentDice[2] + currentDice[3] + currentDice[4];
	}

	private void getTheResults() {
		int result = 0;
		for (int player = 0; player < nPlayers; player++) {
			for (int category = 1; category < SIXES; category++) {
				result = result + scoreboard[player][category];
			}
			scoreboard[UPPER_SCORE - 1][player] = result;
			display.updateScorecard(UPPER_SCORE, player + 1, result);
			getTheUpperBonus(player);
			result = 0;
			for (int category = 8; category < CHANCE; category++) {
				result = result + scoreboard[category][player];
			}
			scoreboard[LOWER_SCORE - 1][player] = result;
			display.updateScorecard(LOWER_SCORE, player + 1, result);
			scoreboard[TOTAL - 1][player] = scoreboard[UPPER_SCORE - 1][player] + scoreboard[UPPER_BONUS - 1][player]
					+ scoreboard[LOWER_SCORE - 1][player];
			display.updateScorecard(TOTAL, player + 1, scoreboard[TOTAL - 1][player]);
			result = 0;
		}

	}

	private void getTheUpperBonus(int player) {
		if (scoreboard[UPPER_SCORE - 1][player] > 63) {
			scoreboard[UPPER_BONUS - 1][player] = 35;
			display.updateScorecard(UPPER_BONUS, player + 1, 35);
		} else if (scoreboard[UPPER_SCORE - 1][player] < 63) {
			scoreboard[UPPER_BONUS - 1][player] = 0;
			display.updateScorecard(UPPER_BONUS, player + 1, 0);
		}
	}

	private void getTheWinner() {
		int[] results = new int[nPlayers];
		int score = 0;
		int winner = 0;
		for (int j = 0; j < nPlayers; j++) {
			if (scoreboard[TOTAL - 1][player] > score) {
				score = scoreboard[TOTAL - 1][player];
				winner = player;
			}
		}
		display.printMessage(
				"Congratulations" + playerNames[winner] + "you are the winner , with a total score of" + score);
	}

	/* Private instance variables */
	private int player;
	private int sameDiceCounter;
	private int nPlayers;
	private String[] playerNames;
	private YahtzeeDisplay display;
	private RandomGenerator rgen = new RandomGenerator();

}