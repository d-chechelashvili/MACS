
/*
 * File: Breakout.java
 * -------------------
 * Name:
 * Section Leader:
 * 
 * This file will eventually implement the game of Breakout.
 */

import acm.graphics.*;
import acm.program.*;
import acm.util.*;

import java.applet.*;
import java.awt.*;
import java.awt.event.*;

public class Breakout extends GraphicsProgram {

	/** Width and height of application window in pixels */
	public static final int APPLICATION_WIDTH = 400;
	public static final int APPLICATION_HEIGHT = 600;

	/** Dimensions of game board (usually the same) */
	private static final int WIDTH = APPLICATION_WIDTH;
	private static final int HEIGHT = APPLICATION_HEIGHT;

	/** Dimensions of the paddle */
	private static final int PADDLE_WIDTH = 60;
	private static final int PADDLE_HEIGHT = 10;

	/** Offset of the paddle up from the bottom */
	private static final int PADDLE_Y_OFFSET = 30;

	/** Number of bricks per row */
	private static final int NBRICKS_PER_ROW = 10;

	/** Number of rows of bricks */
	private static final int NBRICK_ROWS = 10;

	/** Separation between bricks */
	private static final int BRICK_SEP = 4;

	/** Width of a brick */
	private static final int BRICK_WIDTH = (WIDTH - (NBRICKS_PER_ROW - 1) * BRICK_SEP) / NBRICKS_PER_ROW;

	/** Height of a brick */
	private static final int BRICK_HEIGHT = 8;

	/** Radius of the ball in pixels */
	private static final int BALL_RADIUS = 10;

	/** Offset of the top brick row from the top */
	private static final int BRICK_Y_OFFSET = 70;

	/** Number of turns */
	private static final int NTURNS = 3;

	private static final int DELAY = 12;

	private static final int BALL_DIAM = 2 * BALL_RADIUS;

	private RandomGenerator rgen = RandomGenerator.getInstance();
	private GRect paddle;
	private GOval ball;
	private double vx; // horizontal speed of the ball
	private double vy; // vertical speed of the ball
	private int lives = 0;
	private int brickNum = NBRICKS_PER_ROW * NBRICK_ROWS;

	public void init() {
		addMouseListeners();
		setUpBoard();
	}

	public void run() {
		playGame();
	}

	private void setUpBoard() {
		setUpBricks();
		setUpBall();
		setUpPaddle();
	}

	// this method draws bricks BRICK_Y_OFFSET distance away from the top
	private void setUpBricks() {
		// if x0 starts from 0 there is BRICK_SEP length of space leftover.
		// to centre the diagram horizontally :
		int x0 = (WIDTH - (NBRICKS_PER_ROW * BRICK_WIDTH + (NBRICKS_PER_ROW - 1) * BRICK_SEP)) / 2;
		int y0 = BRICK_Y_OFFSET;
		// this for cycle resembles how many rows of bricks are there
		for (int i = 1; i <= NBRICK_ROWS; i++) {
			drawTier(x0, y0, i);
			y0 += BRICK_HEIGHT + BRICK_SEP; // each row of bricks is distanced from neighbours with BRICK_SEP
		}
	}

	// this method draws one row of bricks and sets their colour
	private void drawTier(int x0, int y0, int i) {
		for (int j = NBRICKS_PER_ROW; j > 0; j--) {
			GRect brick = new GRect(x0, y0, BRICK_WIDTH, BRICK_HEIGHT);
			brick.setFilled(true);
			// colours change every two row
			if (i == 1 || i == 2) {
				brick.setColor(Color.RED);
			} else if (i == 3 || i == 4) {
				brick.setColor(Color.ORANGE);
			} else if (i == 5 || i == 6) {
				brick.setColor(Color.YELLOW);
			} else if (i == 7 || i == 8) {
				brick.setColor(Color.GREEN);
			} else if (i == 9 || i == 10) {
				brick.setColor(Color.CYAN);
			}
			add(brick);
			x0 += BRICK_WIDTH + BRICK_SEP;
		}
	}

	// this method creates and adds paddle to the given location
	private void setUpPaddle() {
		double padX0 = WIDTH / 2 - PADDLE_WIDTH / 2;
		double padY0 = HEIGHT - (2 * PADDLE_Y_OFFSET + PADDLE_HEIGHT);
		paddle = new GRect(padX0, padY0, PADDLE_WIDTH, PADDLE_HEIGHT);
		paddle.setFilled(true);
		paddle.setColor(Color.BLACK);
		add(paddle);
	}

	// controls movement of the paddle(centre of the paddle must follow the mouse)
	public void mouseMoved(MouseEvent e) {
		if ((e.getX() >= PADDLE_WIDTH / 2) && (e.getX() <= WIDTH - PADDLE_WIDTH / 2)) {
			paddle.move(e.getX() - (paddle.getX() + PADDLE_WIDTH / 2), 0);
		}
	}

	// creates the ball and adds it in the middle
	private void setUpBall() {
		double ballX0 = WIDTH / 2 - BALL_RADIUS;
		double ballY0 = HEIGHT / 2 - BALL_RADIUS;
		ball = new GOval(ballX0, ballY0, BALL_DIAM, BALL_DIAM);
		ball.setFilled(true);
		ball.setColor(Color.BLACK);
		add(ball);
	}

	// controls the gameplay(rules, ball movement e.t.c)
	private void playGame() {
		waitForClick();
		moveBall();
	}

	// controls the ball movement and collision
	private void moveBall() {
		takeOff();
		while (lives < NTURNS && brickNum > 0) {
			pause(DELAY);
			ball.move(vx, vy);
			setCollisionRules();
		}
		remove(ball);
		displayResult();
	}

	// this method ball speed and movement angle during take off
	private void takeOff() {
		vy = 3.0;
		// ball should start moving with random angle , so vx is random
		vx = rgen.nextDouble(1.0, 3.0);
		if (rgen.nextBoolean(0.5))
			vx = -vx;
	}

	// controls the collision aspect of ball movement
	// every time ball collides it must change directions
	private void setCollisionRules() {
		setBounds();
		checkCollision();
	}

	// keeps the ball on the screen
	private void setBounds() {
		if (ball.getY() <= 0)
			vy = -vy;
		if (ball.getX() <= 0 || ball.getX() >= WIDTH - BALL_DIAM)
			vx = -vx;
		checkFall();
	}

	// checks if ball has fallen below the bottom wall
	private void checkFall() {
		if (ball.getY() > HEIGHT - BALL_DIAM) {
			resetBall();
			lives++;
			if (lives < NTURNS) {
				waitForClick();
			}
		}
	}

	// returns the ball to starting position
	private void resetBall() {
		remove(ball);
		ball.setLocation(WIDTH / 2 - BALL_RADIUS, HEIGHT / 2 - BALL_RADIUS);
		add(ball);
	}

	// defines what happens when ball collides with something other than walls
	private void checkCollision() {
		GObject collider = getCollidingObject();
		// paddle should just change ball's direction
		if (collider == paddle) {
			// when ball hits the paddle from top it should set the ball on top and change
			// its vertical speed
			if (ball.getY() + BALL_DIAM < paddle.getY() + Math.abs(vy)) {
				ball.setLocation(ball.getX(), paddle.getY() - BALL_DIAM);
				vy = -Math.abs(vy);
			}
			// else it means that ball has fallen below the point from where it would go upwards
			// it must fall
			else {
				vy = Math.abs(vy);
				checkSide();
			}
		}
		// if collider isn't paddle it will be a brick.
		// in that case ball should change direction and remove that brick
		else if (collider != null) {
			vy = -vy;
			remove(collider);
			brickNum--;
		}
	}

	// checks if ball has hit the paddle from left or right side and acts
	// accordingly
	private void checkSide() {
		// if ball collided from right
		if (Math.abs(ball.getX() - paddle.getX()) < Math
				.abs((ball.getX() + BALL_DIAM) - (paddle.getX() + PADDLE_WIDTH))) {
			ball.setLocation(Math.max(paddle.getX() - BALL_DIAM, 0), ball.getY());
			vx = -Math.abs(vx);
		}
		// if ball collided from left
		else {
			ball.setLocation(Math.min(paddle.getX() + PADDLE_WIDTH, WIDTH - BALL_DIAM), ball.getY());
			vx = Math.abs(vx);
		}

	}

	// this method checks edges of the rectangle in which ball is drawn
	// and if it collides with something returns that object
	private GObject getCollidingObject() {
		GObject leftTop = getElementAt(ball.getX(), ball.getY());
		GObject leftBottom = getElementAt(ball.getX(), ball.getY() + BALL_DIAM);
		GObject rightTop = getElementAt(ball.getX() + BALL_DIAM, ball.getY());
		GObject rightBottom = getElementAt(ball.getX() + BALL_DIAM, ball.getY() + BALL_DIAM);
		if (leftTop != null) {
			return leftTop;
		} else if (leftBottom != null) {
			return leftBottom;
		} else if (rightTop != null) {
			return rightTop;
		} else if (rightBottom != null) {
			return rightBottom;
		}
		return null;
	}

	// displays whether player lost or won
	private void displayResult() {
		removeAll();
		GLabel result;
		if (lives == NTURNS) {
			result = new GLabel("DEFEAT");
			result.setFont("LONDON-20");
			result.setColor(Color.RED);
			result.setLocation(WIDTH / 2 - result.getWidth() / 2, HEIGHT / 2 + result.getAscent() / 2);
			add(result);
		} else if (brickNum == 0) {
			result = new GLabel("VICTORY");
			result.setFont("LONDON-20");
			result.setColor(Color.BLUE);
			result.setLocation(WIDTH / 2 - result.getWidth() / 2, HEIGHT / 2 + result.getAscent() / 2);
			add(result);
		}
	}
}