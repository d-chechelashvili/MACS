/*
 * File: CheckerboardKarel.java
 * ----------------------------
 * When you finish writing it, the CheckerboardKarel class should draw
 * a checkerboard using beepers, as described in Assignment 1.  You
 * should make sure that your program works for all of the sample
 * worlds supplied in the starter folder.
 */

import stanford.karel.*;


public class CheckerboardKarel extends SuperKarel  {
	
	public void run() {
		repeatFillOneRow();
		
	}
	private void repeatFillOneRow() {
		
		if(frontIsClear()){
			
			while(frontIsClear()) {
				fillOneRow();
				changeCourse(); //to set up another fillOneRow
				goBackToCheck(); //to make it work for every kind of map
		}
			 
		 
		
		}else {
			//it means karel's world has only one Column
			turnLeft();
			fillOneRow();
		}
	}
	private void fillOneRow() {
	
	//precondition - (1/Xmax , y) >/<
	//postcondition - (Xmax/1 , y) >/<  
		
		putBeeper();
		while(frontIsClear()) {
			move();
			  //to prevent OBOB
			if(frontIsClear()) {
				move();
				putBeeper();
				
			}
			
		}
	}
	private void changeCourse() {
	
	//pre - (Xmax/1 ,y) >/<
	//post - (Xmax/1, y+1) </>
		
		if (facingEast()) {
			turnLeft();
			//to prevent OBOB
			if(frontIsClear()) {
				move();
				turnLeft();	
			}
			
		} else  {
			turnRight();
			//to prevent OBOB
			if(frontIsClear()) {
				move();
				turnRight();
			}
			
		}
		
	}
	private void goBackToCheck() {
		
		//beepers should be surrounded only with empty squares
		
			if (facingWest()) {
			turnLeft();
			move();
			turnLeft();
			if(beepersPresent()) {
				changeCourse();
				move();
			}else {
				changeCourse();
			}
			
		}
		
	}
		
}

