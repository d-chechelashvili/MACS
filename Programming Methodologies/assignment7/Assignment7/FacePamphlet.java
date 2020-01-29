
/* 
 * File: FacePamphlet.java
 * -----------------------
 * When it is finished, this program will implement a basic social network
 * management system.
 */

import acm.program.*;
import acm.graphics.*;
import acm.util.*;
import java.awt.event.*;
import java.util.Iterator;

import javax.swing.*;

public class FacePamphlet extends Program implements FacePamphletConstants {
	/**
	 * This method has the responsibility for initializing the interactors in the
	 * application, and taking care of any other initialization that needs to be
	 * performed.
	 */
	public void init() {
		db = new FacePamphletDatabase();
		canv = new FacePamphletCanvas();
		add(canv);
		addInteractors();
		addActionListeners();
	}

	private void addInteractors() {
		addNorthSector();
		addWestSector();
	}

	// adds interactors to applet that belong in north sector
	private void addNorthSector() {
		add(new JLabel("Name"), NORTH);
		tf[0] = new JTextField(TEXT_FIELD_SIZE);
		add(tf[0], NORTH);
		tf[0].addActionListener(this);
		add(new JButton("Add"), NORTH);
		add(new JButton("Delete"), NORTH);
		add(new JButton("Lookup"), NORTH);
	}

	// adds interactors to applet that belong in west sector 4 textfields, 4 buttons
	// and 3 empty labels
	private void addWestSector() {
		for (int i = 0; i < 4; i++) {
			addTextField(i);
			addButton(i);
			addEmptyLabel(i);
		}
	}

	private void addTextField(int i) {
		tf[i + 1] = new JTextField(TEXT_FIELD_SIZE);
		add(tf[i + 1], WEST);
		tf[i + 1].addActionListener(this);
	}

	private void addButton(int i) {
		switch (i) {
		case 0:
			add(new JButton("Change Status"), WEST);
			break;
		case 1:
			add(new JButton("Change Picture"), WEST);
			break;
		case 2:
			add(new JButton("Add Friend"), WEST);
			break;
		case 3:
			add(new JButton("Unfriend"), WEST);
			break;
		}
	}

	private void addEmptyLabel(int i) {
		if (i < 3)
			add(new JLabel(EMPTY_LABEL_TEXT), WEST);
	}

	/**
	 * This class is responsible for detecting when the buttons are clicked or
	 * interactors are used, so you will have to add code to respond to these
	 * actions.
	 */
	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("Add")) {
			addProfile(0);
			tf[0].setText("");
		} else if (e.getActionCommand().equals("Delete")) {
			delete(0);
			tf[0].setText("");
		} else if (e.getActionCommand().equals("Lookup")) {
			lookUp(0);
			tf[0].setText("");
		} else if (e.getSource() == tf[1] || e.getActionCommand().equals("Change Status")) {
			changeStatus(1);
			tf[1].setText("");
		} else if (e.getSource() == tf[2] || e.getActionCommand().equals("Change Picture")) {
			changePicture(2);
			tf[2].setText("");
		} else if (e.getSource() == tf[3] || e.getActionCommand().equals("Add Friend")) {
			addFriend(3);
			tf[3].setText("");
		} else if (e.getSource() == tf[4] || e.getActionCommand().equals("Unfriend")) {
			unfriend(4);
			tf[4].setText("");
		}
	}

	// checks if textfield contains any string
	private boolean notEmpty(int textFieldIndex) {
		if (tf[textFieldIndex].getText().length() > 0) {
			return true;
		}
		return false;
	}

	// creates a new profile with given name and adds it to the database
	private void addProfile(int index) {
		if (notEmpty(index)) {
			if (db.containsProfile(tf[index].getText())) {
				current = db.getProfile(tf[index].getText());
				canv.displayProfile(current);
				canv.showMessage("A profile with the name " + current.getName() + " already exists");
			} else {
				current = new FacePamphletProfile(tf[index].getText());
				db.addProfile(current);
				canv.displayProfile(current);
				canv.showMessage("New profile created");
			}
		}
	}

	// deletes a profile with given name and updates all friendlists accordingly
	private void delete(int index) {
		if (notEmpty(index)) {
			current = null;
			canv.displayProfile(current);
			FacePamphletProfile profile = db.getProfile(tf[index].getText());
			if (profile != null) {
				String name = profile.getName();
				db.deleteProfile(name);
				canv.showMessage("Profile of " + name + " deleted");
			} else {
				canv.showMessage("A profile with the name " + tf[index].getText() + " does not exist");
			}
		}
	}

	// searchs for a profile with given name and if it exists displays it
	private void lookUp(int index) {
		if (notEmpty(index)) {
			current = db.getProfile(tf[index].getText());
			canv.displayProfile(current);
			if (current != null) {
				canv.showMessage("Displaying " + current.getName());
			} else {
				canv.showMessage("A profile with the name " + tf[index].getText() + " does not exist");
			}
		}
	}

	// changes status of current profile
	private void changeStatus(int index) {
		if (notEmpty(index)) {
			if (current != null) {
				current.setStatus(tf[index].getText());
				canv.displayProfile(current);
				canv.showMessage("Status Updated to " + tf[index].getText());
			} else {
				canv.showMessage("Please select a profile to change status");
			}
		}
	}

	// changes picture of current profile
	private void changePicture(int index) {
		if (notEmpty(index)) {
			if (current != null) {
				GImage image = null;
				try {
					image = new GImage(tf[index].getText());
					current.setImage(image);
					canv.displayProfile(current);
					canv.showMessage("Picture updated");
				} catch (ErrorException ex) {
					canv.showMessage("Unable to open image file: " + tf[index].getText());
				}
			} else {
				canv.showMessage("Please select a profile to change image");
			}
		}
	}

	// adds person as a friend if entered text satisfies all conditions
	private void addFriend(int index) {
		if (notEmpty(index)) {
			String friendName = tf[index].getText();
			if (checkCurrent(friendName) && valid(friendName) && !alreadyFriend(friendName)) {
				db.getProfile(friendName).addFriend(current.getName());
				canv.displayProfile(current);
				canv.showMessage(friendName + " added as a friend");
			}
		}
	}

	// checks if there is current profile to add a friend to and also controls that
	// people cannot be friends with themselves
	private boolean checkCurrent(String friendName) {
		if (current != null && !current.getName().equals(friendName)) {
			return true;
		}
		canv.showMessage("Please select a profile to add friend");
		return false;
	}

	// checks if a profile with given name exists
	private boolean valid(String friendName) {
		if (db.containsProfile(friendName))
			return true;
		canv.showMessage("Profile with that name does not exist");
		return false;
	}

	// checks if a person already has a person with entered name in friendlist
	private boolean alreadyFriend(String friendName) {
		if (!current.addFriend(friendName)) {
			canv.showMessage(friendName + " is already friend of " + current.getName());
			return true;
		}
		return false;
	}

	// removes the person with entered name from friendlist
	private void unfriend(int index) {
		if (notEmpty(index) && current != null) {
			String friend = tf[index].getText();
			if (current.removeFriend(friend)) {
				canv.displayProfile(current);
				canv.showMessage(friend + " has been removed from friendlist");
			} else {
				canv.showMessage(friend + " is not in the friendlist");
			}
		}
	}

	private JTextField[] tf = new JTextField[TEXTFIELD_NUM];
	private FacePamphletDatabase db;
	private FacePamphletProfile current;
	private FacePamphletCanvas canv;
}