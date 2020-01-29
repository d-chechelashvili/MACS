
/*
 * File: FacePamphletCanvas.java
 * -----------------------------
 * This class represents the canvas on which the profiles in the social
 * network are displayed.  NOTE: This class does NOT need to update the
 * display when the window is resized.
 */

import acm.graphics.*;
import java.awt.*;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.util.*;

public class FacePamphletCanvas extends GCanvas implements FacePamphletConstants {

	/**
	 * Constructor This method takes care of any initialization needed for the
	 * display
	 */
	public FacePamphletCanvas() {
		message = new GLabel("");
	}

	/**
	 * This method displays a message string near the bottom of the canvas. Every
	 * time this method is called, the previously displayed message (if any) is
	 * replaced by the new message text passed in.
	 */
	public void showMessage(String msg) {
		message.setLabel(msg);
		message.setFont(MESSAGE_FONT);
		double x = getWidth() / 2 - message.getWidth() / 2;
		double y = getHeight() - BOTTOM_MESSAGE_MARGIN;
		add(message, x, y);
	}

	/**
	 * This method displays the given profile on the canvas. The canvas is first
	 * cleared of all existing items (including messages displayed near the bottom
	 * of the screen) and then the given profile is displayed. The profile display
	 * includes the name of the user from the profile, the corresponding image (or
	 * an indication that an image does not exist), the status of the user, and a
	 * list of the user's friends in the social network.
	 */
	public void displayProfile(FacePamphletProfile profile) {
		removeAll();
		if (profile != null) {
			addName(profile);
			addImage(profile);
			addStatus(profile);
			addFriends(profile);
		}
	}

	// adds name to canvas
	private void addName(FacePamphletProfile profile) {
		GLabel name = new GLabel(profile.getName());
		name.setFont(PROFILE_NAME_FONT);
		nameBaseLine = name.getHeight() + TOP_MARGIN;
		add(name, LEFT_MARGIN, nameBaseLine);
	}

	// adds image to canvas
	private void addImage(FacePamphletProfile profile) {
		if (profile.getImage() != null) {
			GImage im = profile.getImage();
			//double sx = IMAGE_WIDTH / im.getWidth();
		//	double sy = IMAGE_HEIGHT / im.getHeight();
			//im.scale(sx, sy);
			im.setBounds(LEFT_MARGIN, TOP_MARGIN + nameBaseLine + IMAGE_MARGIN, IMAGE_WIDTH, IMAGE_HEIGHT);
			add(im, LEFT_MARGIN, TOP_MARGIN + nameBaseLine + IMAGE_MARGIN);
		} else {
			addEmptyIm();
		}
	}

	// if there is no image associated with the current profile adds sign that there
	// is no image there , where image should be
	private void addEmptyIm() {
		GRect emptyIm = new GRect(IMAGE_WIDTH, IMAGE_HEIGHT);
		add(emptyIm, LEFT_MARGIN, TOP_MARGIN + nameBaseLine + IMAGE_MARGIN);
		GLabel noImage = new GLabel("No Image");
		noImage.setFont(PROFILE_IMAGE_FONT);
		double labelX = emptyIm.getX() + emptyIm.getWidth() / 2 - noImage.getWidth() / 2;
		double labelY = emptyIm.getY() + emptyIm.getHeight() / 2 + noImage.getAscent() / 2;
		add(noImage, labelX, labelY);
	}

	// adds status to canvas
	private void addStatus(FacePamphletProfile profile) {
		GLabel status = new GLabel("No current status");
		status.setFont(PROFILE_STATUS_FONT);
		if (profile.getStatus() != "") {
			status.setLabel(profile.getName() + " is " + profile.getStatus());
		}
		double statusY = TOP_MARGIN + nameBaseLine + IMAGE_MARGIN + IMAGE_HEIGHT + STATUS_MARGIN;
		add(status, LEFT_MARGIN, statusY);
	}

	// adds friendlist to canvas
	private void addFriends(FacePamphletProfile profile) {
		double Y = TOP_MARGIN + IMAGE_MARGIN;
		addHeader(Y);
		Iterator<String> it = profile.getFriends();
		for (int i = 1; it.hasNext(); i++) {
			GLabel friend = new GLabel(it.next());
			friend.setFont(PROFILE_FRIEND_FONT);
			friend.setLocation(getWidth() / 2, Y + i * friend.getHeight());
			add(friend);
		}
	}

	// adds header of friendlist - label "Friends:" to canvas
	private void addHeader(double Y) {
		GLabel header = new GLabel("Friends:");
		header.setFont(PROFILE_FRIEND_LABEL_FONT);
		add(header, getWidth() / 2, Y);
	}

	private GLabel message;
	private double nameBaseLine;
}