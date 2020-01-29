
/*
 * File: NameSurferGraph.java
 * ---------------------------
 * This class represents the canvas on which the graph of
 * names is drawn. This class is responsible for updating
 * (redrawing) the graphs whenever the list of entries changes or the window is resized.
 */

import acm.graphics.*;

import java.awt.Color;
import java.awt.Font;
import java.awt.event.*;
import java.util.*;

public class NameSurferGraph extends GCanvas implements NameSurferConstants, ComponentListener {

	/**
	 * Creates a new NameSurferGraph object that displays the data.
	 */
	public NameSurferGraph() {
		addComponentListener(this);
	}

	/**
	 * Clears the list of name surfer entries stored inside this class.
	 */
	public void clear() {
		entryList.clear();
		colorList.clear();
		cl.clear();
	}

	// clears data about the entry user chose to clear
	public void clearSelected(int index) {
		entryList.remove(index);
		cl.clearSelected(colorList.get(index));
		colorList.remove(index);
	}

	/* Method: addEntry(entry) */
	/**
	 * Adds a new NameSurferEntry to the list of entries on the display. Note that
	 * this method does not actually draw the graph, but simply stores the entry;
	 * the graph is drawn by calling update.
	 */
	public void addEntry(NameSurferEntry entry) {
		if (entry != null && !entryList.contains(entry)) {
			entryList.add(entry);
			colorList.add(cl.getAvailableColor());
		}
	}

	/**
	 * Updates the display image by deleting all the graphical objects from the
	 * canvas and then reassembling the display according to the list of entries.
	 * Your application must call update after calling either clear or addEntry;
	 * update is also called whenever the size of the canvas changes.
	 */
	public void update() {
		removeAll();
		drawBorders();
		drawGraphs();
	}

	// draws vertical borders for each decade and top and bottom boundaries of
	// graphs
	private void drawBorders() {
		drawVertical();
		drawHorizontal();
	}

	// draws top and bottom boundaries of graphs
	private void drawVertical() {
		GLine topLine = new GLine(0, GRAPH_MARGIN_SIZE, getWidth(), GRAPH_MARGIN_SIZE);
		GLine bottomLine = new GLine(0, getHeight() - GRAPH_MARGIN_SIZE, getWidth(), getHeight() - GRAPH_MARGIN_SIZE);
		add(topLine);
		add(bottomLine);
	}

	// draws vertical borders for each decade and adds numbers showing which column
	// belongs to which decade
	private void drawHorizontal() {
		double x = 0;
		for (int i = 0; i < NDECADES; i++) {
			GLine line = new GLine(x, 0, x, getHeight());
			add(line);
			drawDecade(x, i);
			x += getWidth() / NDECADES;
		}
	}

	// draws decades in their respective columns
	private void drawDecade(double x, int i) {
		GLabel decade = new GLabel("" + (START_DECADE + i * 10));
		decade.setLocation(x, getHeight());
		decade.setFont(getSizeableFont());
		add(decade);
	}

	// gets the font that will be resized as the applet is resized
	private Font getSizeableFont() {
		int x = (APPLICATION_WIDTH + APPLICATION_HEIGHT) / 18;
		int fontSize = (getHeight() + getWidth()) / x;
		Font f = new Font("TimesRoman", Font.PLAIN, fontSize);
		return f;
	}

	// draws graph for each entry
	private void drawGraphs() {
		for (int i = 0; i < entryList.size(); i++) {
			drawGraph(i);
		}
	}

	// draws graph for single entry
	private void drawGraph(int entry) {
		double x0 = 0;
		for (int i = 0; i < NDECADES - 1; i++) {
			double x = x0 + getWidth() / NDECADES;
			GLine line = new GLine(x0, getY(entry, i), x, getY(entry, i + 1));
			line.setColor(colorList.get(entry));
			add(line);
			addName(entry, i, line);
			x0 = x;
		}
	}

	// calculates the y coordinate according of given name's placement in given
	// decade
	private double getY(int entryNum, int decade) {
		double y = 0;
		NameSurferEntry entry = entryList.get(entryNum);
		if (entry.getRank(decade) == 0) {
			y = getHeight() - GRAPH_MARGIN_SIZE;
		} else {
			System.out.println(entry);
			y = GRAPH_MARGIN_SIZE
					+ ((((double) getHeight() - 2 * GRAPH_MARGIN_SIZE) / MAX_RANK) * entry.getRank(decade));
		}
		return y;
	}

	// adds names as labels for each decade where graph crosses it's border
	// also there is a special case on last decade because graph doesn't cross it
	private void addName(int entryNum, int decade, GLine line) {
		NameSurferEntry entry = entryList.get(entryNum);
		int rank = entry.getRank(decade);
		GLabel name;
		if (rank == 0) {
			name = new GLabel(entry.getName() + "*");
		} else {
			name = new GLabel(entry.getName() + " " + rank);
		}
		name.setLocation(line.getStartPoint().getX(), line.getStartPoint().getY());
		name.setFont(getSizeableFont());
		name.setColor(line.getColor());
		add(name);

		if (decade == NDECADES - 2) {
			addLastLabel(entry, decade, line);
		}
	}

	// adds the name as a label for last decade
	private void addLastLabel(NameSurferEntry entry, int decade, GLine line) {
		int lastRank = entry.getRank(decade + 1);
		GLabel lastLabel;
		if (lastRank == 0) {
			lastLabel = new GLabel(entry.getName() + "*");
		} else {
			lastLabel = new GLabel(entry.getName() + " " + lastRank);
		}
		lastLabel.setLocation(line.getEndPoint().getX(), line.getEndPoint().getY());
		lastLabel.setFont(getSizeableFont());
		lastLabel.setColor(line.getColor());
		add(lastLabel);
	}

	/* Implementation of the ComponentListener interface */
	public void componentHidden(ComponentEvent e) {
	}

	public void componentMoved(ComponentEvent e) {
	}

	public void componentResized(ComponentEvent e) {
		update();
	}

	public void componentShown(ComponentEvent e) {
	}

	private ArrayList<NameSurferEntry> entryList = new ArrayList<>();
	private ArrayList<Color> colorList = new ArrayList<Color>();
	private NameSurferColors cl = new NameSurferColors();
}