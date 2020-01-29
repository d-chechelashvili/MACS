
/*
 * File: NameSurfer.java
 * ---------------------
 * When it is finished, this program will implements the viewer for
 * the baby-name database described in the assignment handout.
 */

import acm.program.*;
import java.awt.event.*;
import javax.swing.*;

public class NameSurfer extends Program implements NameSurferConstants {

	/* Method: init() */
	/**
	 * This method has the responsibility for reading in the data base and
	 * initializing the interactors at the bottom of the window.
	 */
	public void init() {
		db = new NameSurferDataBase(NAMES_DATA_FILE);
		graph = new NameSurferGraph();
		add(graph);
		addInteractors();
		addActionListeners();
		nameField.addActionListener(this);
	}

	// adds interactors on the applet
	private void addInteractors() {
		add(new JLabel("Name"), SOUTH);
		add(nameField, SOUTH);
		add(new JButton("Graph"), SOUTH);
		add(graphChooser, SOUTH);
		add(new JButton("Clear selected"), SOUTH);
		add(new JButton("Clear all"), SOUTH);
	}

	/* Method: actionPerformed(e) */
	/**
	 * This class is responsible for detecting when the buttons are clicked, so you
	 * will have to define a method to respond to button actions.
	 */
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == nameField || e.getActionCommand().equals("Graph")) {
			graph();
		} else if (e.getActionCommand().equals("Clear selected") && graphChooser.getSelectedItem() != null) {
			removeSelected();
		} else if (e.getActionCommand().equals("Clear all")) {
			graph.clear();
			graph.update();
		}
	}

	// draws the graph for entered name if it is in the database
	private void graph() {
		NameSurferEntry entry = db.findEntry(nameField.getText());
		if (entry != null) {
			graph.addEntry(entry);
			graphChooser.addItem(entry.getName());
			graph.update();
		}
		nameField.setText("");
	}

	// removes graph for selected name
	private void removeSelected() {
		int index = graphChooser.getSelectedIndex();
		graph.clearSelected(index);
		graphChooser.removeItemAt(index);
		graphChooser.removeAll();
		graph.update();
	}

	private JTextField nameField = new JTextField(10);
	private NameSurferDataBase db;
	private NameSurferGraph graph;
	private JComboBox<String> graphChooser = new JComboBox<String>();
}