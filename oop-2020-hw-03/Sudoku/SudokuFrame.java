import javafx.scene.layout.HBox;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;
import javax.swing.text.Document;

import java.awt.*;
import java.awt.event.*;


public class SudokuFrame extends JFrame {

    private JTextArea puzzle;
    private JTextArea solution;
    private JButton check;
    private JCheckBox autoChecker;

    public SudokuFrame() {
        super("Sudoku Solver");
        LayoutManager lm = new BorderLayout(4, 4);
        setLayout(lm);

        puzzle = new JTextArea(15, 20);
        add(puzzle, BorderLayout.CENTER);

        solution = new JTextArea(15, 20);
        add(solution, BorderLayout.EAST);

        puzzle.setBorder(new TitledBorder("Puzzle"));
        solution.setBorder(new TitledBorder("Solution"));

        Box controlsHolder = new Box(BoxLayout.X_AXIS);
        add(controlsHolder, BorderLayout.SOUTH);

        check = new JButton("Check");
        check.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                check();
            }
        });
        controlsHolder.add(check);

        autoChecker = new JCheckBox("Auto Check", true);
        controlsHolder.add(autoChecker);
        autoChecker();

        // Could do this:
        //setLocationByPlatform(true);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        pack();
        setVisible(true);
    }

    private void autoChecker() {
        Document doc = puzzle.getDocument();
        doc.addDocumentListener(new DocumentListener() {
            @Override
            public void insertUpdate(DocumentEvent e) {
                if (autoChecker.isSelected()) {
                    check();
                }
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                if (autoChecker.isSelected()) {
                    check();
                }
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                //Plain text components do not fire these events
            }
        });
    }

    private void check() {
        try {
            Sudoku current = new Sudoku(puzzle.getText());
            int solutionsNum = current.solve();
            String solutionText = current.getSolutionText();
            long timeElapsed = current.getElapsed();
            solution.setText(solutionText + "\n" + "solutions:" + solutionsNum + "\n" +
                    "elapsed:" + timeElapsed + "\n");
        } catch (RuntimeException ex) {
            solution.setText("Parsing Problem");
        }
    }


    public static void main(String[] args) {
        // GUI Look And Feel
        // Do this incantation at the start of main() to tell Swing
        // to use the GUI LookAndFeel of the native platform. It's ok
        // to ignore the exception.
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception ignored) {
        }

        SudokuFrame frame = new SudokuFrame();
    }

}
