import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.SQLException;

public class MetropolisFrame extends JFrame implements MetropolisConstants {
    private JTextField metropolisName;
    private JTextField continentName;
    private JTextField populationCount;
    private JButton addButton;
    private JButton searchButton;
    private JComboBox<String> populationRestriction;
    private JComboBox<String> matchType;
    private MetropolisTableModel model;

    public MetropolisFrame(MetropolisTableModel model) {
        super("Metropolis Viewer");
        this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        this.model = model;
        LayoutManager lm = new BorderLayout();
        setLayout(lm);
        fillNorthSector();
        fillEastSector();
        fillCentralSector();
        pack();
        setVisible(true);
    }


    private void fillNorthSector() {
        JPanel northPanel = new JPanel();
        add(northPanel, BorderLayout.NORTH);

        JLabel metropolis = new JLabel("Metropolis: ");
        northPanel.add(metropolis);

        metropolisName = new JTextField(TEXT_FIELD_SIZE);
        northPanel.add(metropolisName);

        JLabel continent = new JLabel("Continent: ");
        northPanel.add(continent);

        continentName = new JTextField(TEXT_FIELD_SIZE);
        northPanel.add(continentName);

        JLabel population = new JLabel("Population: ");
        northPanel.add(population);

        populationCount = new JTextField(TEXT_FIELD_SIZE);
        northPanel.add(populationCount);
    }

    private void fillEastSector() {
        JPanel eastPanel = new JPanel();
        add(eastPanel, BorderLayout.EAST);
        eastPanel.setLayout(new BoxLayout(eastPanel, BoxLayout.PAGE_AXIS));

        addButtons(eastPanel);
        addComboBoxes(eastPanel);
    }

    private void addButtons(JPanel eastPanel) {
        addButton = new JButton("Add");
        addButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    add();
                } catch (SQLException throwables) {
                    throwables.printStackTrace();
                }
            }
        });
        eastPanel.add(addButton);

        searchButton = new JButton("Search");
        searchButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    search();
                } catch (SQLException throwables) {
                    throwables.printStackTrace();
                }
            }
        });
        eastPanel.add(searchButton);
    }


    private void addComboBoxes(JPanel eastPanel) {
        Box searchOptions = new Box(BoxLayout.PAGE_AXIS);
        searchOptions.setBorder(new TitledBorder("Search Options"));

        populationRestriction = new JComboBox<>();
        populationRestriction.addItem("Population Larger Than");
        populationRestriction.addItem("Population Smaller Than");
        populationRestriction.setMaximumSize(
                new Dimension(Integer.MAX_VALUE, populationRestriction.getMinimumSize().height));

        matchType = new JComboBox<>();
        matchType.addItem("Exact Match");
        matchType.addItem("Partial Match");
        matchType.setMaximumSize(
                new Dimension(Integer.MAX_VALUE, matchType.getMinimumSize().height));

        searchOptions.add(populationRestriction);
        searchOptions.add(matchType);
        eastPanel.add(searchOptions);
    }

    private void fillCentralSector() {
        JTable metropolisTable = new JTable(model);

        JScrollPane metropolisTableScroll = new JScrollPane(metropolisTable,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

        add(metropolisTableScroll, BorderLayout.CENTER);
    }

    private void add() throws SQLException {
        long populationCount = -1;
        String populationCountStr = this.populationCount.getText();
        if(populationCountStr.length() > 0) populationCount = Long.parseLong(populationCountStr);

        model.add(metropolisName.getText(), continentName.getText(), populationCount);
    }

    private void search() throws SQLException {
        long populationCount = -1;
        String populationCountStr = this.populationCount.getText();
        if(populationCountStr.length() > 0) populationCount = Long.parseLong(populationCountStr);

        model.search(metropolisName.getText(), continentName.getText(), populationCount,
                populationRestriction.getSelectedIndex(), matchType.getSelectedIndex());
    }
}
