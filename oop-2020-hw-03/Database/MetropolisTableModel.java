import javax.swing.table.AbstractTableModel;
import java.sql.SQLException;

public class MetropolisTableModel extends AbstractTableModel implements MetropolisConstants{
    private MetropolisDatabase db;

    public MetropolisTableModel(MetropolisDatabase db){
        this.db = db;
    }

    @Override
    public int getRowCount() {
        int rowCount = 0;
        try {
            rowCount = db.getRowCount();
        } catch (SQLException throwables) {
            throwables.printStackTrace();
        }
        return rowCount;
    }

    @Override
    public int getColumnCount() {
        return db.getColumnCount();
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
        Object value = null;
        try {
            value = db.getValueAt(rowIndex, columnIndex);
        } catch (SQLException throwables) {
            throwables.printStackTrace();
        }
        return value;
    }

    public void add(String metropolisName, String continentName, long populationCount) throws SQLException {
        db.add(metropolisName, continentName, populationCount);
        fireTableStructureChanged();
    }


    public void search(String metropolisName, String continentName, long populationCount,
                       int selectedIndex, int selectedIndex1) throws SQLException {
        db.search(metropolisName, continentName, populationCount, selectedIndex, selectedIndex1);
        fireTableStructureChanged();
    }
}
