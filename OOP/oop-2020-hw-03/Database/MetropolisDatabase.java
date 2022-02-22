
import java.sql.*;


public class MetropolisDatabase implements MetropolisConstants {
    private Connection dbCon;
    private ResultSet lastResultSet;

    public MetropolisDatabase(String dbName) {
        try {
            Class.forName("com.mysql.cj.jdbc.Driver");
            dbCon = DriverManager.getConnection(
                    "jdbc:mysql://localhost",
                    "root",
                    "my_password");
            Statement useDbStm = dbCon.createStatement();
            useDbStm.execute("USE " + dbName);
        } catch (SQLException | ClassNotFoundException e) {
            throw new RuntimeException("");
        }
    }


    public int getRowCount() throws SQLException {
        int lastRowIndex = 0;
        if (lastResultSet == null) return 0;
        int currRowIndex = lastResultSet.getRow();
        lastResultSet.last();
        lastRowIndex = lastResultSet.getRow();
        lastResultSet.absolute(currRowIndex);
        return lastRowIndex;
    }

    public int getColumnCount() {
        if (lastResultSet == null) return 0;
        return DATABASE_COLUMN_COUNT;
    }

    public Object getValueAt(int rowIndex, int columnIndex) throws SQLException {
        Object value = null;
        int currRowIndex = lastResultSet.getRow();
        lastResultSet.absolute(rowIndex + 1);
        value = lastResultSet.getObject(columnIndex + 1);
        return value;
    }

    public void add(String metropolisName, String continentName, long populationCount) throws SQLException {

        if (populationCount < 0 || metropolisName.equals("") || continentName.equals("")) return;

        Statement addStatement = dbCon.createStatement();
        addStatement.executeUpdate("INSERT INTO metropolises VALUES(" +
                "\"" + metropolisName + "\"" + "," + "\"" + continentName + "\"" + "," +
                "\"" + populationCount + "\"" + ");");
        search(metropolisName, continentName,
                populationCount, POPULATION_PULLDOWN_GREATER_INDEX, MATCHTYPE_PULLDOWN_EXACT_INDEX);

    }


    public void search(String metropolisName, String continentName, long populationCount,
                       int populatuionIndex, int matchTypeIndex) throws SQLException {

        Statement searchStatement = dbCon.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE, ResultSet.CONCUR_READ_ONLY);
        StringBuilder query = new StringBuilder("SELECT * FROM metropolises");
        addWhereClauses(query, metropolisName, continentName, populationCount, populatuionIndex, matchTypeIndex);
        query.append(";");
        lastResultSet = searchStatement.executeQuery(query.toString());

    }

    private void addWhereClauses(StringBuilder query, String metropolisName, String continentName, long populationCount,
                                 int populatuionIndex, int matchTypeIndex) {

        if (metropolisName.length() > 0 && continentName.length() > 0 && populationCount < 0) return;

        boolean whereUsed = false;
        whereUsed = restrictPopulationColumn(whereUsed, populatuionIndex, populationCount, query);
        whereUsed = restrictStringColumns(whereUsed, matchTypeIndex, "metropolis", metropolisName, query);
        restrictStringColumns(whereUsed, matchTypeIndex, "continent", continentName, query);
    }

    private boolean restrictPopulationColumn(boolean whereUsed, int populatuionIndex,
                                             long populationCount, StringBuilder query) {
        if (populationCount >= 0) {
            query.append(" WHERE population ");
            whereUsed = true;
            if (populatuionIndex == POPULATION_PULLDOWN_GREATER_INDEX) query.append(">= ");
            else query.append("< ");
            query.append(populationCount);
        }
        return whereUsed;
    }

    private boolean restrictStringColumns(boolean whereUsed, int matchTypeIndex,
                                          String columnName, String variableName, StringBuilder query) {
        if (variableName.length() > 0) {
            if (!whereUsed) {
                query.append(" WHERE ");
                whereUsed = true;
            } else query.append(" AND ");
            query.append(columnName).append(" ");
            if (matchTypeIndex == MATCHTYPE_PULLDOWN_EXACT_INDEX)
                query.append("= \"").append(variableName).append("\"");
            else
                query.append("LIKE \"%").append(variableName).append("%\"");
        }
        return whereUsed;
    }


//    private void print(){
//        try {
//            while (lastResultSet.next()) {
//                String name = lastResultSet.getString("metropolis");
//                String continent = lastResultSet.getString("continent");
//                long pop = lastResultSet.getLong("population");
//                System.out.println(name + "\t" + continent + "\t" + pop);
//            }
//        }catch (SQLException throwables) {
//            throwables.printStackTrace();
//        }
//        System.out.println(getRowCount());
//        System.out.println(getColumnCount());
//    }
}
