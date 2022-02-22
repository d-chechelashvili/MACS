

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInstance;

import java.sql.*;

import static org.junit.jupiter.api.Assertions.*;

@TestInstance(TestInstance.Lifecycle.PER_CLASS)
class MetropolisDatabaseTest implements MetropolisConstants{

    private Connection conn;
    private MetropolisDatabase db;

    @BeforeAll
    void setUp() throws ClassNotFoundException, SQLException {
        Class.forName("com.mysql.cj.jdbc.Driver");
        conn = DriverManager.getConnection(
                "jdbc:mysql://localhost",
                "root",
                "my_password");
        Statement stm = conn.createStatement();
        stm.execute("CREATE DATABASE metropolis_test_db");
        stm.execute("USE metropolis_test_db");
        stm.execute(
                "CREATE TABLE metropolises (" +
                        "metropolis CHAR(64)," +
                        "continent  CHAR(64)," +
                        "population BIGINT" +
                        ");");

        String metropolisName = "Tbilisi";
        String continentName = "Europe";
        long populationCount = 2000000;
        stm.executeUpdate("INSERT INTO metropolises VALUES(" +
                "\"" + metropolisName + "\"" + "," + "\"" + continentName + "\"" +"," +
                "\"" + populationCount + "\"" +");");

        metropolisName = "Prague";
        continentName = "Europe";
        populationCount = 1000000;
        stm.executeUpdate("INSERT INTO metropolises VALUES(" +
                "\"" + metropolisName + "\"" + "," + "\"" + continentName + "\"" +"," +
                "\"" + populationCount + "\"" +");");

        metropolisName = "Bangkok";
        continentName = "Asia";
        populationCount = 1200000;
        stm.executeUpdate("INSERT INTO metropolises VALUES(" +
                "\"" + metropolisName + "\"" + "," + "\"" + continentName + "\"" +"," +
                "\"" + populationCount + "\"" +");");


        db = new MetropolisDatabase("metropolis_test_db");
        assertEquals(0, db.getRowCount());
        assertEquals(0, db.getColumnCount());
    }


    @Test
    void testConstructor(){
        assertThrows(RuntimeException.class, ()->{
            MetropolisDatabase invalidDb = new MetropolisDatabase("INVALID_DB_NAME");
        });

    }

    @Test
    void testSearch1() throws SQLException {
        String metropolisName = "Prague";
        String continentName = "Europe";
        long populationCount = 1000000;
        db.search(metropolisName, continentName, populationCount,
                POPULATION_PULLDOWN_GREATER_INDEX, MATCHTYPE_PULLDOWN_EXACT_INDEX);

        assertEquals(3, db.getColumnCount());
        assertEquals(1, db.getRowCount());
        assertEquals(metropolisName, db.getValueAt(0, 0));
        assertEquals(continentName, db.getValueAt(0, 1));
        assertEquals(populationCount, db.getValueAt(0, 2));


    }

    @Test
    void testSearch2() throws SQLException {
        db.search("","", 3000000,
                POPULATION_PULLDOWN_LESSER_INDEX, MATCHTYPE_PULLDOWN_PARTIAL_INDEX);

        assertEquals(3, db.getColumnCount());
        assertEquals(3, db.getRowCount());

        db.search("Tbilisi","", -1,
                POPULATION_PULLDOWN_LESSER_INDEX, MATCHTYPE_PULLDOWN_PARTIAL_INDEX);

        assertEquals(1, db.getRowCount());
        db.search("","asia", -1,
                POPULATION_PULLDOWN_LESSER_INDEX, MATCHTYPE_PULLDOWN_PARTIAL_INDEX);

        assertEquals(1, db.getRowCount());
    }

    @Test
    void testAdd() throws SQLException {
        db.add("New York", "North America", 13232131);
        assertEquals(1, db.getRowCount());
        assertEquals(3, db.getColumnCount());
        db.add("New York", "", 13232131);
        db.add("New York", "", -1);
        db.add("", "", -1);
        db.add("", "North America", -1);
        assertEquals(1, db.getRowCount());
        assertEquals(3, db.getColumnCount());


    }

    @Test
    void testGetValueAt() throws SQLException {
        db.search("","", 3000000,
                POPULATION_PULLDOWN_LESSER_INDEX, MATCHTYPE_PULLDOWN_PARTIAL_INDEX);

        assertEquals(3, db.getColumnCount());
        assertEquals(3, db.getRowCount());
        assertEquals("Tbilisi", db.getValueAt(0, 0));
        assertEquals("Prague", db.getValueAt(1, 0));
        assertEquals("Europe", db.getValueAt(1, 1));
    }

    @Test
    void testRowCount() throws SQLException {
        db.search("","", 13232132,
                POPULATION_PULLDOWN_LESSER_INDEX, MATCHTYPE_PULLDOWN_PARTIAL_INDEX);
        assertEquals(4, db.getRowCount());
    }

    @AfterAll
    void tearDown() throws SQLException {
        Statement stm = conn.createStatement();
        stm.execute("DROP DATABASE metropolis_test_db");
        conn.close();
    }
}