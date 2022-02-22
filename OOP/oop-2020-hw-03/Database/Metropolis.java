import javax.swing.*;

public class Metropolis implements MetropolisConstants {
    public Metropolis(String databaseName){
        MetropolisDatabase db = new MetropolisDatabase(databaseName);
        MetropolisTableModel model = new MetropolisTableModel(db);
        MetropolisFrame frame = new MetropolisFrame(model);
    }
}
