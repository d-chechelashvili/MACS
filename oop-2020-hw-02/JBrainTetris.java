
import javax.swing.*;
import java.awt.*;

public class JBrainTetris extends JTetris{
    private JCheckBox brainMode;
    private JCheckBox animateFallMode;
    private JPanel little;
    private JSlider adversary;
    private JLabel statusLabel;
    private DefaultBrain brain;
    private Brain.Move bestMove;
    private int pieceCount;

    /**
     * Creates a new JTetris where each tetris square
     * is drawn with the given number of pixels.
     *
     * @param pixels pix
     */
    JBrainTetris(int pixels) {
        super(pixels);
        brain = new DefaultBrain();
    }


    @Override
    public JComponent createControlPanel() {
        JPanel panel = (JPanel)super.createControlPanel();
        panel.add(new JLabel("Brain:"));
        brainMode = new JCheckBox("Brain active");
        panel.add(brainMode);
        animateFallMode = new JCheckBox("Animated fall active");
        panel.add(animateFallMode);
        animateFallMode.setSelected(true);
        little = new JPanel();
        little.add(new JLabel("Adversary:"));
        adversary = new JSlider(0, 100, 0); // min, max, current
        adversary.setPreferredSize(new Dimension(100,15));
        little.add(adversary);
        statusLabel = new JLabel("ok");
        little.add(statusLabel);
        panel.add(little);
        return panel;
    }

    @Override
    public void tick(int verb) {
        if(brainMode.isSelected()){
            if(pieceCount != super.count) {
                super.board.undo();
                pieceCount = super.count;
                bestMove = brain.bestMove(super.board, super.currentPiece, HEIGHT, bestMove);
            }

            if(bestMove != null && verb == JTetris.DOWN){
                if(bestMove.x < super.currentX) super.tick(JTetris.LEFT);
                else if(bestMove.x > super.currentX) super.tick(JTetris.RIGHT);
                if(!bestMove.piece.equals(currentPiece)) super.tick(JTetris.ROTATE);
                if(bestMove.piece.equals(currentPiece)
                        && bestMove.x == super.currentX && !animateFallMode.isSelected() && bestMove.y < super.currentY){
                    super.tick(JTetris.DROP);
                }
            }
        }

        super.tick(verb);
    }

    @Override
    public Piece pickNextPiece() {
        int randomNum = random.nextInt(100);

        if(randomNum >= adversary.getValue()){
            statusLabel.setText("ok");
            return super.pickNextPiece();
        }

        statusLabel.setText("*ok*");
        return getWorstPiece();
    }

    private Piece getWorstPiece() {
        Brain.Move currMove = new Brain.Move();
        int maxScoreIndex = 0;
        double maxScore = 0;
        for (int i = 0; i < super.pieces.length; i++) {
            currMove = brain.bestMove(super.board, super.pieces[i], HEIGHT, currMove);
            if(currMove == null) continue;
            double currScore = currMove.score;
            if(currScore > maxScore){
                maxScore = currScore;
                maxScoreIndex = i;
            }
        }
        return super.pieces[maxScoreIndex];
    }

    public static void main(String[] args) {
        // Set GUI Look And Feel Boilerplate.
        // Do this incantation at the start of main() to tell Swing
        // to use the GUI LookAndFeel of the native platform. It's ok
        // to ignore the exception.
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception ignored) { }

        JBrainTetris tetris = new JBrainTetris(16);
        JFrame frame = JBrainTetris.createFrame(tetris);
        frame.setVisible(true);
    }
}
