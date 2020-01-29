import acm.util.RandomGenerator;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

public class NameSurferColors {

	/**
	 * creates a new NameSurferColors object which is responsible for painting every
	 * graph with different colour
	 */
	public NameSurferColors() {

	}

	/**
	 * by mixing red green and blue gets a colour which is useful for painting
	 * graphs by using special algorithm
	 */
	public Color getAvailableColor() {
		int[] RGB = new int[3];
		Color colour = getColor(RGB);
		boolean available = checkAvailability(colour, RGB);
		if (!available) {
			colour = getNew(colour, RGB);
		}
		return colour;
	}

	// gets a colour by randomly mixing red green and blue
	private Color getColor(int[] RGB) {
		for (int i = 0; i < 3; i++) {
			int x = rgen.nextInt(256);
			RGB[i] = x;
		}
		Color colour = new Color(RGB[0], RGB[1], RGB[2]);
		return colour;
	}

	// checks if retrieved colour is useful for us by checking whether it is used
	// or is too light for us to use
	private boolean checkAvailability(Color colour, int[] RGB) {
		if (!(usedColors.containsKey(colour) || tooLight(RGB))) {
			usedColors.put(colour, RGB);
			return true;
		}
		return false;
	}

	// calculates the luminance and checks if the retrieved colour is too light for
	// us to use
	// if colour is too light it will be hard to see on white background so its
	// not useful
	private boolean tooLight(int[] RGB) {
		int R = RGB[0];
		int G = RGB[1];
		int B = RGB[2];
		double luminance = 0.2126 * R + 0.7152 * G + 0.0722 * B;
		if (luminance > 185) {
			return true;
		}
		return false;
	}

	// gets the new colour with a special algorithm for both cases
	private Color getNew(Color colour, int[] RGB) {
		if (usedColors.containsKey(colour)) {
			colour = changeCoincident(colour, RGB);
		} else {
			colour = changeLuminance(colour, RGB);
		}
		return colour;
	}

	// if there is a same colour already used it should have same red green and blue
	// values
	// this method changes one of them randomly and tests if new colour is useful
	private Color changeCoincident(Color colour, int[] RGB) {
		int i = rgen.nextInt(3);
		RGB[i] = rgen.nextInt(256);
		colour = new Color(RGB[0], RGB[1], RGB[2]);
		boolean available = checkAvailability(colour, RGB);
		if (!available) {
			colour = getNew(colour, RGB);
		}
		return colour;
	}

	// this method changes luminance so that the retrieved colour is easy to see on
	// white background
	// because value of blue has very small effect on luminance this method only
	// chooses from values of red and green
	// also because 0.7152 of luminance comes from green's value changing it will
	// have the most effect, so this method has
	// 0.7152 chance of changing green's value
	private Color changeLuminance(Color colour, int[] RGB) {
		boolean RG = rgen.nextBoolean(0.7152);
		if (RG) {
			RGB[1] = rgen.nextInt(256);
		} else {
			RGB[0] = rgen.nextInt(256);
		}
		colour = new Color(RGB[0], RGB[1], RGB[2]);
		boolean available = checkAvailability(colour, RGB);
		if (!available) {
			colour = getNew(colour, RGB);
		}
		return colour;
	}

	/**
	 * Clears the list of used colours stored inside this class.
	 */
	public void clear() {
		usedColors.clear();
	}

	/**
	 * Clears the colour stored inside this class that is associated with the entry
	 * user removed.
	 */
	public void clearSelected(Color color) {
		usedColors.remove(color);
	}

	private RandomGenerator rgen = new RandomGenerator();
	private Map<Color, int[]> usedColors = new HashMap<Color, int[]>();
}