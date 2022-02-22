import nn


#autograder run times: 01:37:51, 01:33:92, 01:37:31, 02:10:75, 01:42:77, 02:12:05, 02:29:23, 01:52:17, 02:36:27, 01:37:53
#my time record with lower training bounds: 00:58:11

class PerceptronModel(object):
    def __init__(self, dimensions):
        """
        Initialize a new Perceptron instance.

        A perceptron classifies data points as either belonging to a particular
        class (+1) or not (-1). `dimensions` is the dimensionality of the data.
        For example, dimensions=2 would mean that the perceptron must classify
        2D points.
        """
        self.w = nn.Parameter(1, dimensions)

    def get_weights(self):
        """
        Return a Parameter instance with the current weights of the perceptron.
        """
        return self.w

    def run(self, x):
        """
        Calculates the score assigned by the perceptron to a data point x.

        Inputs:
            x: a node with shape (1 x dimensions)
        Returns: a node containing a single number (the score)
        """
        "*** YOUR CODE HERE ***"
        return nn.DotProduct(self.get_weights(), x)

    def get_prediction(self, x):
        """
        Calculates the predicted class for a single data point `x`.

        Returns: 1 or -1
        """
        "*** YOUR CODE HERE ***"
        prediction = 1.0
        dot_product = nn.as_scalar(self.run(x))
        if dot_product < 0:
            prediction = -1.0
        return prediction

    def train(self, dataset):
        """
        Train the perceptron until convergence.
        """
        "*** YOUR CODE HERE ***"
        converged = False
        while not converged:
            converged = True
            for x, y in dataset.iterate_once(1):
                y_star = self.get_prediction(x)
                y = nn.as_scalar(y)
                if y_star != y:
                    converged = False
                    self.w.update(x, y)


class RegressionModel(object):
    """
    A neural network model for approximating a function that maps from real
    numbers to real numbers. The network should be sufficiently large to be able
    to approximate sin(x) on the interval [-2pi, 2pi] to reasonable precision.
    """

    def __init__(self):
        # Initialize your model parameters here
        "*** YOUR CODE HERE ***"
        self.batch_size = 20
        self.feature_num = 1
        self.learning_rate = -0.08
        self.output_size = 1  # approximation of sin(x)
        self.layer_sizes = [10, 15, self.output_size]
        self.layer_num = len(self.layer_sizes)
        self.network_weights = []
        self.network_biases = []
        # initializes network weights and biases
        for i in range(self.layer_num):
            feature_num = self.feature_num if i == 0 else self.layer_sizes[i - 1]
            layer_size = self.layer_sizes[i]
            layer_weights = nn.Parameter(feature_num, layer_size)
            layer_biases = nn.Parameter(1, layer_size)
            self.network_weights.append(layer_weights)
            self.network_biases.append(layer_biases)

    def run(self, x):
        """
        Runs the model for a batch of examples.

        Inputs:
            x: a node with shape (batch_size x 1)
        Returns:
            A node with shape (batch_size x 1) containing predicted y-values
        """
        "*** YOUR CODE HERE ***"
        predictions = None
        for i in range(self.layer_num):
            layer_weights = self.network_weights[i]
            layer_biases = self.network_biases[i]
            xw = nn.Linear(x, layer_weights)
            predictions = nn.AddBias(xw, layer_biases)
            x = nn.ReLU(predictions)

        return predictions

    def get_loss(self, x, y):
        """
        Computes the loss for a batch of examples.

        Inputs:
            x: a node with shape (batch_size x 1)
            y: a node with shape (batch_size x 1), containing the true y-values
                to be used for training
        Returns: a loss node
        """
        "*** YOUR CODE HERE ***"
        return nn.SquareLoss(self.run(x), y)

    def train(self, dataset):
        """
        Trains the model.
        """
        "*** YOUR CODE HERE ***"
        total_loss = 0.03
        while total_loss > 0.02:
            for x, y in dataset.iterate_once(self.batch_size):
                curr_loss = self.get_loss(x, y)
                parameters = []
                for i in range(self.layer_num):
                    parameters.append(self.network_weights[i])
                    parameters.append(self.network_biases[i])
                grad = nn.gradients(curr_loss, parameters)
                for i in range(self.layer_num):
                    self.network_weights[i].update(grad[2 * i], self.learning_rate)
                    self.network_biases[i].update(grad[2 * i + 1], self.learning_rate)

            total_loss = nn.as_scalar(self.get_loss(nn.Constant(dataset.x), nn.Constant(dataset.y)))


class DigitClassificationModel(object):
    """
    A model for handwritten digit classification using the MNIST dataset.

    Each handwritten digit is a 28x28 pixel grayscale image, which is flattened
    into a 784-dimensional vector for the purposes of this model. Each entry in
    the vector is a floating point number between 0 and 1.

    The goal is to sort each digit into one of 10 classes (number 0 through 9).

    (See RegressionModel for more information about the APIs of different
    methods here. We recommend that you implement the RegressionModel before
    working on this part of the project.)
    """

    def __init__(self):
        # Initialize your model parameters here
        "*** YOUR CODE HERE ***"
        self.batch_size = 200  # 60
        self.feature_num = 784
        self.learning_rate = -0.525
        self.output_size = 10  # number of digits
        self.layer_sizes = [220, 150, self.output_size]
        # 170, 100, 10 50 -0.1 # 200 -0.5 180 120 10 # 200 200 150 10 -0.5
        self.layer_num = len(self.layer_sizes)
        self.network_weights = []
        self.network_biases = []
        # initializes network weights and biases
        for i in range(self.layer_num):
            feature_num = self.feature_num if i == 0 else self.layer_sizes[i - 1]
            layer_size = self.layer_sizes[i]
            layer_weights = nn.Parameter(feature_num, layer_size)
            layer_biases = nn.Parameter(1, layer_size)
            self.network_weights.append(layer_weights)
            self.network_biases.append(layer_biases)

    def run(self, x):
        """
        Runs the model for a batch of examples.

        Your model should predict a node with shape (batch_size x 10),
        containing scores. Higher scores correspond to greater probability of
        the image belonging to a particular class.

        Inputs:
            x: a node with shape (batch_size x 784)
        Output:
            A node with shape (batch_size x 10) containing predicted scores
                (also called logits)
        """
        "*** YOUR CODE HERE ***"
        predictions = None
        for i in range(self.layer_num):
            layer_weights = self.network_weights[i]
            layer_biases = self.network_biases[i]
            xw = nn.Linear(x, layer_weights)
            predictions = nn.AddBias(xw, layer_biases)
            x = nn.ReLU(predictions)

        return predictions

    def get_loss(self, x, y):
        """
        Computes the loss for a batch of examples.

        The correct labels `y` are represented as a node with shape
        (batch_size x 10). Each row is a one-hot vector encoding the correct
        digit class (0-9).

        Inputs:
            x: a node with shape (batch_size x 784)
            y: a node with shape (batch_size x 10)
        Returns: a loss node
        """
        "*** YOUR CODE HERE ***"
        return nn.SoftmaxLoss(self.run(x), y)

    def train(self, dataset):
        """
        Trains the model.
        """
        "*** YOUR CODE HERE ***"
        validation_accuracy = 0
        while validation_accuracy < 0.975:
            for x, y in dataset.iterate_once(self.batch_size):
                curr_loss = self.get_loss(x, y)
                parameters = []
                for i in range(self.layer_num):
                    parameters.append(self.network_weights[i])
                    parameters.append(self.network_biases[i])
                grad = nn.gradients(curr_loss, parameters)
                for i in range(self.layer_num):
                    self.network_weights[i].update(grad[2 * i], self.learning_rate)
                    self.network_biases[i].update(grad[2 * i + 1], self.learning_rate)
            validation_accuracy = dataset.get_validation_accuracy()


class LanguageIDModel(object):
    """
    A model for language identification at a single-word granularity.

    (See RegressionModel for more information about the APIs of different
    methods here. We recommend that you implement the RegressionModel before
    working on this part of the project.)
    """

    def __init__(self):
        # Our dataset contains words from five different languages, and the
        # combined alphabets of the five languages contain a total of 47 unique
        # characters.
        # You can refer to self.num_chars or len(self.languages) in your code
        self.num_chars = 47
        self.languages = ["English", "Spanish", "Finnish", "Dutch", "Polish"]

        # Initialize your model parameters here
        "*** YOUR CODE HERE ***"
        self.batch_size = 30
        self.feature_num = self.num_chars
        self.learning_rate = -0.05
        self.output_size = len(self.languages)
        self.layer_sizes = [250, 250, self.output_size] # 60, 0.05, 250, 250
        self.layer_num = len(self.layer_sizes)
        self.network_weights = []
        self.network_biases = []
        # initializes network weights and biases
        for i in range(self.layer_num):
            feature_num = self.feature_num if i == 0 else self.layer_sizes[i - 1]
            layer_size = self.layer_sizes[i]
            layer_weights = nn.Parameter(feature_num, layer_size)
            layer_biases = nn.Parameter(1, layer_size)
            self.network_weights.append(layer_weights)
            self.network_biases.append(layer_biases)

    def run(self, xs):
        """
        Runs the model for a batch of examples.

        Although words have different lengths, our data processing guarantees
        that within a single batch, all words will be of the same length (L).

        Here `xs` will be a list of length L. Each element of `xs` will be a
        node with shape (batch_size x self.num_chars), where every row in the
        array is a one-hot vector encoding of a character. For example, if we
        have a batch of 8 three-letter words where the last word is "cat", then
        xs[1] will be a node that contains a 1 at position (7, 0). Here the
        index 7 reflects the fact that "cat" is the last word in the batch, and
        the index 0 reflects the fact that the letter "a" is the inital (0th)
        letter of our combined alphabet for this task.

        Your model should use a Recurrent Neural Network to summarize the list
        `xs` into a single node of shape (batch_size x hidden_size), for your
        choice of hidden_size. It should then calculate a node of shape
        (batch_size x 5) containing scores, where higher scores correspond to
        greater probability of the word originating from a particular language.

        Inputs:
            xs: a list with L elements (one per character), where each element
                is a node with shape (batch_size x self.num_chars)
        Returns:
            A node with shape (batch_size x 5) containing predicted scores
                (also called logits)
        """
        "*** YOUR CODE HERE ***"
        # Batch_size x feature_size * feature_size x layer_size = batch_size x layer_size
        latest_prediction = self.getInitialPrediction(xs[0])
        latest_prediction = nn.ReLU(latest_prediction)
        for x in xs[1:]:
            # (Batch_size x feature_size * feature_size x layer_size) +
            # + (batch_size x layer_size * self.network_weights[1]) =
            # = batch_size x layer_size because for add function both matrices should have same dimensions
            # from that we can see that (batch_size x layer_size * self.network_weights[1]) = batch_size x layer_size
            # and self.network_weights[1] should be layer_size x layer_size
            # after that there are two ways to return batch_size x 5.
            # first way is that first and second layer weights should be 47x5 and 5x5
            # that means that we cannot have neither more nor less than 5 perceptrons in each layer which isnt optimal
            # second way is that we have one more layer that converts previous matrices to batch_size x 5
            # in this way we can have first layer weights of size 47xN, second layer NxN and third layer Nx5
            initial_prediction = self.getInitialPrediction(x)
            latest_prediction = nn.Add(initial_prediction, nn.Linear(latest_prediction, self.network_weights[1]))
            latest_prediction = nn.ReLU(latest_prediction)

        final_prediction = nn.Linear(latest_prediction, self.network_weights[2])
        return final_prediction

    def get_loss(self, xs, y):
        """
        Computes the loss for a batch of examples.

        The correct labels `y` are represented as a node with shape
        (batch_size x 5). Each row is a one-hot vector encoding the correct
        language.

        Inputs:
            xs: a list with L elements (one per character), where each element
                is a node with shape (batch_size x self.num_chars)
            y: a node with shape (batch_size x 5)
        Returns: a loss node
        """
        "*** YOUR CODE HERE ***"
        return nn.SoftmaxLoss(self.run(xs), y)

    def train(self, dataset):
        """
        Trains the model.
        """
        "*** YOUR CODE HERE ***"
        validation_accuracy = 0
        while validation_accuracy < 0.87:
            for x, y in dataset.iterate_once(self.batch_size):
                curr_loss = self.get_loss(x, y)
                parameters = []
                for i in range(self.layer_num):
                    parameters.append(self.network_weights[i])
                grad = nn.gradients(curr_loss, parameters)
                for i in range(self.layer_num):
                    self.network_weights[i].update(grad[i], self.learning_rate)
            validation_accuracy = dataset.get_validation_accuracy()

    def getInitialPrediction(self, x):
        layer_weights = self.network_weights[0]
        layer_biases = self.network_biases[0]
        xw = nn.Linear(x, layer_weights)
        predictions = nn.AddBias(xw, layer_biases)
        return predictions
