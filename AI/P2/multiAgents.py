# multiAgents.py
# --------------
# Licensing Information:  You are free to use or extend these projects for
# educational purposes provided that (1) you do not distribute or publish
# solutions, (2) you retain this notice, and (3) you provide clear
# attribution to UC Berkeley, including a link to http://ai.berkeley.edu.
# 
# Attribution Information: The Pacman AI projects were developed at UC Berkeley.
# The core projects and autograders were primarily created by John DeNero
# (denero@cs.berkeley.edu) and Dan Klein (klein@cs.berkeley.edu).
# Student side autograding was added by Brad Miller, Nick Hay, and
# Pieter Abbeel (pabbeel@cs.berkeley.edu).


import random
import util

from game import Agent


def euclideanDistance(xy1, xy2):
    return ((xy1[0] - xy2[0]) ** 2 + (xy1[1] - xy2[1]) ** 2) ** 0.5


class ReflexAgent(Agent):
    def getGhostMultiplier(self, newPos, newGhostStates, newScaredTimes):
        danger = 35
        distToGhost = euclideanDistance(newPos, newGhostStates[0].getPosition())
        if (newScaredTimes[0] > 0):
            danger += 0
        # danger -= newScaredTimes[0] * (35 - manhattanDistance(newPos, newGhostStates[0].getPosition()))
        elif (distToGhost <= 4):
            danger += 160 * (8 - distToGhost)
        return danger

    def getFoodMultiplier(self, foodGrid, position):
        foodCount = 0
        closestDistance = foodGrid.height + foodGrid.width + 1
        for i in range(foodGrid.height - 1, -1, -1):
            for j in range(foodGrid.width):
                if foodGrid[j][i]:
                    foodCount += 1
                    currPoint = (j, i)
                    currDistance = util.manhattanDistance(position, currPoint)
                    if currDistance < closestDistance:
                        closestDistance = currDistance

        return 100 / max(closestDistance, 1)

    """
    A reflex agent chooses an action at each choice point by examining
    its alternatives via a state evaluation function.

    The code below is provided as a guide.  You are welcome to change
    it in any way you see fit, so long as you don't touch our method
    headers.
    """

    def getAction(self, gameState):
        """
        You do not need to change this method, but you're welcome to.

        getAction chooses among the best options according to the evaluation function.

        Just like in the previous project, getAction takes a GameState and returns
        some Directions.X for some X in the set {NORTH, SOUTH, WEST, EAST, STOP}
        """
        # Collect legal moves and successor states
        legalMoves = gameState.getLegalActions()

        # Choose one of the best actions
        scores = [self.evaluationFunction(gameState, action) for action in legalMoves]
        bestScore = max(scores)
        bestIndices = [index for index in range(len(scores)) if scores[index] == bestScore]
        chosenIndex = random.choice(bestIndices)  # Pick randomly among the best

        "Add more of your code here if you want to"

        return legalMoves[chosenIndex]

    def evaluationFunction(self, currentGameState, action):
        """
        Design a better evaluation function here.

        The evaluation function takes in the current and proposed successor
        GameStates (pacman.py) and returns a number, where higher numbers are better.

        The code below extracts some useful information from the state, like the
        remaining food (newFood) and Pacman position after moving (newPos).
        newScaredTimes holds the number of moves that each ghost will remain
        scared because of Pacman having eaten a power pellet.

        Print out these variables to see what you're getting, then combine them
        to create a masterful evaluation function.
        """
        # Useful information you can extract from a GameState (pacman.py)
        successorGameState = currentGameState.generatePacmanSuccessor(action)
        newPos = successorGameState.getPacmanPosition()
        newFood = successorGameState.getFood()
        newGhostStates = successorGameState.getGhostStates()
        newScaredTimes = [ghostState.scaredTimer for ghostState in newGhostStates]

        # food: +10
        # scared ghost: +200
        # win: +500
        # lose: -500
        # time penalty: -1
        "*** YOUR CODE HERE ***"
        probableScore = self.getFoodMultiplier(currentGameState.getFood(), newPos) \
                        - self.getGhostMultiplier(newPos, newGhostStates, newScaredTimes)
        # print(probableScore)
        return probableScore


def scoreEvaluationFunction(currentGameState):
    """
    This default evaluation function just returns the score of the state.
    The score is the same one displayed in the Pacman GUI.

    This evaluation function is meant for use with adversarial search agents
    (not reflex agents).
    """
    return currentGameState.getScore()


class MultiAgentSearchAgent(Agent):
    """
    This class provides some common elements to all of your
    multi-agent searchers.  Any methods defined here will be available
    to the MinimaxPacmanAgent, AlphaBetaPacmanAgent & ExpectimaxPacmanAgent.

    You *do not* need to make any changes here, but you can if you want to
    add functionality to all your adversarial search agents.  Please do not
    remove anything, however.

    Note: this is an abstract class: one that should not be instantiated.  It's
    only partially specified, and designed to be extended.  Agent (game.py)
    is another abstract class.
    """

    def __init__(self, evalFn='scoreEvaluationFunction', depth='2'):
        self.index = 0  # Pacman is always agent index 0
        self.evaluationFunction = util.lookup(evalFn, globals())
        self.depth = int(depth)


class MinimaxAgent(MultiAgentSearchAgent):
    """
    Your minimax agent (question 2)
    """

    def maxValue(self, gameState, depth):
        import math
        val = -math.inf
        legalActions = gameState.getLegalActions(0)
        bestAction = legalActions[0]
        for action in legalActions:
            successor = gameState.generateSuccessor(0, action)
            successorValue = self.value(successor, depth, 1)
            if successorValue > val:
                bestAction = action
                val = successorValue
        return val, bestAction

    def minValue(self, gameState, depth, index):
        import math
        val = math.inf

        newIndex = (index + 1) % gameState.getNumAgents()
        newDepth = depth
        if newIndex == 0:
            newDepth -= 1

        legalActions = gameState.getLegalActions(index)
        for action in legalActions:
            successor = gameState.generateSuccessor(index, action)
            successorValue = self.value(successor, newDepth, newIndex)
            val = min(val, successorValue)
        return val

    def value(self, gameState, depth, index):
        if depth == 0 or gameState.isWin() or gameState.isLose():
            return self.evaluationFunction(gameState)
        if index == 0:
            return self.maxValue(gameState, depth)[0]
        return self.minValue(gameState, depth, index)

    def getAction(self, gameState):
        """
        Returns the minimax action from the current gameState using self.depth
        and self.evaluationFunction.

        Here are some method calls that might be useful when implementing minimax.

        gameState.getLegalActions(agentIndex):
        Returns a list of legal actions for an agent
        agentIndex=0 means Pacman, ghosts are >= 1

        gameState.generateSuccessor(agentIndex, action):
        Returns the successor game state after an agent takes an action

        gameState.getNumAgents():
        Returns the total number of agents in the game

        gameState.isWin():
        Returns whether or 0not the game state is a winning state

        gameState.isLose():
        Returns whether or not the game state is a losing state
        """
        "*** YOUR CODE HERE ***"
        return self.maxValue(gameState, self.depth)[1]
        # util.raiseNotDefined()


class AlphaBetaAgent(MultiAgentSearchAgent):
    """
    Your minimax agent with alpha-beta pruning (question 3)
    """

    def maxValue(self, gameState, depth, a, b):
        import math
        val = -math.inf
        legalActions = gameState.getLegalActions(0)
        bestAction = legalActions[0]
        for action in legalActions:
            successor = gameState.generateSuccessor(0, action)
            successorValue = self.value(successor, depth, 1, a, b)
            if successorValue > val:
                bestAction = action
                val = successorValue
            if val > b:
                return val, bestAction
            a = max(a, val)
        return val, bestAction

    def minValue(self, gameState, depth, index, a, b):
        import math
        val = math.inf

        newIndex = (index + 1) % gameState.getNumAgents()
        newDepth = depth
        if newIndex == 0:
            newDepth -= 1

        legalActions = gameState.getLegalActions(index)
        for action in legalActions:
            successor = gameState.generateSuccessor(index, action)
            successorValue = self.value(successor, newDepth, newIndex, a, b)
            val = min(val, successorValue)
            if val < a:
                return val
            b = min(b, val)
        return val

    def value(self, gameState, depth, index, a, b):
        if depth == 0 or gameState.isWin() or gameState.isLose():
            return self.evaluationFunction(gameState)
        if index == 0:
            return self.maxValue(gameState, depth, a, b)[0]
        return self.minValue(gameState, depth, index, a, b)

    def getAction(self, gameState):
        """
        Returns the minimax action using self.depth and self.evaluationFunction
        """
        "*** YOUR CODE HERE ***"
        import math
        a = -math.inf
        b = math.inf
        return self.maxValue(gameState, self.depth, a, b)[1]
        # util.raiseNotDefined()


class ExpectimaxAgent(MultiAgentSearchAgent):
    """
      Your expectimax agent (question 4)
    """

    def maxValue(self, gameState, depth):
        import math
        val = -math.inf
        legalActions = gameState.getLegalActions(0)
        bestAction = legalActions[0]
        for action in legalActions:
            successor = gameState.generateSuccessor(0, action)
            successorValue = self.value(successor, depth, 1)
            if successorValue > val:
                bestAction = action
                val = successorValue
        return val, bestAction

    def expValue(self, gameState, depth, index):
        val = 0

        newIndex = (index + 1) % gameState.getNumAgents()
        newDepth = depth
        if newIndex == 0:
            newDepth -= 1

        legalActions = gameState.getLegalActions(index)
        for action in legalActions:
            successor = gameState.generateSuccessor(index, action)
            probability = 1 / len(legalActions)
            successorValue = self.value(successor, newDepth, newIndex)
            val += probability * successorValue
        return val

    def value(self, gameState, depth, index):
        if depth == 0 or gameState.isWin() or gameState.isLose():
            return self.evaluationFunction(gameState)
        if index == 0:
            return self.maxValue(gameState, depth)[0]
        return self.expValue(gameState, depth, index)

    def getAction(self, gameState):
        """
        Returns the expectimax action using self.depth and self.evaluationFunction

        All ghosts should be modeled as choosing uniformly at random from their
        legal moves.
        """
        "*** YOUR CODE HERE ***"
        return self.maxValue(gameState, self.depth)[1]
        # util.raiseNotDefined()


def getBetterFoodMultiplier(currentGameState):
    #ეს სთეითი პოტენციურად ღირს 500 + ქურსქორი + 10 * გეთნამფუდი - რამდენიც მოძრაობაც დარჩა ამ ფუდების ასაღებად
    #foodScore = 500
    #foodScore += 10 * getNumFood(currentGameState.getFood())
    # problem  = FoodSearchProblem(currentGameState)
    # stepCount = len(breadthFirstSearch(problem))
    # foodScore -= stepCount #ცოტა მეტი უნდა დავაკლო რადგან გოსთი შეუშლის ხელს?

    #my birthday 13/12 passes barely - 1009.5
    #operation overlord 1945/6/6 - 1041.9
    #everything between 57390 - 73035.335093050925934 gives 1063.8
    foodScore = 73035.335093050925934/max(1, currentGameState.getNumFood())
    return foodScore


def getBetterGhostMultiplier(currentGameState):
    ghostStates = currentGameState.getGhostStates()
    scaredTimes = [ghostState.scaredTimer for ghostState in ghostStates]
    pos = currentGameState.getPacmanPosition()

    danger = 0
    euclideanDistToGhost = euclideanDistance(pos, ghostStates[0].getPosition())
    if scaredTimes[0] > 0:
        if scaredTimes[0] > euclideanDistToGhost and euclideanDistToGhost <= 5:
            danger -= max(215 - 50 * euclideanDistToGhost, 10)
    elif euclideanDistToGhost <= 4:
        danger = 115 * (5.5 - euclideanDistToGhost)
    return danger


def betterEvaluationFunction(currentGameState):
    """
    Your extreme ghost-hunting, pellet-nabbing, food-gobbling, unstoppable
    evaluation function (question 5).

    DESCRIPTION: <write something here so we know what you did>
    evaluation function should give us estimate of score we will get from this point by playing optimally
    """
    "*** YOUR CODE HERE ***"
    # food: +10
    # scared ghost: +200
    # win: +500
    # lose: -500
    # time penalty: -1
    stateValue = currentGameState.getScore() + getBetterFoodMultiplier(currentGameState) \
                 - getBetterGhostMultiplier(currentGameState)
    return stateValue
    # util.raiseNotDefined()

better = betterEvaluationFunction
