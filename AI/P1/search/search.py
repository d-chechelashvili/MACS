# search.py
# ---------
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


"""
In search.py, you will implement generic search algorithms which are called by
Pacman agents (in searchAgents.py).
"""

import util


class SearchProblem:
    """
    This class outlines the structure of a search problem, but doesn't implement
    any of the methods (in object-oriented terminology: an abstract class).

    You do not need to change anything in this class, ever.
    """

    def getStartState(self):
        """
        Returns the start state for the search problem.
        """
        util.raiseNotDefined()

    def isGoalState(self, state):
        """
          state: Search state

        Returns True if and only if the state is a valid goal state.
        """
        util.raiseNotDefined()

    def getSuccessors(self, state):
        """
          state: Search state

        For a given state, this should return a list of triples, (successor,
        action, stepCost), where 'successor' is a successor to the current
        state, 'action' is the action required to get there, and 'stepCost' is
        the incremental cost of expanding to that successor.
        """
        util.raiseNotDefined()

    def getCostOfActions(self, actions):
        """
         actions: A list of actions to take

        This method returns the total cost of a particular sequence of actions.
        The sequence must be composed of legal moves.
        """
        util.raiseNotDefined()


def tinyMazeSearch(problem):
    """
    Returns a sequence of moves that solves tinyMaze.  For any other maze, the
    sequence of moves will be incorrect, so only use this for tinyMaze.
    """
    from game import Directions
    s = Directions.SOUTH
    w = Directions.WEST
    return [s, s, w, s, w, w, s, w]


def depthFirstSearch(problem):
    """
    Search the deepest nodes in the search tree first.

    Your search algorithm needs to return a list of actions that reaches the
    goal. Make sure to implement a graph search algorithm.

    To get started, you might want to try some of these simple commands to
    understand the search problem that is being passed in:

    print("Start:", problem.getStartState())
    print("Is the start a goal?", problem.isGoalState(problem.getStartState()))
    print("Start's successors:", problem.getSuccessors(problem.getStartState()))
    """
    "*** YOUR CODE HERE ***"
    visitedStates = list()
    fringe = util.Stack()
    fringe.push((problem.getStartState(), []))  # second argument for memorizing path

    # neighbor[0] and currState[0] are just coordinates for respective states
    # currState[1] is full path from start to that state
    # neighbor[1] is direction where that neighbor is
    while not fringe.isEmpty():
        currState = fringe.pop()
        if problem.isGoalState(currState[0]):
            return currState[1]
        if not visitedStates.__contains__(currState[0]):
            visitedStates.append(currState[0])
            for neighbor in problem.getSuccessors(currState[0]):
                fringe.push((neighbor[0], currState[1] + [neighbor[1]]))
    return []


def breadthFirstSearch(problem):
    """Search the shallowest nodes in the search tree first."""
    "*** YOUR CODE HERE ***"

    visitedStates = list()
    fringe = util.Queue()
    fringe.push((problem.getStartState(), []))  # second argument for memorizing path

    # neighbor[0] and currState[0] are just coordinates for respective states
    # currState[1] is full path from start to that state
    # neighbor[1] is direction where that neighbor is
    while not fringe.isEmpty():
        currState = fringe.pop()
        if problem.isGoalState(currState[0]):
            return currState[1]
        if not visitedStates.__contains__(currState[0]):
            visitedStates.append(currState[0])
            for neighbor in problem.getSuccessors(currState[0]):
                fringe.push((neighbor[0], currState[1] + [neighbor[1]]))
    return []


def uniformCostSearch(problem):
    """Search the node of least total cost first."""
    "*** YOUR CODE HERE ***"
    visitedStates = list()
    fringe = util.PriorityQueue()
    fringe.push((problem.getStartState(), [], 0),
                0)  # second argument for memorizing path, third argument for path cost

    # neighbor[0] and currState[0] are just coordinates for respective states
    # currState[1] is full path from start to that state
    # neighbor[1] is direction where that neighbor is
    # currState[2] is total cost of the currState[1] path
    while not fringe.isEmpty():
        currState = fringe.pop()
        if problem.isGoalState(currState[0]):
            return currState[1]
        if not visitedStates.__contains__(currState[0]):
            visitedStates.append(currState[0])
            for neighbor in problem.getSuccessors(currState[0]):
                totalPathCost = currState[2] + neighbor[2]
                fringe.push((neighbor[0], currState[1] + [neighbor[1]], totalPathCost), totalPathCost)
    return []


def nullHeuristic(state, problem=None):
    """
    A heuristic function estimates the cost from the current state to the nearest
    goal in the provided SearchProblem.  This heuristic is trivial.
    """
    return 0


def aStarSearch(problem, heuristic=nullHeuristic):
    """Search the node that has the lowest combined cost and heuristic first."""
    "*** YOUR CODE HERE ***"
    visitedStates = list()
    fringe = util.PriorityQueue()
    fringe.push((problem.getStartState(), [], 0),
                0)  # second argument for memorizing path, third argument for path cost

    # neighbor[0] and currState[0] are just coordinates for respective states
    # currState[1] is full path from start to that state
    # neighbor[1] is direction where that neighbor is
    # currState[2] is total cost of the currState[1] path
    while not fringe.isEmpty():
        currState = fringe.pop()
        if problem.isGoalState(currState[0]):
            return currState[1]
        if not visitedStates.__contains__(currState[0]):
            visitedStates.append(currState[0])
            for neighbor in problem.getSuccessors(currState[0]):
                totalPathCost = currState[2] + neighbor[2]
                fringe.push((neighbor[0], currState[1] + [neighbor[1]], totalPathCost),
                            totalPathCost + heuristic(neighbor[0], problem))
    return []


# Abbreviations
bfs = breadthFirstSearch
dfs = depthFirstSearch
astar = aStarSearch
ucs = uniformCostSearch
