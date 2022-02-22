# valueIterationAgents.py
# -----------------------
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


# valueIterationAgents.py
# -----------------------
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


import mdp, util

from learningAgents import ValueEstimationAgent
import collections


class ValueIterationAgent(ValueEstimationAgent):
    """
        * Please read learningAgents.py before reading this.*

        A ValueIterationAgent takes a Markov decision process
        (see mdp.py) on initialization and runs value iteration
        for a given number of iterations using the supplied
        discount factor.
    """

    def __init__(self, mdp, discount=0.9, iterations=100):
        """
          Your value iteration agent should take an mdp on
          construction, run the indicated number of iterations
          and then act according to the resulting policy.

          Some useful mdp methods you will use:
              mdp.getStates()
              mdp.getPossibleActions(state)
              mdp.getTransitionStatesAndProbs(state, action)
              mdp.getReward(state, action, nextState)
              mdp.isTerminal(state)
        """
        self.mdp = mdp
        self.discount = discount
        self.iterations = iterations
        self.values = util.Counter()  # A Counter is a dict with default 0
        self.runValueIteration()

    def runValueIteration(self):
        # Write value iteration code here
        "*** YOUR CODE HERE ***"
        i = 0
        while i < self.iterations:
            self.computeValues()
            i += 1

    def computeValues(self):
        states = self.mdp.getStates()
        newValues = util.Counter()
        for state in states:
            newValues[state] = self.computeHighestQValue(state)
        self.values = newValues

    def computeHighestQValue(self, state):
        if self.mdp.isTerminal(state):
            return 0

        import math
        actions = self.mdp.getPossibleActions(state)
        bestValue = -math.inf
        for action in actions:
            currValue = self.computeQValueFromValues(state, action)
            if currValue > bestValue:
                bestValue = currValue
        return bestValue

    def getValue(self, state):
        """
          Return the value of the state (computed in __init__).
        """
        return self.values[state]

    def computeQValueFromValues(self, state, action):
        """
          Compute the Q-value of action in state from the
          value function stored in self.values.
        """
        "*** YOUR CODE HERE ***"
        if action is None:
            return 0

        transitions = self.mdp.getTransitionStatesAndProbs(state, action)
        QValue = 0
        for transition in transitions:
            nextState = transition[0]
            probability = transition[1]
            QValue += probability * (self.mdp.getReward(state, action, nextState) +
                                     self.discount * self.getValue(nextState))

        return QValue
        # util.raiseNotDefined()

    def computeActionFromValues(self, state):
        """
          The policy is the best action in the given state
          according to the values currently stored in self.values.

          You may break ties any way you see fit.  Note that if
          there are no legal actions, which is the case at the
          terminal state, you should return None.
        """
        "*** YOUR CODE HERE ***"
        import math
        actions = self.mdp.getPossibleActions(state)
        bestAction = None
        bestValue = -math.inf
        for action in actions:
            currValue = self.computeQValueFromValues(state, action)
            if currValue > bestValue:
                bestValue = currValue
                bestAction = action

        return bestAction
        # util.raiseNotDefined()

    def getPolicy(self, state):
        return self.computeActionFromValues(state)

    def getAction(self, state):
        "Returns the policy at the state (no exploration)."
        return self.computeActionFromValues(state)

    def getQValue(self, state, action):
        return self.computeQValueFromValues(state, action)


class AsynchronousValueIterationAgent(ValueIterationAgent):
    """
        * Please read learningAgents.py before reading this.*

        An AsynchronousValueIterationAgent takes a Markov decision process
        (see mdp.py) on initialization and runs cyclic value iteration
        for a given number of iterations using the supplied
        discount factor.
    """

    def __init__(self, mdp, discount=0.9, iterations=1000):
        """
          Your cyclic value iteration agent should take an mdp on
          construction, run the indicated number of iterations,
          and then act according to the resulting policy. Each iteration
          updates the value of only one state, which cycles through
          the states list. If the chosen state is terminal, nothing
          happens in that iteration.

          Some useful mdp methods you will use:
              mdp.getStates()
              mdp.getPossibleActions(state)
              mdp.getTransitionStatesAndProbs(state, action)
              mdp.getReward(state)
              mdp.isTerminal(state)
        """
        ValueIterationAgent.__init__(self, mdp, discount, iterations)

    def runValueIteration(self):
        "*** YOUR CODE HERE ***"
        i = 0
        states = self.mdp.getStates()
        while i < self.iterations:
            state = states[i % len(states)]
            self.values[state] = self.computeHighestQValue(state)
            i += 1


class PrioritizedSweepingValueIterationAgent(AsynchronousValueIterationAgent):
    """
        * Please read learningAgents.py before reading this.*

        A PrioritizedSweepingValueIterationAgent takes a Markov decision process
        (see mdp.py) on initialization and runs prioritized sweeping value iteration
        for a given number of iterations using the supplied parameters.
    """

    def __init__(self, mdp, discount=0.9, iterations=100, theta=1e-5):
        """
          Your prioritized sweeping value iteration agent should take an mdp on
          construction, run the indicated number of iterations,
          and then act according to the resulting policy.
        """
        self.theta = theta
        ValueIterationAgent.__init__(self, mdp, discount, iterations)

    def runValueIteration(self):
        "*** YOUR CODE HERE ***"
        predecessors = self.computePredecessors()
        queue = util.PriorityQueue()
        states = self.mdp.getStates()

        for state in states:
            if self.mdp.isTerminal(state):
                continue
            diff = abs(self.values[state] - self.computeHighestQValue(state))
            queue.push(state, -diff)

        i = 0
        while i < self.iterations:
            if queue.isEmpty():
                break
            currState = queue.pop()
            self.values[currState] = self.computeHighestQValue(currState)
            for p in predecessors[currState]:
                diff = abs(self.values[p] - self.computeHighestQValue(p))
                if diff > self.theta:
                    queue.update(p, -diff)
            i += 1

    def computePredecessors(self):
        allPredecessors = {}
        states = self.mdp.getStates()
        for state in states:
            actions = self.mdp.getPossibleActions(state)
            for action in actions:
                transitions = self.mdp.getTransitionStatesAndProbs(state, action)
                successorStates = [transition[0] for transition in transitions]
                for successorState in successorStates:
                    allPredecessors.setdefault(successorState, set()).add(state)
        return allPredecessors



    # def computePredecessors(self):
    #     allPredecessors = {}
    #     visited = set()
    #     states = self.mdp.getStates()
    #     for state in states:
    #         if state not in visited:
    #             currPredecessors = set()
    #             self.computePredecessorsUtil(state, currPredecessors, visited, allPredecessors)
    #     return allPredecessors
    #
    # def computePredecessorsUtil(self, state, currPredecessors, visited, allPredecessors):
    #     if self.mdp.isTerminal(state) or state in visited:
    #         return
    #     visited.add(state)
    #     allPredecessors.setdefault(state, set()).union(currPredecessors)
    #     currPredecessors.add(state)
    #     actions = self.mdp.getPossibleActions(state)
    #     for action in actions:
    #         transitions = self.mdp.getTransitionStatesAndProbs(state, action)
    #         successorStates = [transition[0] for transition in transitions]
    #         for successorState in successorStates:
    #             self.computePredecessorsUtil(successorState, currPredecessors, visited, allPredecessors)
