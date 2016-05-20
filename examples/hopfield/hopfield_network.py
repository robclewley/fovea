import numpy as np
import random

class HopfieldNetwork:
    """
    (C) Daniel McNeela, 2016

    Implements the Hopfield Network, a recurrent neural network developed by John Hopfield
    circa 1982.

    c.f. https://en.wikipedia.org/wiki/Hopfield_Network
    """
    def __init__(self, num_neurons):
        """
        Instantiates a Hopfield Network comprised of "num_neurons" neurons.
        
        num_neurons         The number of neurons in the network.
        _weights            The network's weight matrix.
        _trainers           A dictionary containing the methods available for 
                            training the network.
        _vec_activation     A vectorized version of the network's activation function.
        """
        self.num_neurons = num_neurons
        self._weights = np.zeros((self.num_neurons, self.num_neurons), dtype=np.int_)
        self._trainers = {"hebbian": self._hebbian, "storkey": self._storkey}
        self._vec_activation = np.vectorize(self._activation)

    def weights(self):
        """
        Getter method for the network's weight matrix.
        """
        return self._weights

    def reset(self):
        """
        Reset's the network's matrix to the matrix which is identically zero.

        Useful for retraining the network from scratch after an initial round
        of training has already been completed.
        """
        self._weights = np.zeros((self.num_neurons, self.num_neurons), dtype=np.int_)

    def train(self, patterns, method="hebbian", threshold=0):
        """
        The wrapper function for the network's various training methods stored in
        self._trainers.

        patterns        A list of the on which to train the network. Patterns are
                        bipolar vectors of the form 

                        [random.choice([-1, 1]) for i in range(self.num_neurons)].

                        Example of properly formatted input for a Hopfield Network
                        containing three neurons:

                            [[-1, 1, 1], [1, -1, 1]]

        method          The training algorithm to be used. Defaults to "hebbian".
                        Look to self._trainers for a list of the available options.
        threshold       The threshold value for the network's activation function.
                        Defaults to 0.
        """
        try:
            return self._trainers[method](patterns, threshold)
        except KeyError:
            print(method + " is not a valid learning method.")

    def learn(self, patterns, steps=None):
        """
        To be used after training the network.

        Given 'patterns', learn(patterns) classifies these patterns based on those
        which the network has already seen.
        """
        if steps:
            for i in range(steps):
                learn = np.dot(patterns, self._weights)
            return self._vec_activation(learn)
        else:
            pre_learn = patterns
            while True:
                post_learn = [np.dot(pattern, self._weights) for pattern in patterns]
                if np.array_equal(pre_learn, post_learn):
                    return self._vec_activation(post_learn)
                pre_learn = post_learn

    def _activation(self, value, threshold=0):
        """
        The network's activation function.

        Defaults to the sign function.
        """
        if value < threshold:
            return -1
        return 1

    def _hebbian(self, patterns, threshold=0):
        """
        Implements Hebbian learning.
        """
        for pattern in patterns:
            self._weights += np.outer(pattern, pattern)
        self._weights = self._vec_activation(self._weights, threshold)
        np.fill_diagonal(self._weights, 0)

    def _storkey(self, patterns):
        """
        Implements Storkey learning.
        """
        pass
