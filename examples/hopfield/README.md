Included in this directory is a visualization of a Hopfield Network applied to an example in character recognition.

The file `hopfield_network.py` contains the definition of the `HopfieldNetwork` class. 
Using this class, Hopfield Networks containing any desired number of constituent neurons can be constructed from scratch. 
The given network implementation provides the following methods intended for the end user:

	weights() 						Returns the network's current weight matrix.
	
	reset()							Resets the network's weight matrix so that the network can be retrained on new data.

	train(patterns,					Trains the network on a series of input states passed as "patterns." Two training methods
		  method="hebbian",			are provided: hebbian and storkey.
		  threshold=0,
		  inject=lambda x: None)

	learn(patterns,					Causes the network to classify the states given by "patterns" based on the states upon
		  steps=None,				which it was trained. Two learning methods are provided: synchronous and asynchronous.
		  mode="asynchronous",		Synchronous learning is the faster of the two, but is not guaranteed to converge.
		  inject=lambda x: None)	Asynchronous learning is rather slow, but is guaranteed to converge to a local low-energy state.

	energy(state)					Calculates the energy associated with the given state.

The file `visuals.py` contains the code for running the network visualization and all associated helper functions for drawing
individual components of that visualization to the Matplotlib canvas. The primary definition of the file is that of the
`VisualHopfield` class. This defines a "visual" Hopfield Network that subclasses the implementation given in hopfield_network.py.

The only method from this class which end users should concern themselves with is `run_visualization(training_data, learning_data=None)`.
This method calls on all the internally-defined helper methods to run a full visualization of the network training on the provided
`training_data` and learning the provided `learning_data`. Thus to run the visualization, do the following...

	$ python -i visuals.py
	
	>>> myNet = VisualHopfield(num_neurons)			# Replace num_neurons with your desired number of neurons.

	>>> training_data = [your_data_here]			# Construct a list of training_data

	>>> learning_data = [learning_data_here]		# Optionally, construct a list of learning data.

	>>>	myNet.run_visualization(training_data, learning_data)

Alternatively, to view the visualization as applied to an example in character recognition, run the ocr.py script.

There are five subplots associated with the network visualization:

1. The first of these is the main network diagram. This displays each neuron as a circle connected to every other
neuron in the network as is the convention within the Hopfield paradigm. A green line between two neurons
represents a connection of weight 0, a blue line a connection of weight 1, and a red line a connection of weight -1.
During the training portion of the visualization, the connections that are in the process of being altered in response
to training data are highlighted by way of a thicker linewidth and change in color in tandem with the changing of the
network's weight matrix.

2. Second is the energy function diagram. This provides a plot of the network's energy landscape after it has been trained.
   Since for nearly all cases the network state vectors are of dimension greater than three, the energy landscape is calculated
   by applying the network's energy function to the two-dimensional PCA axes of the network's training vectors. The energy
   function provides a measure of the error in the  
