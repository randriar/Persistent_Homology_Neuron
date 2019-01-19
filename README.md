# Persistent_Homology_Neuron
Independent Research with the Mathematics Department at Lafayette College.
I had worked with a professor along this project (Summer 2017-Spring 2018).<br>
I presented this research at the LVAIC poster Symposium at Moravian College, PA on November 2017.


TITLE: Using Persistent Homology to study Neuron’s morphology
<br><br>
ABSTRACT:<br>
Our project aims to use topology to study the shape of a neuron and hopefully gain insight into the neuron’s biological function. We use persistent homology, a popular tool from topological data analysis that can keep track of features across all possible scales, to analyze the morphology of the cell. Currently, our project aims to implement a MatLab program that will compute the persistent homology of a neuron from a sampling of its spatial coordinates.
<br><br>
From re-centering the cell, computing the Kernel Density Estimator (KDE), modeling spherical shells by KDEs, and creating explicit simplex streams, our program computes barcodes at each interval to record the persistent features that intersect each spherical shell. The homology and the intervals of persistence are computed by using the ‘Javaplex’ package in MatLab. The barcodes are determined by the intervals of persistence which represent an important feature of the neuron: the scale at which a component of the neuron is born and dies. We have also explored using the histogram method to measure the density of the neuron in place of the traditional KDE. This method allows us to focus our attention more directly on the neuron, without over-smoothing.
<br><br>
We also examine the relationship between our work and the traditional Sholl Analysis approach to morphology. We looked at the Schoenen Ramification Index and Branching Index. These types of Sholl Analysis are computed from the number of spines that intersect each spherical radius. This data can be directly recovered from our presentation of the persistent homology.
