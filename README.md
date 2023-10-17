# Quantum Monte Carlo

## Background

Quantum Monte Carlo (QMC) uses the Markov Chain Monte Carlo (MCMC) algorithm, in particular, the Metropolis-Hastings algorithm, to simulate a quantum system by modelling it as a chain of coupled classical particles of length $L$. QMC builds off of Feynman's path integral formalism, which highlights the link between classical and quantum mechanics. 

The basis of the project used the following equation: 

$$S(x) = \epsilon \sum_{j=0}^{L-1} \left[\frac{m}{2} \left\frac{x_{j+1}-x_{j}}{\epsilon}\right^{2} + \frac{1}{2}kx_{j}^{2} \right]$$

where $\epsilon = \frac{1}{LT}$, $m$ is the mass, $x_{j}$ is the initial position of the particle and $k$ is the spring constant. $S(x)$ is known as the action, defining the chain of coupled particles. THe first term is derived from the kinetic part of the Hamiltonian and the second term denotes the harmonic 'spring' coupling between the particles. $L$ refers to how many particles are in the chain. 

## Results

The results produced from the program are expected, and match known theory, with the $L=20$ case trending towards the $L=\infty$ graph. 

## Conclusion and Future Work

The MCMC algorithm can be used to find the statistics of a quantum particle with good agreement with theory. An increase in $L$ and the number of iterations in the algorithm should produce results with less error and agree even closer with theory. Further study could result in the quantum anharmonic oscillator being used instead of basing the project off the quantum harmonic oscilator. 

## What This Repository Contains

The main code file, Quantum Monte Carlo Project [Final Code].py is all that is required to run the code. The results folder includes a variety of plots showcasing the different values of $L$. 
