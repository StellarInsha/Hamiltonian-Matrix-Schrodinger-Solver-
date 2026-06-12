
## Why I Built This

In introductory quantum mechanics, many systems are solved analytically because they appear to be mathematically convenient. Real physical systems are often not easy to solve by hands.

This project revolves around how computers can solve the one-dimensional time-independent Schrödinger equation numerically by constructing the Hamiltonian matrix and finding its eigenvalues & eigenvectors.




### Harmonic Oscillator

A numerical solution of the quantum harmonic oscillator, used as a benchmark because its exact analytical solution is known.

### Anharmonic Oscillator

A harmonic oscillator with an additional x⁴ term. The energy levels are no longer equally spaced as seen in harmonic oscillator.

### Double Well Potential

A system with two potential minima separated by a barrier. The lowest energy states exhibit quantum tunneling.

### Finite Square Well

A finite-depth potential well that supports bound states with negative energies.



## Method Used

The continuous Schrödinger equation is discretized on a spatial grid.

The kinetic energy operator is represented using a finite-difference approximation to the second derivative, while the potential energy is represented as a diagonal matrix.

Together they form the Hamiltonian matrix:

H = T + V

Diagonalizing the Hamiltonian gives:

* Energy eigenvalues
* Wavefunctions (eigenstates)



## What I Learned from this

* Finite-difference methods
* Hamiltonian matrix formulation
* Eigenvalue problems
* Numerical solutions of quantum systems
* The connection between linear algebra and quantum mechanics




This project was completed as part of my self-study of quantum mechanics and computational physics.
