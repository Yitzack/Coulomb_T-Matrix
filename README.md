# Coulomb_T-Matrix
Code base for a project to calculate the spectral function and T-matrix of the coulomb potential.

The simple point of this project is to calculate the spectral function of two fermion particles attracted to each other by the Coulomb potential with as few other effects as possible at arbitrary momentum. The grand scheme of physics research, this isn't particularly useful because we can do the same thing with Schrodinger's Equation and with much less work. But that is what is so great. We know what the answer is supposed to be. The question isn't what the answer ought be, but how do we get that answer. That is the actual point of this project. The point of this project is to determine the mathematical mechanisms required to do this for arbitrary central potentials.

## Around.h
This a C++ header and object that is to function as a quantity with uncertainty.

## Elements.h
I'm not sure I need this, but it enables the programmed parallel manipulation of an array of data. It doesn't get a parallel speed up on purpose. However, the compiler may make a parallel optimization for the CPU.

## Interpolation.h
This is an object for the storage and evaluation of a bicubic spline interpolation. It only needs the data and the size and it assumes that it is evenly spaced. It takes `Interpolation(<T>**Data, int XLength, int YLength)` and will return a bicubic interpolation spline object of the data with type `<T>`.
