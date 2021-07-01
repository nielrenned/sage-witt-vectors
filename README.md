# Sage Witt Vectors

This is a SageMath implementation of the algorithms in Luís Finotti's paper *Computations with Witt Vectors and the Greenberg Transform*. There is an implementation in Magma [here](https://github.com/lrfinotti/witt).

# Documentation

We add two classes: `WittRing` and `WittVector`, which behave like usual Sage objects. More documentation to come.

# State of the Code Base

This is still in the early stages. We have yet to optimize the operations, e.g. convert some functions to Cython, or profile the code. We also have a goal of using multiprocessing to get even more improvements.

We have some unit tests, which are comparisons with the Magma codebase linked above.