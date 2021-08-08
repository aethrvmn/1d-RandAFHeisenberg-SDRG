# Ongoing Thesis Project

The goal of this code is to model the RG transformation for the 1d Random Antiferromagnetic Heisenberg Quantum Spin Chain.

# WIP
# Numerical Dasgupta-Ma RG

## Introduction to the Model

In order to better grasp the ideas presented in the previous chapter,
and also to appreciate the complexity of the Strong Disorder RG, we
created a numerical model to simulate the process using Python.  
If we wanted to simulate the nearest-neighbour spin chain we would need
two rows, one for each left spin and one for each right spin.
$$\\begin{aligned}
    \\text{left spin} :=\\ &\[1\\ 2\\ 3\\ 4\\ \\dots\\ N-1\]\\\\
    \\text{right spin} :=\\ &\[2\\ 3\\ 4\\ 5\\ \\dots\\ N\]    \\end{aligned}$$
where the bonds would then form a matrix
