# Computable-Performance-Bounds

This repository contains the code and report of my project for the course CS 754: Advanced Image Processing. It is based on the paper titled "Computable Performance Bounds on Sparse Recovery", by Gongguo Tang and Arye Nehorai.

# Why? And what?
In relation to sparse recovery based on l1 minimization, existing performance bounds are either hard to verify (such as the Restricted Isometry Property and the Nullspace Property) or give very loose bounds (as in the case of Mutual Coherence-based bounds). This work introduces a novel set of bounds which are simultaneously 
a) easily verifiable (can be calculated using 'n' convex optimizations where the matrix used for compression is of dimension m x n).

b) highly tight (upto an order of magnitude tighter than RIP-based bounds, which were estimated using Monte Carlo approximations: the MC approximations yield a bound _tighter_ than the actual RIP-based bound)

c) Widely applicable (the necessary conditions for the bounds to hold are more encompassing than those for RIP and MC-based bounds)

This work uses the spgl1 and CVX solvers for convex optimization.
