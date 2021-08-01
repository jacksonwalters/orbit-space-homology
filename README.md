# orbit-space-homology
Given a convex polyhedral space X and a finite linear group G, computes the (simplicial) homology of the orbit space X/G. Computes a Dirichlet fundamental domain and refines triangulation by including extra vertices which are fixed by the action of G. 

This code was originally written for the case of unlabled networks with non-negative, real edge weights where X=(R_+)^N, N = n choose 2 = #edges, G=\Sigma_n \subset \Sigma_N, the edge permutations induced by permuting vertices. 

REQUIRED SOFTWARE:

Sage - http://www.sagemath.org/

lrslib - http://cgm.cs.mcgill.ca/~avis/C/lrs.html

SETUP:

This code is written in Python for the Sage mathematics software system. 
You will need to download, install, and configure the latest version of Sage, available at http://www.sagemath.org/.

The Polyhedra functions in Sage are able to make use of the following optional Sage package 'lrslib'.
This is most easily installed within Sage using the following terminal command: 

sage -i lrslib

Using this optional package will increase the speed these computations immensely, very roughly by a factor
of 10.

RUNTIME:

The main file is 'fund_domain.sage'. Sample Frechet means can be computed using 'sample_frechet_mean.sage'. To get started, change directories to the one containing this file. Enter the command

load('fund_domain.sage')

in the Sage environment. You will receive a prompt to set the global variable

n = your desired number of vertices

The global variable N = n choose 2. The global variable G=\Sigma_n will then be computed as a subgroup of the permutation group \Sigma_N. To construct a fundamental domain (for n=4, N=6) centered at the distinct vector l=[1,2,3,4,5,6.1], simply enter

F=fund_domain([1,2,3,4,5,6.1])[1]; F

F is a Sage polyhedron object representing a fundamental Dirichlet domain. Refer to the Sage documentation for the many available options for handling polyhedra:

http://doc.sagemath.org/html/en/reference/geometry/sage/geometry/polyhedron/constructor.html

To compute the sample Frechet mean of two networks,

load('sample_frechet_mean.sage')
s1=[1,2,3,4,5,6]
s2=[4,5,6,3,2,1]
mu=sample_frechet_mean(s1,s2); mu

The output will be a list of vectors.
