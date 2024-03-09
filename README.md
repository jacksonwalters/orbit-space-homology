# orbit-space-homology
Given a convex polyhedral space X and a finite linear group G, computes the (simplicial) (co)homology of the orbit space X/G. Computes a Dirichlet fundamental domain and refines triangulation by including extra vertices which are fixed by the action of G. 

This code was originally written for the case of unlabled networks with non-negative, real edge weights where $X=R_{\ge 0}^N$, $N = \binom{n}{2}$ = #edges, $G=\Sigma_n \subset \Sigma_N$, the edge permutations induced by permuting vertices. 

REQUIRED SOFTWARE:

Sage - http://www.sagemath.org/

lrslib - http://cgm.cs.mcgill.ca/~avis/C/lrs.html

SETUP:

This code is written in Python for the SageMath software system. 
You will need to install and configure the latest version of SageMath, available at http://www.sagemath.org/.

The Polyhedra functions in Sage are able to make use of the following optional Sage package 'lrslib'.
This is most easily installed within Sage using the following terminal command: 

sage -i lrslib

Using this optional package will increase the speed these computations immensely, very roughly by a factor
of 10.

RUNTIME:

To get started, open

```orbit_space_simplicial_homology.ipynb```

in the Sage environment, i.e. a Jupyter notebook with the Sage kernel. You will receive a prompt to set the global variable

n = your desired number of vertices

The global variable $N = \binom{n}{2}$. The global variable $G=\Sigma_n$ will then be computed as a subgroup of the permutation group $\Sigma_N$. To construct a fundamental domain (for $n=4$, $N=6$) centered at the distinct vector $l=[1,2,3,4,5,6.1]$, simply enter

$$F=\text{fund_domain}([1,2,3,4,5,6.1])[1]; F$$

$F$ is a Sage polyhedron object representing a fundamental Dirichlet domain. Refer to the Sage documentation for the many available options for handling polyhedra:

http://doc.sagemath.org/html/en/reference/geometry/sage/geometry/polyhedron/constructor.html

To compute the cohomology, run the rest of the notebook. It will compute the $k$-faces for each $k$ from $1,\dots,N$ (very slow for $N \ge 6$). It will then assemble the boundary maps, and build the chain complex. The homology is $\text{im}(\partial_{k+1})/\text{ker}(\partial_{k})$.

G-CW complexes: https://math.mit.edu/research/undergraduate/urop-plus/documents/2016/Liu.pdf
