#basis for 1-dim irrep.
v_b=vector([1,1,1,1,1,1])
v_o=(1/sqrt(6))*vector([1,1,1,1,1,1])

#bases for 2-dim irrep.
h1=vector([1,-1,0,0,-1,1])
h2=vector([0,1,-1,-1,1,0])
subB_2=matrix([h1,h2]).transpose()

#gram-schmidt to orthonormalize above
h_bar_1=vector([1/2,-1/2,0,0,-1/2,1/2])
h_bar_2=vector([1/(2*sqrt(3)),1/(2*sqrt(3)),-1/sqrt(3),-1/sqrt(3),1/(2*sqrt(3)),1/(2*sqrt(3))])
subO_2=matrix([h_bar_1,h_bar_2]).transpose()

#basiss for 3-dim irrep.
l1=vector([0,1,1,-1,-1,0])
l2=vector([1,-1,0,0,1,-1])
l3=vector([0,-1,1,-1,1,0])
subB_3=matrix([l1,l2,l3]).transpose()

#gram-schmidt to orthonormalize above
l_bar_1=vector([0,1/2,1/2,-1/2,-1/2,0])
l_bar_2=vector([1/sqrt(3),-1/(2*sqrt(3)),1/(2*sqrt(3)),-1/(2*sqrt(3)),1/(2*sqrt(3)),-1/sqrt(3)])
l_bar_3=vector([-(1/sqrt(6)),-(1/sqrt(6)),(1/sqrt(6)),-(1/sqrt(6)),(1/sqrt(6)),(1/sqrt(6))])
subO_3=matrix([l_bar_1,l_bar_2,l_bar_3]).transpose()

#basis for whole space
B=matrix([v_b,h1,h2,l1,l2,l3]).transpose()
O=matrix([v_o,h_bar_1,h_bar_2,l_bar_1,l_bar_2,l_bar_3]).transpose()

#sigma acts as matrix in basis B on a vector v
def act_vect(sig,v,B=identity_matrix(N)): return (B.inverse()*sig.matrix()*B)*v

#matrix representation of perm. sig in given basis
def mat_rep(sig,B=identity_matrix(N)): return B.inverse()*sig.matrix()*B

#helper function
def all_mats(B=identity_matrix(N)): return [mat_rep(sig,B) for sig in G]

#build subgroup H
sig0=G.identity()
sig1=SymmetricGroup(4).list()[21]
sig2=SymmetricGroup(4).list()[3]
sig3=SymmetricGroup(4).list()[13]

H=PermutationGroup([embed(sig0,4),embed(sig1,4),embed(sig2,4),embed(sig3,4)])

#coset representatives for H
cosets=G.cosets(H)

#coset representatives for G/H
coset_reps = [coset_rep[0] for coset_rep in cosets]

#return all matrix reps. in V_2 wrt basis
def all_mats_subgroup(B=subB_2): return [mat_rep(coset_rep,B) for coset_rep in coset_reps]



	
	
