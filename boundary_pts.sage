#return orbit rep of x lying in F
def fund_domain_rep(x, F):
	orb = [act_vect(g,vector(x)) for g in G]
	for gx in orb:
		if F.contains(gx): return gx.list()

#return some point on the boundary of F
def boundary_pt(F,ind=0):
	rays=F.faces(N-1)[ind].rays()
	return sum(vector(ray) for ray in rays)/(N-1)
	
#push a boundary point a distance eps in the direction normal to its
#containing facet
def near_boundary_pt(eps,F,ind=0):
	#obtain vector normal to subspace containing the N-1 dim. facet
	rays=F.faces(N-1)[ind].rays()
	vec_rays=[vector(ray) for ray in rays]
	A=matrix(vec_rays).transpose()
	perp=A.kernel()
	n=perp.basis()[0]
	n_hat=n/n.norm()
	
	#push boundary point distance eps in +-n_hat direction
	#and take rep in F
	p=boundary_pt(F,ind)
	s_plus=p+eps*n_hat
	s_minus=p-eps*n_hat
	if F.contains(s_plus): return s_plus
	if F.contains(s_minus): return s_minus
	print "Failed."
	
#test if p is on the boundary of F
def test_bndry(p,F): return F.contains(p) and not(F.interior_contains(p))
	
#check if p is in the interior of its boundary face
def face_int(p,F,ind=0): return F.faces(N-1)[ind].as_polyhedron().relative_interior_contains(p)


def check_all_ints(eps,F):
	num_facets=len(F.faces(N-1))
	s_list=[near_boundary_pt(eps,F,i) for i in range(num_facets)]
	F_list=[fund_domain(vector(s,RDF),B=identity_matrix(6),br=QQ)[1] for s in s_list]
	R_list=[F.intersection(F_i) for F_i in F_list]
	int_list=[[R_i.intersection(R_j) for R_j in R_list] for R_i in R_list]
	for i in range(num_facets):
		for j in range(num_facets):
			avg=vector(s_list[i])+vector(s_list[j])
			R=R_list[i].intersection(R_list[j])
			if not R.contains(avg): return [[i,j],int_list]
	return int_list
	
	
	
