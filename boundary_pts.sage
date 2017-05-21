#return orbit rep of x lying in F
def fund_domain_rep(x, F):
	orb = [act_vect(g,vector(x)) for g in G]
	for gx in orb:
		if F.contains(gx): return gx.list()

#return some point on the boundary of F
def boundary_pt(F):
	rays=F.faces(N-1)[0].rays()
	return sum(vector(ray) for ray in rays)/(N-1)
	
def test_bndry(p,F): return F.contains(p) and not(F.interior_contains(p))
	
#check if p is in the interior of its boundary face
def face_int(p,F): return F.faces(N-1)[0].as_polyhedron().relative_interior_contains(p)
	
def near_boundary_pt(eps,F):
	#obtain vector normal to subspace containing the N-1 dim. facet
	rays=F.faces(N-1)[0].rays()
	vec_rays=[vector(ray) for ray in rays]
	A=matrix(vec_rays).transpose()
	perp=A.kernel()
	n=perp.basis()[0]
	n_hat=n/n.norm()
	
	#push boundary point distance eps in n_hat direction
	#and take rep in F
	p=boundary_pt(F)
	s_tilde=p+eps*n_hat
	return fund_domain_rep(s_tilde,F)
	