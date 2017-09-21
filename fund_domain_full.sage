#initialize
n = eval(raw_input("Please specify n: "))
N=binomial(n,2)

#list of unordered pairs of elements in {1,...,n}
def pairs(n): return flatten([[[i+1,j+1] for j in range(i+1,n)] for i in range(n-1)],max_level=1)

#sigma acts on an unordered pair
def act(sigma, p): return sorted([sigma(p[0]),sigma(p[1])])

#sigma acts as matrix in basis B on a vector v
def act_vect(sig,v,B=identity_matrix(N)): return (B.inverse()*sig.matrix()*B)*v

#permutation in Sigma_N which induced from perm in Sigma_n
def embed(sigma,n): return perm_from_sort([act(sigma,p) for p in pairs(n)]).to_cycles()

#find permutation which orders a list L	
def perm_from_sort(L): return Permutation([pair[0] for pair in sorted(enumerate(L, 1), key=lambda x: x[1])])

#image of generators of Sigma_n under embedding into Sigma_N
def embed_gens(n): return [embed(gen,n) for gen in SymmetricGroup(n).gens()]

#generate subgroup of permutation group G
G = PermutationGroup(embed_gens(n))

#fundamental domain given arbitrary distinct vector l and basis B
def fund_domain(l=[i for i in range(N)],B=identity_matrix(N),br=QQ):

	#augmented matrix for half-plane ineqs
	A = [ [l[j] - act_vect(G[i],vector(l),B)[j] for j in range(N)] for i in range(G.order())]
	b_1 = [0 for i in range(G.order())]
	bA = matrix(br,b_1).transpose().augment(matrix(br,A))

	#augmented matrix for positivity ineqs
	pos = B.inverse().transpose()
	b_2 = [0 for i in range(N)]
	bPos = matrix(br,b_2).transpose().augment(pos)

	#augmented matrix for fund. domain
	aug_l = bA.stack(bPos)
	
	#find convex polyhedral region, i.e. intersection of all ineqs
	poly = Polyhedron(ieqs=aug_l,base_ring=br)
	
	return([aug_l,poly])
	
#compactify a convex cone by adding a plane normal to x
def compactify(region,comp_dir=[1 for i in range(N)],chi=1,br=None):
	if br==None: br=region.base_ring()
	P_chi = [-comp_dir_i for comp_dir_i in comp_dir]
	P_chi.insert(0,chi)
	ieqs_list = [list(ieqs_i) for ieqs_i in region.inequalities()]
	ieqs_list.append(P_chi)
	return Polyhedron(ieqs=ieqs_list,base_ring=br)

#slice the fundamental domain with the plane x_1 + ... + x_N = 1
def cross_section(region,slice_dir=[1 for i in range(N)],chi=1,br=None):
	if br==None: br=region.base_ring()
	return Polyhedron(ieqs=region.inequalities(),eqns=[[-1]+slice_dir],base_ring=br)

#compactify fund. domain over QQ in direction comp_dir w/ cutoff chi
#and compute volume using optional Sage package 'lrslib'
def lrs_vol(region,comp_dir=[1 for i in range(N)],chi=1): return compactify(region,comp_dir,chi,QQ).volume(engine='lrs')

#orbit of point x under grounp action			
def orb(x): return list(set([tuple(act_vect(g,vector(x))) for g in G]))

#return orbit rep of x lying in F
def fund_domain_rep(x, F):
	for gx in orb(x):
		if F.contains(gx): return gx.list()

#size of stabilizer of point x
def stab_size(x): return G.order()/len(orb(x))

#returns interior point of polyhedron = average of vertices and rays 
#whose convex hull is polyhedron
def int_point(poly): 
	vert_rays = poly.vertices_list()+poly.rays_list()
	p = (1/len(vert_rays))*sum(vector(vr) for vr in vert_rays)
	#check
	if poly.relative_interior_contains(p): return p
	else: return "Failed."

#check if two faces of a polyhedron are glued together by group action
def faces_glued(f1,f2):
	p1 = int_point(f1.as_polyhedron())
	poly2 = f2.as_polyhedron()
	for o1 in orb(p1):
		if poly2.relative_interior_contains(o1): return True
	return False

#group element which glues face f1 to face f2
def gluing_elmt(f1,f2):
	p1 = int_point(f1.as_polyhedron())
	poly2 = f2.as_polyhedron()
	elmts = []
	for g in G:
		if poly2.relative_interior_contains(act_vect(g,vector(p1))): elmts.append(g)
	return elmts
	
#def faces_glued(f1,f2):
#	poly1 = f1.as_polyhedron()
#	poly2 = f2.as_polyhedron()
#	vr2 = [tuple(v) for v in poly2.rays()+poly2.vertices()]
#	for g in G:
#		act1 = [tuple(act_vect(g,vector(v))) for v in poly1.rays()+poly1.vertices()]
#		if set(act1) == set(vr2): return True
#	return False

#def gluing_elmt(f1,f2):
#	poly1 = f1.as_polyhedron()
#	poly2 = f2.as_polyhedron()
#	vr2 = [tuple(v) for v in poly2.rays()+poly2.vertices()]
#	elmts = []
#	for g in G:
#		act1 = [tuple(act_vect(g,vector(v))) for v in poly1.rays()+poly1.vertices()]
#		if set(act1) == set(vr2): elmts.append(g)
#	return elmts
	
#return list of lists of k-dim faces of polyhedron F after accounting for gluing by group
#action
def glued_face_lattice(F):
	glued_face_lat = []
	for k in range(F.dimension()+1):
		new_k_faces = []
		#faces which are glued are considered duplicate
		for face in F.faces(k):
			#check if face is already in new_k_faces
			contains_face = False
			for new_face in new_k_faces: 
				if(faces_glued(face,new_face)): contains_face = True
			if not contains_face: new_k_faces.append(face)
		glued_face_lat.append(new_k_faces)
	return glued_face_lat

#returns the Euler characteristic of a space
def euler_char(F): return sum([(-1)^k*len(F.faces(k)) for k in range(F.dim()+1)])

#returns orbifold Euler characteristic given fundamental domain F for orbifold X/G.
#takes into account glued faces as to not include extra terms in alternating sum.
def orb_euler_char(glued_lat):
	#initialize alternating sum
	alt_sum = 0
	#loop over faces of each dim. k
	for k in range(len(glued_lat)):
		for f in glued_lat[k]:
			#convert face to polyhedron object
			f_poly = f.as_polyhedron()
			#choose arbitrary interior point of face
			p = int_point(f_poly)
			#add summand
			alt_sum += (-1)^k*(1/stab_size(p))
	return alt_sum
	
#construct vector space containing polyhedron face
def tangent_vs(face):
	face_rays = [vector(v) for v in face.rays()]
	face_verts = [vector(v) for v in face.vertices()]
	diff = [[v2-v1 for v2 in face_verts] for v1 in face_verts]
	flat_diff = [item for sublist in diff for item in sublist]
	return span(face_rays+flat_diff)
	
#induced orientation of a boundary face
def outward_normal(bndry_face,face):
	#construct vector space containing face
	V_face=tangent_vs(face)
	
	#construct vector space containing bndry_face
	V_bndry=tangent_vs(bndry_face)
	
	#construct unit normal vector n as orthogonal complement of V_bndry in V_face
	n=V_bndry.basis_matrix().right_kernel().intersection(V_face).basis()[0]
	
	#choose outward pointing normal
	eps=10^(-32)
	p=int_point(bndry_face.as_polyhedron())
	s_plus= p + eps*n
	s_minus= p - eps*n
	if face.as_polyhedron().relative_interior_contains(s_plus): return -n
	if face.as_polyhedron().relative_interior_contains(s_minus): return n
	
	print 'Could not find outward pointing normal vector.'

#returns the induced orientation of bndry_face in face using an
#outward pointing normal vector wrt face
def induced_orient(bndry_face,face):
	#construct v.s. containing 'face'
	V_face=tangent_vs(face)
	#construct pos. oriented basis for bndry_face wrt face
	V_bndry=tangent_vs(bndry_face)
	n=outward_normal(bndry_face,face)
	bndry_basis=V_bndry.basis()
	
	coord_vecs = [V_face.coordinate_vector(n)]+[V_face.coordinate_vector(b) for b in bndry_basis]
	if matrix(coord_vecs).det() > 0: 
		if len(bndry_basis)==0: orient_basis = 1
		else: orient_basis = bndry_basis
	if matrix(coord_vecs).det() < 0: 
		if len(bndry_basis)==0: orient_basis = -1
		else: 
			orient_basis = [-bndry_basis[0]] + [bndry_basis[i] for i in range(1,len(bndry_basis))]
	
	return orient_basis
	
#compute boundary of a face before gluing
def bndry_up(face,F):
	face_bndry=[]
	for t in face.as_polyhedron().faces(k-1):
		for s in F.faces(k-1):
			if t.as_polyhedron()==s.as_polyhedron():
				face_bndry.append(s)
	return face_bndry

#returns boundary of polyhedron face inside F
def bndry(face,glued_lat,F):
	k=face.dim()
	
	#compute boundary of face before gluing
	face_bndry=bndry_up(face,F)
	
	bndry_vec=[0]*len(glued_lat[k-1])
	for f in face_bndry:
		for g in glued_lat[k-1]:
			#compute all group elements which map f to g
			sigma=gluing_elmt(f,g)
			if len(sigma) > 0:
				#tangent space to representative g of [f]=[g]
				V_g = tangent_vs(g)
				#tangent space to face f
				V_f = tangent_vs(f)
				
				#check if gluing is orientation preserving or reversing
				#by mapping the arbitrary basis assigned to f to the arbitrary
				#basis assigned to g. need to do this for each \sigma that maps
				#f to g. there are #stab([f]) of these elements. we then divide out
				#by this number to 'average' in the 'orbifold' sense
				orient_pres = 0
				for sig in sigma:
					glued_basis = [act_vect(sig,v) for v in V_f.basis()]
					glued_coord = [V_g.coordinate_vector(b) for b in glued_basis]
					orient_pres += sign(matrix(glued_coord).det())
				orient_pres = orient_pres/len(sigma)
					
				
				#check whether the induced orientation of f matches the arbitrary
				#orientation assigned to f. dimension 0 is a special case.
				induced_f_orient = induced_orient(f,face)
				if f.dim() == 0: 
					induced_match = induced_f_orient
				else: 
					f_coord = [V_f.coordinate_vector(b) for b in induced_f_orient]
					induced_match = sign(matrix(f_coord).det())
				
				#update the number of times the glued face g occurs in the boundary
				#accounting for signs
				indx = glued_lat[k-1].index(g)
				bndry_vec[indx] += orient_pres*induced_match
	
	return bndry_vec
	
#boundary map
def	bndry_map(k,glued_lat,F): return matrix(ZZ,[bndry(f,glued_lat,F) for f in glued_lat[k]]).transpose()

#cellular homology complex
def cell_chain_cmplx(glued_lat,F): return ChainComplex({k:bndry_map(k,glued_lat,F) for k in range(1,len(glued_lat))},degree=-1)

#boundary map in Z/2Z
def	mod2_bndry_map(k,glued_lat,F): return matrix(Integers(2),[bndry(f,glued_lat,F) for f in glued_lat[k]]).transpose()

#cellular homology in Z/2Z
def mod2_cell_chain_cmplx(glued_lat,F): return ChainComplex({k:mod2_bndry_map(k,glued_lat,F) for k in range(1,len(glued_lat))},degree=-1)
	

	









		
