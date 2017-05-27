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
def embed(sigma,n): return perm_from_sort([act(sigma,p) for p in pairs(n)])

#find permutation which orders a list L	
def perm_from_sort(L): return Permutation([pair[0] for pair in sorted(enumerate(L, 1), key=lambda x: x[1])]).to_cycles()

#image of generators of Sigma_n under embedding into Sigma_N
def embed_gens(n): return [embed(gen,n) for gen in SymmetricGroup(n).gens()]

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

	
#translation decomposition of fundamental domain F_l at x	
def fund_decomp(x,l=[i for i in range(N)],br=QQ):

	VR = []
	
	for k in range(G.order()):
		
		#fund. domain centered at l
		F_l = fund_domain(l)[0]
		
		#fund. domain centered at sigma.x
		sig_x = [x[G[k].inverse()(i+1)-1] for i in range(len(x))]
		F_sig_x = fund_domain(sig_x)[0]
		
		#intersection of F_l and F_sig_x
		int = F_l.stack(F_sig_x)
		
		#find convex polyhedral region corresponding to intersection
		poly = Polyhedron(ieqs=int,base_ring=br)
		
		#check uniqueness
		isUnique = True
		for i in range(len(VR)):
			if poly == VR[i][0]: 
				isUnique = False
				VR[i][1].append(G[k])
		if isUnique == True:
			VR.append([poly,[G[k]]])
	
				
	return(VR)
	
#compactify a convex cone by adding a plane normal to x
def compactify(region,comp_dir=[1 for i in range(N)],chi=10,br=None):
	if br==None: br=region.base_ring()
	P_chi = [-comp_dir_i for comp_dir_i in comp_dir]
	P_chi.insert(0,chi)
	ieqs_list = [list(ieqs_i) for ieqs_i in region.inequalities()]
	ieqs_list.append(P_chi)
	return Polyhedron(ieqs=ieqs_list,base_ring=br)

#compactify fund. domain over QQ in direction comp_dir w/ cutoff chi
#and compute volume using optional Sage package 'lrslib'
def lrs_vol(region,comp_dir=[1 for i in range(N)],chi=10): return compactify(region,comp_dir,chi,QQ).volume(engine='lrs')

#given Vrepresentation of fund. domain computed in QQ, approximate  
#normalized rays over RR and remove duplicates
def vRep_approx(F_x,prec=53): 
	ray_approx=list(set([tuple(vector(RealField(prec),v)/norm(vector(RealField(prec),v))) for v in F_x.ray_generator()]))
	vert_approx=[vector(RealField(prec),v) for v in F_x.vertices_list()]
	return ray_approx+vert_approx


#compare two fundamental domains
def compare_domains(F_1,F_2):
	
	#domains
	print('Domains: \n F_1: %r \n F_2: %r' % (F_1, F_2))
	
	#check to see if they're the same
	if F_1 == F_2: print('Same Regions in R^%r' % N)
	else: print('Different Regions in R^%r' % N)
	
	#volumes
	print('Volumes: \n V_1=%r \n V_2=%r' % (lrs_vol(F_1),lrs_vol(F_2)))
	
	#number of planes in H-rep
	print('#H-rep: \n #F_1=%r \n #F_2=%r' % (F_1.n_Hrepresentation(),F_2.n_Hrepresentation()))
			
#return orbit rep of x lying in F
def fund_domain_rep(x, F):
	orb = [act_vect(g,vector(x)) for g in G]
	for gx in orb:
		if F.contains(gx): return gx.list()


#extract info from decomp vector
def dim_vector(decomp): return sorted([R[0].dim() for R in decomp])
def polys(decomp): return [R[0] for R in decomp]
def perms(decomp): return [R[1] for R in decomp]
def facet_graphs(decomp): return [Graph(R[0].facet_adjacency_matrix()) for R in decomp]
def top_dim_cells(decomp): 
	L=[]; 
	for i in range(len(decomp)):
		if decomp[i][0].dim() == N: L.append(decomp[i]) 
	return L
		
	









		
