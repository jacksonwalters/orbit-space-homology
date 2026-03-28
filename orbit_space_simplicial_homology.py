n = 4
N=binomial(n,2)
#N=4

#list of unordered pairs of elements in {1,...,n}
def pairs(n): return flatten([[[i+1,j+1] for j in 
                            range(i+1,n)] for i in range(n-1)],max_level=1)

#sigma acts on an unordered pair
def act(sigma, p): return sorted([sigma(p[0]),sigma(p[1])])

#sigma acts as matrix in basis on a vector v
def act_vect(g,v,basis=identity_matrix(N)): 
    return (basis.inverse()*g.matrix()*basis)*vector(v)

#permutation in Sigma_N which induced from perm in Sigma_n
def embed(sigma,n): return perm_from_sort([act(sigma,p) for p in pairs(n)])

#find permutation which orders a list L	
def perm_from_sort(L): return Permutation([pair[0] for pair in 
                                           sorted(enumerate(L, 1), key=lambda x: x[1])]).to_cycles()

#image of generators of Sigma_n under embedding into Sigma_N
def embed_gens(n): return [embed(gen,n) for gen in SymmetricGroup(n).gens()]

# generators of Sigma_n as transpositions
def swaps(n): return [SymmetricGroup(n)([(i,i+1)]) for i in range(1,n)]
def embedded_swaps(n): return [SymmetricGroup(binomial(n,2))(embed(gen,n)) for gen in swaps(n)]

#return n when N=(n choose 2)
def invert_binomial(N):
    for n in range(floor(sqrt(2*N)),floor(sqrt(2*N)+2)):
        if binomial(n,2) == N:
            return n
        
#angle between two vectors
def angle(x,y): return float(vector(x)*vector(y)/(vector(x).norm()*vector(y).norm()))

#naive minimum finding acting over entire group
def find_min(x,y):
    assert len(x) == len(y)
    print(angle(x,y))
    n = invert_binomial(len(x))
    G = PermutationGroup(embed_gens(n))
    for sigma in G:
        if vector(x)*(sigma.matrix()*vector(y)) >= vector(x)*vector(y):
            y = sigma.matrix()*vector(y)
    print(angle(x,y))
    return y

#fundamental domain given arbitrary distinct vector l and basis B
def fund_domain(center=[i for i in range(N)],basis=identity_matrix(N),br=QQ,group=MatrixGroup(identity_matrix(N))):

    #augmented matrix for half-plane ineqs
    A = [ [center[j] - act_vect(group[i],vector(center),basis)[j] for j in range(N)] for i in range(group.order())]
    b_1 = group.order()*[0]
    bA = matrix(br,b_1).transpose().augment(matrix(br,A))

    #augmented matrix for positivity ineqs
    pos = basis.inverse().transpose()
    b_2 = [0 for i in range(N)]
    bPos = matrix(br,b_2).transpose().augment(pos)

    #augmented matrix for fund. domain
    aug_l = bA.stack(bPos)
    
    #find convex polyhedral region, i.e. intersection of all ineqs
    poly = Polyhedron(ieqs=aug_l,base_ring=br)

    return([aug_l,poly])

#slice the fundamental domain with the plane x_1 + ... + x_N = 1
def cross_section(region,slice_dir=[1 for i in range(N)],chi=1,br=None):
    if br is None: br=region.base_ring()
    return Polyhedron(ieqs=region.inequalities(),eqns=[[-1]+slice_dir],base_ring=br)

#find fixed points under G in polyhedron F to include as new vertices
#find intersection of fixed point subspace with F as polyhedron, and return those vertices
def fixed_verts(F,br=QQ):
    new_verts=set()
    for g in G:
        #get equations for fixed point subspace
        A=g.matrix() #matrix associated to group element g
        B=A-identity_matrix(N) #eqns defining fixed pt subspace are (A-Id)x == 0
        b = N*[0] #0 vector
        bB = matrix(br,b).transpose().augment(matrix(br,B))#form augmented matrix for equations defining subspace
        #get equations for intersection
        eqns=matrix(br,[list(eq) for eq in F.equations()]) #list of equations defining F
        eqns=eqns.stack(bB) #augmented matrix including fixed point subspace equations
        #get list of inequalities for fundamental domain as convex polyhedron
        ieqs=matrix(br,[list(ieq) for ieq in F.inequalities()])
        #form intersection of fund_domain with fixed point subspace
        intersection_F_fixed_pt_subspace=Polyhedron(ieqs=ieqs,eqns=eqns)
        verts=intersection_F_fixed_pt_subspace.vertices()
        for vert in verts:
            new_verts.add(tuple(vert))
    return sorted(list(new_verts))

#action of group element g on a face
#if g.face is not in fund_domain defined by vertex list, return None
def act_face(g,face,vertices):
    try:
        return [vertices.index(tuple(act_vect(g,vertices[i]))) for i in face]
    except ValueError:
        return None
    
#determine if face1 is glued to face2 by the action of G
#faces are given by a list vertex indices
def faces_glued(face1,face2,G,vertices):
    for g in G:
        g_face1 = act_face(g,face1,vertices)
        if g_face1 is not None and set(g_face1) == set(face2):
            return True
    return False

#for two glued faces face1, face2, determine if the gluing map preserves orientation
#if a gluing map g1 preserves orientation, then all do since g1*g2^(-1) is a self map and must be trivial
def gluing_preserves_orientation(face1,face2,vertices):
    S=SymmetricGroup(range(len(vertices))) #symmetric group on vertex indices
    for g in G:
        g_face1 = act_face(g,face1,vertices)
        if g_face1 is not None and set(g_face1) == set(face2):
            sigma_g=perm_from_sort(g_face1)
            g_cycle=PermutationGroupElement(sigma_g,parent=S) #get permutation which orders g.v_{i_0}...g.v_{i_k}
            return g_cycle.sign()
    return None

#determine if there is group element which glues face non-trivially to itself
def trivial_self_gluing(face,G,vertices):
    trivial = True
    for g in G:
        #act on vertices defining face
        g_face = act_face(g,face,vertices)
        if g_face is not None and set(g_face) == set(face):
            #check if the gluing is non-trivial
            if [g_face.index(i) for i in face] != [i for i in range(len(face))]:
                trivial = False
    return trivial

#toss out any facets with non-trivial self-gluings such as [0,3,4] or [1,6]
#build list of faces for each dim
import itertools
def faces(k): return [face for face in itertools.combinations(range(len(vertices)), k+1) if trivial_self_gluing(face,G,vertices)]

#find classes of vertices which are glued
def glued_faces(k):
    glued_verts = []
    for i in faces(k):
        found = False
        for equiv_class in glued_verts:
            if len(equiv_class) >= 1:
                if faces_glued(i,equiv_class[0],G,vertices):
                    equiv_class.append(i)
                    found=True
        if not found:
            glued_verts.append([i])
    return glued_verts

#compute boundary of each face. keep track of orientation
def boundary(face,glued_faces,vertices):
    dim = len(face)-1 #face/simplex dimension is number of vertices in face-1
    boundary = [0]*len(glued_faces[dim-1]) #initialize vector to count occurence boundary faces
    for i in range(len(face)):
        face_remove_i = face[:i] + face[i+1:] #remove vertex at index i
        sign = (-1)^i
        #find representative to which face\{i} is glued to
        for glued_face in glued_faces[dim-1]:
            face_rep = glued_face[0] #use first face in glued_face list as representative
            face_index = None
            if face_remove_i in glued_face: #check if face_remove_i is in gluing class
                face_index = glued_faces[dim-1].index(glued_face) #get index of face in glued (k-1)-faces
                orient_preserve = gluing_preserves_orientation(face_rep,face_remove_i,vertices) #determine if gluing map preserves orientation
            if face_index is not None:
                boundary[face_index] += orient_preserve*sign #add up using index
    return boundary

#compute the boundary map matrix over \ZZ
def boundary_map(k,glued_face_list,vertices): return matrix(ZZ,[boundary(face[0],glued_face_list,vertices) for face in glued_face_list[k]]).transpose()

#define a chain complex, optionally truncated at max_degree
def chain_complex(max_degree=1): return ChainComplex({k:boundary_map(k,glued_face_list,vertices) for k in range(1,max_degree+1)},degree=-1)

#define the finite group
G = PermutationGroup(embed_gens(n)) #symmetric group \sigma_n as a subgroup of \sigma_N, with 2 generators
#G = MatrixGroup(matrix(QQ,[[0,0,0,1],[1,0,0,0],[0,1,0,0],[0,0,1,0]]))
#G=SymmetricGroup(4)

#compute fundamental domain to get list of vertices
F=fund_domain(group=G,center=[1,2,3,5,6,8]); F[1]

#compute a cross section of the fundamental domain (sum of components is constant)
F2=cross_section(F[1])
[tuple(v) for v in F2.vertices()]

# compute the fixed points under G in polyhedron F to include as new vertices
vertices=fixed_verts(F2,QQ); vertices

#compute glued 0-faces
zero_faces=glued_faces(0); zero_faces

#compute glued 1-faces
one_faces=glued_faces(1); one_faces

two_faces=glued_faces(2); two_faces

# takes far too long on a laptop to computer three_faces
#three_faces=glued_faces(3); three_faces

# put glued faces together in list indexed by dimension
glued_face_list=[zero_faces,one_faces,two_faces]

# compute the boundary map d_1: C_1 -> C_0
d1=boundary_map(1, glued_face_list, vertices); d1

# compute boundary map d_2: C_2 -> C_1
d2=boundary_map(2, glued_face_list, vertices); d2

#d3 takes way too long to compute on a laptop
#d3=boundary_map(3, glued_face_list, vertices); d3

#form the chain complex with truncation at degree 2
cc = chain_complex(max_degree=2)

#compute the homology of the chain complex
cc.homology()

list(map(lambda x: len(x), glued_face_list))