#computes the sample frechet mean of two points using the formula from [JW]
def sample_frechet_mean(s1,s2):
	#compute a fixed reference fundamental domain
	F_l=fund_domain()[1]
	
	#compute orbits of all s_i
	orb1 = list(set([tuple(act_vect(g,vector(s1))) for g in G]))
	orb2 = list(set([tuple(act_vect(g,vector(s2))) for g in G]))
	
	orb1 = [vector(o1) for o1 in orb1]
	orb2 = [vector(o2) for o2 in orb2]
	
	q1=len(orb1)
	q2=len(orb2)
	
	#compute domains for each orbit representative
	domains1 = [fund_domain(o1)[1] for o1 in orb1]
	domains2 = [fund_domain(o2)[1] for o2 in orb2]
	
	avgs = []
	
	#compute orbit averages, ensuring they are contained
	#in the appropriate region R_i1,...,ir as per [JW]
	for i in range(q1):
		for j in range(q2):
			R = domains1[i].intersection(domains2[j])
			avg=(1/2)*(orb1[i]+orb2[j])
			if R.contains(avg): avgs.append(avg)
	
	max_norm = max([av.norm() for av in avgs])
	maxes=[]
	for av in avgs: 
		if av.norm() == max_norm: maxes.append(av)
	
	#find a fund. domain representative for each Frechet mean
	domain_reps=[tuple(fund_domain_rep(m,F_l)) for m in maxes]
	
	return [vector(s) for s in set(domain_reps)]
	