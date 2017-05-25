def proc_dist(s1,s2): return min([(vector(act_vect(g,vector(s2)))-vector(s1)).norm() for g in G])

#computes the sample frechet mean of two points using the formula from [JW]
def sample_frechet_mean(s1,s2):
	F_l=fund_domain()[1]
	
	orb1 = list(set([tuple(act_vect(g,vector(s1))) for g in G]))
	orb2 = list(set([tuple(act_vect(g,vector(s2))) for g in G]))
	
	orb1 = [vector(o) for o in orb1]
	orb2 = [vector(o) for o in orb2]
	
	q1=len(orb1)
	q2=len(orb2)
	
	domains1 = [fund_domain(o1)[1] for o1 in orb1]
	domains2 = [fund_domain(o2)[1] for o2 in orb2]
	
	avgs = []
	
	for i in range(q1):
		for j in range(q2):
			R = domains1[i].intersection(domains2[j])
			avg=(1/2)*factorial(n)*(orb1[i]/q1+orb2[j]/q2)
			if R.contains(avg): avgs.append(avg)
	
	max_norm = max([av.norm() for av in avgs])
	maxes=[]
	for av in avgs: 
		if av.norm() == max_norm: maxes.append(av)
	
	domain_reps=[tuple(fund_domain_rep(m,F_l)) for m in maxes]
	
	return [vector(s) for s in set(domain_reps)]
	