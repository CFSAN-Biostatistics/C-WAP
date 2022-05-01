
import cvxpy
import argparse
import pandas
import numpy
import ray
import os
import warnings

os.environ["RAY_DISABLE_IMPORT_WARNING"]="1"

rng = numpy.random.default_rng()
parser = argparse.ArgumentParser()
parser.add_argument('-m','--markers', help='variant groups marker table', required=True)
parser.add_argument('-p','--pools', help='pools variants table', required=True)
parser.add_argument('-o','--output', help='decomposed samples matrix output', required=True)
parser.add_argument('-s','--status', help='decomposed samples status output', required=True)
parser.add_argument('--min-marker-freq', help='minimum freq for a marker variant', type=float, default=0.8)
parser.add_argument('--max-noise-freq', help='marker variants below this freq will be considered zero', type=float, default=0.1)
parser.add_argument('--solver', help='optimization solver', type=str, default='ECOS', choices=['ECOS','SCS'])
parser.add_argument('--objective', help='optimization objective', type=str, default='NLL', choices=['LAE','NLL'])
parser.add_argument('--nboot', help='number of bootstrap replications', type=int, default=100)
parser.add_argument('--threads', help='number of threads', type=int, default=1)
args = parser.parse_args()

ray.init(num_cpus=args.threads, local_mode=False, num_gpus=0, logging_level=30)

assert args.min_marker_freq > args.max_noise_freq

eps = numpy.nextafter(0,1) # smallest positive float

def negLLH(f, vars_snps, pm, pw):
	ps = cvxpy.matmul(vars_snps, f)
	ps += eps # avoid Log(0)
	nllh = (-1)*cvxpy.sum(cvxpy.multiply(pm,cvxpy.log(ps))+cvxpy.multiply(pw,cvxpy.log(1-ps)))
	return nllh

@ray.remote
def solve(vs, ps_m, ps_dp, n_vars):
	warnings.filterwarnings(action='ignore', message=r'.*Solution may be inaccurate. Try another solver.*')
	warnings.filterwarnings(action='ignore', message=r'.*invalid value encountered in log.*')
	x = cvxpy.Variable(n_vars,nonneg=True)
	if args.objective == 'NLL':
		# negative log likehood
		objective = cvxpy.Minimize(negLLH(x, vs, ps_m, ps_dp - ps_m))
	elif args.objective == 'LAE':
		# least absolute error
		ps = ps_m/ps_dp
		objective = cvxpy.Minimize(cvxpy.sum(cvxpy.abs(vs@x - ps)))
	
	constraints = [0 <= x, x <= 1, sum(x) <= 1]
	prob = cvxpy.Problem(objective, constraints)
	try:
		result = prob.solve(verbose=False, solver=args.solver, warm_start=True)
	except cvxpy.error.SolverError:
		# re-try with SCS
		if args.solver != "SCS":
			try:
				result = prob.solve(verbose=False, solver="SCS", warm_start=True)
			except cvxpy.error.SolverError:
				# error again
				result = None
		else:
			result = None
	return (prob, result, x)

def gen_bootstrap_samples(vec_ref, vec_alt, n_replicates):
	total = sum(vec_ref) + sum(vec_alt)
	probs = numpy.concatenate((vec_ref, vec_alt))/total
	bts = rng.multinomial(n=total, pvals=probs, size=n_replicates)
	return numpy.split(bts, [len(vec_ref)], axis=1)

@ray.remote
def solve_bootstrap(vs, ps_m, ps_dp, n_vars, n_replicates):
	bts_ref, bts_alt = gen_bootstrap_samples(ps_dp-ps_m, ps_m, n_replicates)
	preds = []
	works = []
	for i in range(n_replicates):
		wid = solve.remote(vs, bts_alt[i], bts_ref[i]+bts_alt[i], n_vars)
		works.append(wid)

	for wid in works:
		prob, result, x = ray.get(wid)
		preds.append(x.value)

	preds = numpy.stack(preds)
	means = preds.mean(axis=0)
	std_err = preds.std(axis=0)
	return means, std_err

def solve_samples(markers_df, pools, marker_snps_list, tags, pool_tags):
	# markers matrix
	vars_snps = markers_df[markers_df.mut.isin(marker_snps_list)].pivot_table(index='mut', columns='sample', values='af', fill_value=0)
	assert (vars_snps.columns == tags).all()
	assert (vars_snps.index == marker_snps_list).all()
	vars_snps = vars_snps.to_numpy()
	
	# pools matrix
	pools_snps_m = pools[pools.mut.isin(marker_snps_list)].pivot_table(index='mut', columns='sample', values='adalt').to_numpy()
	pools_snps_dp = pools[pools.mut.isin(marker_snps_list)].pivot_table(index='mut', columns='sample', values='dp').to_numpy()
	pools_snps_af = pools[pools.mut.isin(marker_snps_list)].pivot_table(index='mut', columns='sample', values='af').to_numpy()
	
	status = open(args.status, 'wt')
	print("\t".join(['sample','status','fit_error','solve_time','n_markers']), file=status)
	out = open(args.output, 'wt')
	
	if args.nboot > 0:
		print("\t".join(['sample', 'variant_group', 'proportion', 'mean', 'std_error']), file=out)
	else:
		print("\t".join(['sample', 'variant_group', 'proportion']), file=out)
	
	works = {}
	bts_works = {}
	for i,s in enumerate(pool_tags):
		zero_idxs = numpy.argwhere(numpy.logical_or(pools_snps_dp[:,i]==0, numpy.isnan(pools_snps_dp[:,i])))
		ps_dp = numpy.delete(pools_snps_dp[:,i], zero_idxs)
		ps_m = numpy.delete(pools_snps_m[:,i], zero_idxs)
		vs = numpy.delete(vars_snps, zero_idxs, axis=0)
		
		wid = solve.remote(vs, ps_m, ps_dp, len(tags))
		works[i] = wid
		
		if args.nboot > 0:
			bts_wid = solve_bootstrap.remote(vs, ps_m, ps_dp, len(tags), args.nboot)
			bts_works[i] = bts_wid

	for i, wid in works.items():
		s = pool_tags[i]
		prob, result, x = ray.get(wid)

		if result is None:
			print("\t".join([s, 'error', "NA", "NA", str(len(vs))]), file=status)
			for j in range(len(tags)):
				if args.nboot > 0:
					print("\t".join([s, tags[j], "NA", "NA", "NA"]), file=out)
				else:
					print("\t".join([s, tags[j], "NA"]), file=out)
			print('warning: pool '+s+' failed to solve')
		else:
			print("\t".join([ s, prob.status, str(result), str(prob.solver_stats.solve_time), str(len(vs))]), file=status)
			
			if args.nboot > 0:
				means, std_err = ray.get(bts_works[i])
			
			for j in range(len(tags)):
				if args.nboot > 0:
					print("\t".join([s, tags[j], str(max(0,x.value[j])), str(means[j]), str(std_err[j]) ]), file=out)
				else:
					print("\t".join([s, tags[j], str(max(0,x.value[j])) ]), file=out)

def main():
	# markers data
	markers_df = pandas.read_csv(args.markers, sep='\t')
	markers_df['mut'] = markers_df.apply(lambda r: str(r['pos'])+":"+r['ref']+">"+r["alt"], axis=1)
	markers_df['af'] = markers_df['adalt']/markers_df['dp']
	markers_df.loc[markers_df.af<args.max_noise_freq, 'af'] = 0
	marker_snps_list = sorted(markers_df[markers_df.af>=args.min_marker_freq].mut.unique())
	tags = sorted(markers_df['sample'].unique())
	
	# pool data
	pools = pandas.read_csv(args.pools, sep='\t')
	pool_tags = sorted(pools['sample'].unique())
	pools['mut'] = pools.apply(lambda r: str(r['pos'])+":"+r['ref']+">"+r["alt"], axis=1)
	pools['af'] = pools['adalt']/pools['dp']
	
	# check for samples with zero coverage
	cov = pools.groupby('sample').agg('sum')
	zero_cov_samples = cov[cov.dp==0].index
	if len(zero_cov_samples) > 0:
		print('warning: dropping samples with zero coverage: '+str(zero_cov_samples))
		pools = pools[ ~( pools['sample'].isin(zero_cov_samples) ) ]
		pool_tags = sorted(pools['sample'].unique())
	
	# check missing muts
	missing_muts = set(marker_snps_list) - set(pools.mut)
	if len(missing_muts) > 0:
		print('warning: dropping '+str(len(missing_muts))+' markers because they are not present on the pools')
		marker_snps_list = sorted(list(set(marker_snps_list) - missing_muts))
	
	print('using '+str(len(marker_snps_list))+' markers to estimate '+str(len(tags))+' variant groups')
	
	solve_samples(markers_df, pools, marker_snps_list, tags, pool_tags)

if __name__ == "__main__":
	main()
