import find_mk_params, sys, os, re, cjson
mkparams = []  # A list of tuples (Ds, Pn, Dn, Ps)
mkout_dir = find_mk_params.run_name+'/mktest_out/'
for f in os.listdir(mkout_dir):
	f = mkout_dir + f
	fixedAS = re.findall('#Fixed\s+(\d+)\s+(\d+)', open(f).read())
	polyAS  = re.findall('#Poly\s+(\d+)\s+(\d+)', open(f).read())
	if len(fixedAS) != 1:
		sys.stderr.write("File %s has %i #Fixed lines\n"%(f,len(fixedAS)))
		continue
	if len(polyAS) != 1:
		sys.stderr.write("File %s has %i #Poly lines\n"%(f,len(polyAS)))
		continue
	Dn, Ds = fixedAS[0]
	Pn, Ps = polyAS[0]
	mkparams.append(map(int,(Ds,Pn,Dn,Ps)))
mkparams = filter(lambda (Ds,Pn,Dn,Ps): Dn*Ps != 0, mkparams)
mkparams = filter(lambda (Ds,Pn,Dn,Ps): Ds>10 and Pn>10 and Dn>10 and Ps>10, mkparams)
avalues = map(lambda (Ds,Pn,Dn,Ps): 1-(Ds*Pn*1.0)/(Dn*Ps), mkparams)
outdir = {'mkparams' : mkparams, 'avalues' : avalues}
open(find_mk_params.run_name+'/mkout.cjson','w').write(cjson.encode(outdir))
open(find_mk_params.run_name+'/alpha_values.txt','w').write(' '.join(str(av) for av in avalues))

