import find_mk_params, sys, os, re, cjson
mkout_dir = find_mk_params.run_name+'/mktest_out/'
gm = cjson.decode(open(find_mk_params.run_name+'/genemaps.json').read())
sn = cjson.decode(open(find_mk_params.run_name+'/species_names.json').read())
geneFamilies = gm['geneFamilies']
geneToSpecies= gm['geneToSpecies']
outf=open(find_mk_params.run_name+'/alpha_values.txt','w')
outf.write("K12Gene\tGene\tAlpha\tDn\tDs\tPn\tPs\n")
for gf in geneFamilies:
	gf.sort()
	outfname = mkout_dir+str(repr(gf).__hash__())+'.txt'
	try:
		outfd = open(outfname)
		mktest_results = outfd.read()
	except Exception as e:
		sys.stderr.write("Couldn't read output for %s (%s)\n"%(outfname,e))
		continue
	fixedAS = re.findall('#Fixed\s+(\d+)\s+(\d+)', mktest_results)
	polyAS  = re.findall('#Poly\s+(\d+)\s+(\d+)', mktest_results)
	if len(fixedAS) != 1:
		sys.stderr.write("File %s has %i #Fixed lines\n"%(outfname,len(fixedAS)))
		continue
	if len(polyAS) != 1:
		sys.stderr.write("File %s has %i #Poly lines\n"%(outfname,len(polyAS)))
		continue
	Dn, Ds = fixedAS[0]
	Pn, Ps = polyAS[0]
	Dn, Ds, Pn, Ps = map(int, (Dn,Ds,Pn,Ps))
	k12gene = ''
	for g in gf:
		if geneToSpecies[g] == 'Escherichia coli K-12':
			k12gene = g
			break
	gene = k12gene if k12gene!='' else gf[0]
	if Dn*Ps == 0:
		continue
	a = 1-(Ds*Pn*1.0)/(Dn*Ps)
	outf.write("%s\t%s\t%f\t%i\t%i\t%i\t%i\n"%(k12gene,gene,a,Dn,Ds,Pn,Ps))
outf.close()

