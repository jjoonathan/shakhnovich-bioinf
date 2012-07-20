import sys, os, re, sqlite3

import eventlet.greenpool
import eventlet.green.subprocess as gsub

run_names = []
for s in sys.argv[1:]:
	if re.match(r'[a-zA-Z]+(-[a-zA-Z]+)+$',s):
		run_names.append(s)

for run_name in run_names:
	rows = []
	is_header = True
	for line in open(run_name+'/ortho_ids'):
		if is_header:
			is_header = False
			continue
		fields = line.split('\t')
		rows.append(fields)
	msg = "Running MKtest for %s... "%run_name
	sys.stdout.write(msg)
	sys.stdout.flush()
	if not os.path.isdir(run_name+'/MKout'):
		os.mkdir(run_name+'/MKout')
	def mktest_row(row):
		ortho_id = row[0]
		num_lspecies = row[1]
		ifname = '%s/clustalout/%s.fasta'%(run_name,ortho_id)
		ofname = '%s/MKout/%s'%(run_name,ortho_id)
		if os.path.isfile(ofname):
			return
		of=open(ofname,'w')
		errf=open(ofname+'_err','w')
		gsub.call(('MKtest','-n',num_lspecies,'-i',ifname),stdout=of,stderr=errf)
		of.close()
		errf.close()
	mktest_pool = eventlet.greenpool.GreenPool(size=16)
	i=0
	for result in mktest_pool.imap(mktest_row, rows):
		i+=1
		sys.stdout.write('\x1B[%iG%i%% (%i/%i)'%(len(msg)+1,i*100./len(rows),i,len(rows)))
		sys.stdout.flush()
	print "Done."
	msg = "Creating alpha map for %s... "%run_name
	sys.stdout.write(msg)
	sys.stdout.flush()
	mkout_dir = run_name+'/MKout/'
	outf=open(run_name+'/alpha_values.tsv','w')
	outf.write("GroupID\tLCount\tRCount\tAlpha\tDn\tDs\tPn\tPs\tLShortnames\tRShortnames\n")
	for row in rows:
		lshortnames, rshortnames = row[3], row[4]
		lcount, rcount = row[1], row[2]
		ortho_id = row[0]
		mkout_fname = mkout_dir+ortho_id
		mktest_results = open(mkout_fname).read()
		fixedAS = re.findall('#Fixed\s+(\d+)\s+(\d+)', mktest_results)
		polyAS  = re.findall('#Poly\s+(\d+)\s+(\d+)', mktest_results)
		if len(fixedAS) != 1:
			sys.stderr.write("Group %s has %i #Fixed lines\n"%(ortho_id,len(fixedAS)))
			continue
		if len(polyAS) != 1:
			sys.stderr.write("Group %s has %i #Poly lines\n"%(ortho_id,len(polyAS)))
			continue
		Dn, Ds = fixedAS[0]
		Pn, Ps = polyAS[0]
		Dn, Ds, Pn, Ps = map(int, (Dn,Ds,Pn,Ps))
		if Dn*Ps != 0:
			a = 1-(Ds*Pn*1.0)/(Dn*Ps)
		else:
			a = ''
		outf.write("{ortho_id}\t{lcount}\t{rcount}\t{a}\t{Dn}\t{Ds}\t{Pn}\t{Ps}\t{lshortnames}\t{rshortnames}\n".format(**locals()))
	outf.close()


