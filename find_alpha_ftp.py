print "Loading libraries..."
import sys, os, lxml.html, xml.sax, sqlite3, itertools, subprocess, random, cookielib
import urllib, os.path, mmap, zlib, re, xml.etree.cElementTree, code

import cjson
import eventlet.greenpool
import eventlet.green.urllib2 as urllib2
import eventlet.green.subprocess as gsub

class AlphaFinder():
	def __init__(self):
		pass
	def run(self):
		# self.conn = MySQLdb.Connection('localhost')
		self.errcount = 0
		if not os.path.isdir(self.run_name):
			os.mkdir(self.run_name)
		self.conn = sqlite3.Connection('data.db', isolation_level=None)
		self.conn.text_factory = str
		self.cursor = self.conn.cursor()
		# self.cursor.execute('USE a;')
		self.update_roundup_species_list()
		assert len(self.l_taxon_ids())
		assert len(self.r_taxon_ids())
		self.create_runlist()
		self.fetch_roundup_orthology_assignment()
		self.download_full_genomes()
		self.load_genomes_into_db()
		self.load_orthology_assignment_into_db()
		self.create_indexes()
		self.output_clustalin()
	def update_roundup_species_list(self, force=False):
		sys.stdout.write("Updating roundup species...")
		sys.stdout.flush()
		if self.table_exists('roundup_species') and force != True:
			print " Done. Used cache in data.db."
			return
		self.cursor.execute('BEGIN')
		self.cursor.execute("""
			CREATE TABLE roundup_species (
				ncbi_taxon INTEGER PRIMARY KEY,
				name TEXT,
				num_genes INTEGER,
				proteome_in_db INTEGER DEFAULT 0,
				blacklisted INTEGER DEFAULT 0
			);""")
		gene_doc = urllib2.urlopen('http://roundup.hms.harvard.edu/genomes/').read()
		gene_tree = lxml.html.document_fromstring(gene_doc)
		taxon_rows = gene_tree.xpath('//table[@id="genome_descs"]/tr')[1:]
		# <tr>   <th>TCBI Taxon</th>  <th>Name</th>  <th>Category</th>  <th>Number of Genes</th>   </th>
		for row in taxon_rows:
			children = row.getchildren()
			ncbi_taxon = int(children[0].text_content())
			name = children[1].text_content()
			gene_count = int(children[3].text_content())
			self.cursor.execute('INSERT INTO roundup_species (ncbi_taxon, name, num_genes) VALUES (?,?,?)', (ncbi_taxon, name, gene_count))
		self.cursor.execute('COMMIT')
		print " Done."
	def create_runlist(self):
		self.cursor.execute('BEGIN')
		self.cursor.execute("""
			CREATE TABLE IF NOT EXISTS runlist (
				run_name TEXT,
				ncbi_taxon INTEGER,
				is_rspecies INTEGER,
				UNIQUE(run_name, ncbi_taxon)
			);""")
		for taxonid in self.l_taxon_ids():
			self.cursor.execute("""
				INSERT OR IGNORE INTO runlist (run_name, ncbi_taxon, is_rspecies)
				VALUES (?,?,?)""",            (self.run_name, taxonid, 0))
		for taxonid in self.r_taxon_ids():
			self.cursor.execute("""
				INSERT OR IGNORE INTO runlist (run_name, ncbi_taxon, is_rspecies)
				VALUES (?,?,?)""",            (self.run_name, taxonid, 1))
		self.cursor.execute('COMMIT')
	def download_full_genomes(self):
		failed_taxons = []
		sys.stdout.write("Downloading full genomes...\n")
		sys.stdout.flush()
		for ncbi_taxon in itertools.chain(self.l_taxon_ids(), self.r_taxon_ids()):
			cachefile_name = 'full_proteomes/%i.xml'%(ncbi_taxon)
			if os.path.isfile(cachefile_name):
				continue
			url = 'http://www.uniprot.org/uniprot/?query=taxonomy%%3A%i+AND+keyword%%3A181&format=xml'%(ncbi_taxon)
			# url = 'http://www.uniprot.org/uniprot/?query=organism%%3a%i+keyword%%3a181&format=xml'%(ncbi_taxon,)
			args = ('curl', '--create-dirs', '-o', cachefile_name+'_tmp', url)
			print "\n"+' '.join(args)+"\n"+str(ncbi_taxon)+"\n"
			subproc = subprocess.Popen(args, stdout=sys.stdout, stderr=sys.stderr)
			retval = subproc.wait()
			if retval == 0 and os.path.isfile(cachefile_name+'_tmp'):
				os.rename(cachefile_name+'_tmp', cachefile_name)
			else:
				failed_taxons.append(ncbi_taxon)
		if failed_taxons:
			print "Done fetching genomes (failed taxon ids:%s)"%str(failed_taxons)
		print "Done fetching full genomes."
	def table_exists(self, tname):
		self.cursor.execute('SELECT COUNT(name) FROM sqlite_master WHERE type="table" AND name=?',[tname])
		#self.cursor.execute('SELECT COUNT(*) FROM information_schema.tables WHERE table_schema="a" AND table_name=?',[tname])
		return bool(self.cursor.fetchone()[0])
	def fetch_roundup_orthology_assignment(self):
		sys.stdout.write("Fetching roundup orthology assignments...")
		sys.stdout.flush()
		if os.path.isfile('%s/roundup.xml'%(self.run_name)):
			print " Done. Used %s/roundup.xml."%self.run_name
			return True
		# First grab the CSRF, cookies with a GET request
		cjar = cookielib.CookieJar()
		opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cjar))
		request = urllib2.Request('http://roundup.hms.harvard.edu/retrieve/')
		response = opener.open(request).read()
		rtree = lxml.html.document_fromstring(response)
		key = rtree.xpath('//input[@name="csrfmiddlewaretoken"]/@value')[0]
		taxons = self.l_taxon_names()
		taxons.extend(self.r_taxon_names())
		taxons.extend('')
		# Now use the session ID to submit the form
		params = []
		params.append(('csrfmiddlewaretoken',key))
		params.append(('genomes_filter','E'))
		params.append(('genomes_filter','B'))
		params.append(('genomes_filter','A'))
		params.append(('genomes_filter','V'))
		params.append(('genome_choices','Escherichia coli 536'))
		params.append(('genomes',"\r\n".join(taxons)))
		params.append(('divergence','0.8'))
		params.append(('evalue','1e-5'))
		params.append(('distance_lower_limit',''))
		params.append(('distance_upper_limit',''))
		# params.append(('tc_only','on'))
		loc = 'http://roundup.hms.harvard.edu/retrieve/'
		dat = urllib.urlencode(params)
		hdrs = {}
		hdrs['User-Agent'] = 'Mozilla/5.0 (X11; U; Linux i686) Gecko/20071127 Firefox/2.0.0.11'
		hdrs['Referrer'] = 'http://roundup.hms.harvard.edu/retrieve/'
		req = urllib2.Request('http://roundup.hms.harvard.edu/retrieve/', data=dat, headers=hdrs)
		response = opener.open(req)
		response_str = response.read()
		open('/tmp/roundup_intermediate_response.html','w').write(response_str)
		rtree = lxml.html.document_fromstring(response_str)
		links = rtree.xpath('//a/@href')
		xmllinks = filter(lambda s: s.find('tt=xml')!=-1, links)
		if len(xmllinks)!=1:
			print " ERROR: bad number of XML links (%i xmllinks of %i links)"%(len(xmllinks), len(links))
			open('/tmp/roundup.html','w').write(response_str)
			exit(1)
		# Fetch the XML file of results
		req = urllib2.Request('http://roundup.hms.harvard.edu'+xmllinks[0])
		response = opener.open(req).read()
		outfile = open('%s/roundup_tmp.xml'%(self.run_name),'w')
		outfile.write(response)
		outfile.close()
		os.rename('%s/roundup_tmp.xml'%(self.run_name), '%s/roundup.xml'%(self.run_name))
		print " Done."
		return True
	def load_orthology_assignment_into_db(self):
		if self.table_exists("roundup_orthologGroups"):
			self.cursor.execute('SELECT COUNT(*) FROM roundup_orthologGroups WHERE run_name=? LIMIT 1',[self.run_name])
			if self.cursor.fetchone()[0]:
				print "Loading orthology assignments... Done, using existing."
				return
		print "Loading orthology assignments..."
		self.cursor.execute('BEGIN')
		self.cursor.execute("""
			CREATE TABLE IF NOT EXISTS roundup_orthologGroups (
				ortho_id INTEGER PRIMARY KEY,
				avgdist FLOAT,
				run_name TEXT
			);""")
		self.cursor.execute("""
			CREATE TABLE IF NOT EXISTS roundup_orthologAssignments (
				ortho_id INTEGER,
				accession TEXT,
				ncbi_taxon INTEGER
			);""")
		roundup_str = open('%s/roundup.xml'%(self.run_name)).read()
		roundup_tree = xml.etree.cElementTree.fromstring(roundup_str)
		print "Parsng ortholog assignments..."
		geneIdToAccession = {}
		for species in roundup_tree.findall('{http://orthoXML.org/2011/}species'):
			ncbi_taxon = int(species.get('NCBITaxId'))
			for gene in species.findall('.//{http://orthoXML.org/2011/}gene'):
				geneIdToAccession[gene.get('id')] = (ncbi_taxon, gene.get('protId'))
		ortho_id = self.cursor.execute('SELECT MAX(ortho_id) FROM roundup_orthologGroups').fetchone()[0]
		if ortho_id == None:
			ortho_id = 0
		insert_orthologAssns = []
		insert_orthologGroups = []
		for orthologGroup in roundup_tree.findall('.//{http://orthoXML.org/2011/}orthologGroup'):
			avgdist = 0
			ortho_id += 1
			for score in orthologGroup.findall('{http://orthoXML.org/2011/}score'):
				if score.get('id') == 'avgdist':
					avgdist = float(score.get('value'))
			insert_orthologGroups.append((ortho_id, avgdist, self.run_name))
			for geneRef in orthologGroup.findall('{http://orthoXML.org/2011/}geneRef'):
				ncbi_taxon, accession = geneIdToAccession[geneRef.get('id')]
				insert_orthologAssns.append((ortho_id, accession, ncbi_taxon))
		sys.stdout.write("Inserting %i ortholog assignments, %i groups... "%(len(insert_orthologAssns),len(insert_orthologGroups)))
		sys.stdout.flush()
		self.cursor.executemany('INSERT OR IGNORE INTO roundup_orthologGroups (ortho_id, avgdist, run_name) VALUES (?,?,?)', insert_orthologGroups)
		self.cursor.executemany('INSERT OR IGNORE INTO roundup_orthologAssignments (ortho_id, accession, ncbi_taxon) VALUES (?,?,?)', insert_orthologAssns)
		self.cursor.execute('COMMIT')
		print "Done."
	def create_tables(self):
		self.cursor.execute('BEGIN')
		self.cursor.execute("""
			CREATE TABLE IF NOT EXISTS uniprot_accession (
				accession VARCHAR(10),
				shortname VARCHAR(255),
				ncbi_taxon INTEGER
			);""")
		self.cursor.execute("""
			CREATE TABLE IF NOT EXISTS uniprot_dbReference (
				shortname VARCHAR(255),
				ncbi_taxon INTEGER,
				type VARCHAR(255),
				id VARCHAR(255),
				aa_seq_id TEXT DEFAULT NULL,
				nuc_seq_id TEXT DEFAULT NULL,
				ncbi_gene_gi INTEGER
			);""")
		self.cursor.execute("""
			CREATE TABLE IF NOT EXISTS ncbi_seq (
				shortname VARCHAR(255),
				ncbi_taxon INTEGER,
				fasta TEXT
			);""")
		self.cursor.execute("""
			CREATE TABLE IF NOT EXISTS gene2accession (
				protein_accession TEXT PRIMARY KEY,
				protein_gi INTEGER,
				nucleotide_accession TEXT,
				nucleotide_gi INTEGER
			);""")
		self.cursor.execute('COMMIT')
	def create_indexes(self):
			self.cursor.execute("CREATE INDEX IF NOT EXISTS uniprot_accession_idx0 ON uniprot_accession (accession);")
			self.cursor.execute("CREATE INDEX IF NOT EXISTS uniprot_accession_idx1 ON uniprot_accession (shortname,ncbi_taxon);")
			self.cursor.execute("CREATE INDEX IF NOT EXISTS uniprot_dbReference_idx0 ON uniprot_dbReference (shortname,ncbi_taxon);")
			self.cursor.execute("CREATE INDEX IF NOT EXISTS uniprot_dbReference_idx1 ON uniprot_dbReference (type);")
			self.cursor.execute("CREATE INDEX IF NOT EXISTS uniprot_dbReference_idx3 ON uniprot_dbReference (id);")
			self.cursor.execute("CREATE INDEX IF NOT EXISTS ncbi_seq_idx0 ON ncbi_seq (shortname,ncbi_taxon);")
			self.cursor.execute("CREATE INDEX IF NOT EXISTS roundup_orthologAssignments_idx0 ON roundup_orthologAssignments (ortho_id);")
			self.cursor.execute("CREATE INDEX IF NOT EXISTS runlist_idx0 ON runlist (ncbi_taxon);")
	def load_genomes_into_db(self):
		msg = "Loading genomes into db...\n"
		sys.stdout.write(msg)
		sys.stdout.flush()
		taxons = self.l_taxon_ids()
		taxons.extend(self.r_taxon_ids())
		self.create_tables()
		for ncbi_taxon in taxons:
			self.cursor.execute('SELECT proteome_in_db FROM roundup_species WHERE ncbi_taxon=?',[ncbi_taxon])
			if self.cursor.fetchone()[0] == 1:
				continue
			file_name = 'full_proteomes/%i.xml'%(ncbi_taxon)
			with open(file_name,'r+b') as f:
				print "\tProcessing: %s"%file_name
				fmap = mmap.mmap(f.fileno(), 0)
				insert_uniprot = []
				insert_accession = []
				insert_dbref = []
				for matchobj in re.finditer(r"(<entry.{100,}?</entry>)", fmap, flags=re.MULTILINE+re.DOTALL):
					entry_str = matchobj.group(0)
					tree = xml.etree.cElementTree.fromstring(entry_str)
					shortname = tree.findall('name')[0].text
					for accession in tree.findall('accession'):
						insert_accession.append((accession.text,shortname,ncbi_taxon))
					for dbref in tree.findall('dbReference'):
						db_type, db_id, aa_seq_id, nuc_seq_id = dbref.get('type'), dbref.get('id'), None, None
						for prop in dbref.findall('property'):
							if prop.get('type') == 'protein sequence ID':
								aa_seq_id = prop.get('value')
							elif prop.get('type') == 'nucleotide sequence ID':
								nuc_seq_id = prop.get('value')
						insert_dbref.append((shortname, ncbi_taxon, db_type, db_id, aa_seq_id, nuc_seq_id))
				self.cursor.execute('BEGIN')
				self.cursor.executemany("""
					INSERT INTO uniprot_accession (accession, shortname, ncbi_taxon)
					VALUES (?,?,?)""", insert_accession)
				self.cursor.executemany("""
					INSERT INTO uniprot_dbReference (shortname, ncbi_taxon, type, id, aa_seq_id, nuc_seq_id)
					VALUES (?,?,?,?,?,?);""", insert_dbref)
				self.cursor.execute('COMMIT')
			self.cursor.execute('BEGIN')
			self.cursor.execute('UPDATE roundup_species SET proteome_in_db=1 WHERE ncbi_taxon=?',[ncbi_taxon])
			self.cursor.execute('COMMIT')
		print "Creating indices..."
	def load_gene2accession(self):
		if self.cursor.execute('SELECT COUNT(*) FROM gene2accession').fetchone()[0] > 100:
			return
		self.cursor.execute('BEGIN')
		valid_taxons = set(int(row[0]) for row in self.cursor.execute('SELECT ncbi_taxon FROM roundup_genomes').fetchall())
		badstart = "NO LINE STARTS WITH THIS"
		records = []
		for line in open('gene2accession'):
			if line.startswith(badstart):
				continue
			line = line.split('\t')
			taxid = int(line[0])
			if taxid not in valid_taxons:
				badstart = str(taxid)
				continue
			tax_id, GeneID, status, rna_accession, rna_gi, prot_accession,\
			prot_gi, nuc_accession, nuc_gi, nuc_start, nuc_end, nuc_orientation, assembly = line
			records.append((tax_id, GeneID, status, 
		self.cursor.execute('COMMIT')
	def output_clustalin(self):
		clustalin_tar = self.run_name+'.tar.bz2'
		if os.path.isfile(clustalin_tar):
			return
		self.output_clustalin_dir()
		self.output_bsub()
		print "\nPreparing clustalin tarfile..."
		clustalin_tar_f = open(clustalin_tar,'w')
		files = ['ortho_ids', 'out', 'clustalin', 'bsub.sh', 'bsub1.sh', 'clustalout']
		files = [self.run_name+'/'+s for s in files]
		options = ['tar','-cj']
		options.extend(files)
		subprocess.call(options,stdout=clustalin_tar_f)
		clustalin_tar_f.close()
	def output_clustalin_dir(self):
		msg = "Writing clustal input files... "
		sys.stdout.write(msg)
		sys.stdout.flush()
		clustalin_dir = self.run_name + '/clustalin'
		if not os.path.isdir(clustalin_dir):
			os.mkdir(clustalin_dir)
		ortho_grpids = [row[0] for row in self.cursor.execute("""
			SELECT DISTINCT oa.ortho_id
			FROM roundup_orthologAssignments AS oa
			INNER JOIN roundup_orthologGroups AS og ON (oa.ortho_id = og.ortho_id)
			WHERE og.run_name = ?""", [self.run_name]).fetchall()]
		self.num_ortho_grpids = len(ortho_grpids)
		ortho_id_list_file = open(self.run_name+'/ortho_ids','w')
		ortho_id_list_file.write("ortho_id\tl_count\tr_count\tlshortnames\trshortnames\n")
		try:
			existing_files = set(int(s.split('.')[0]) for s in os.listdir(clustalin_dir))
			ortho_grpids = set(ortho_grpids) - existing_files
		except Exception as e:
			pass
		i, num = 0, len(ortho_grpids)
		for ortho_grp_id in ortho_grpids:
			i += 1
			sys.stdout.write("\x1B[%iG%i%% (%i/%i)  "%(len(msg)+1,i*100./num,i,num))
			sys.stdout.flush()
			rseqs, lseqs = [], []
			rshortnames, lshortnames = [], []
			seqs = self.cursor.execute("""
				SELECT DISTINCT runlist.is_rspecies, ua.shortname, seq.nuc
				FROM roundup_orthologAssignments AS oa
				INNER JOIN uniprot_accession AS ua ON (ua.accession=oa.accession)
				INNER JOIN uniprot_seq AS seq ON (seq.shortname=ua.shortname)
				INNER JOIN runlist ON (oa.ncbi_taxon = runlist.ncbi_taxon)
				WHERE oa.ortho_id = ?;""", [ortho_grp_id]).fetchall()
			for is_rspecies, shortname, nuc in seqs:
				nuc = zlib.decompress(nuc)
				if is_rspecies:
					rshortnames.append(shortname)
					rseqs.append(nuc)
				else:
					lshortnames.append(shortname)
					lseqs.append(nuc)
			ortho_id_list_file.write( str(ortho_grp_id) +'\t%i\t%i\t'%(len(lseqs),len(rseqs)) + ' '.join(lshortnames) +'\t'+ ' '.join(rshortnames) + '\n')
			grpfile = open(clustalin_dir+'/%i.fasta'%ortho_grp_id, 'w')
			for seq in itertools.chain(lseqs,rseqs):
				grpfile.write(seq.strip())
				grpfile.write('\n')
			grpfile.close()
		ortho_id_list_file.close()
	def output_bsub(self):
		if not os.path.isdir(self.run_name+'/out'):
			os.mkdir(self.run_name+'/out')
		if not os.path.isdir(self.run_name+'/clustalout'):
			os.mkdir(self.run_name+'/clustalout')
		bsubname = self.run_name+'/bsub.sh'
		f = open(bsubname,'w')
		f.write("""#!/bin/bash
pushd `dirname $0` > /dev/null
SCRIPTPATH=`pwd`
RUN_NAME=`basename $SCRIPTPATH`
popd > /dev/null
echo Submitting $RUN_NAME
bsub <<HERE_DOC
#BSUB -u jjoonathan@gmail.com
#BSUB -J $RUN_NAME[1-%i]
#BSUB -e out/err_%%I
#BSUB -o out/out_%%I
#BSUB -cwd $SCRIPTPATH
#BSUB -q short_serial
./bsub1.sh \\${LSB_JOBINDEX}
HERE_DOC"""%(self.num_ortho_grpids))
		f.close()
		bsub1name = self.run_name+'/bsub1.sh'
		g = open(bsub1name,'w')
		g.write("""#!/bin/bash
ORTHOGRP=`awk "NR==$1+1{print "'$1}' ortho_ids`;
clustalw2 -INFILE=clustalin/${ORTHOGRP}.fasta -OUTFILE=clustalout/${ORTHOGRP}.fasta -OUTPUT=FASTA
""")
		g.close()
		subprocess.call(('chmod','+x',bsubname))
		subprocess.call(('chmod','+x',bsub1name))






class SpeciesSubstrAlphaFinder(AlphaFinder):
	def __init__(self, l_species, r_species):
		self.run_name = l_species.split(' ')[0] + '-' + r_species.split(' ')[0]
		self.l_species = l_species
		self.r_species = r_species
	def l_taxon_names(self):
		searchstr = "%%%s%%"%self.l_species
		self.cursor.execute('SELECT name FROM roundup_species WHERE name LIKE ? AND blacklisted=0', [searchstr])
		return [row[0] for row in self.cursor.fetchall()]
	def l_taxon_ids(self):
		searchstr = "%%%s%%"%self.l_species
		self.cursor.execute('SELECT ncbi_taxon FROM roundup_species WHERE name LIKE ? AND blacklisted=0', [searchstr])
		return [row[0] for row in self.cursor.fetchall()]
	def r_taxon_names(self):
		searchstr = "%%%s%%"%self.r_species
		self.cursor.execute('SELECT name FROM roundup_species WHERE name LIKE ? AND blacklisted=0', [searchstr])
		return [row[0] for row in self.cursor.fetchall()]
	def r_taxon_ids(self):
		searchstr = "%%%s%%"%self.r_species
		self.cursor.execute('SELECT ncbi_taxon FROM roundup_species WHERE name LIKE ? AND blacklisted=0', [searchstr])
		return [row[0] for row in self.cursor.fetchall()]
params = cjson.decode(sys.stdin.readline())
print "Running..."
alphafinder = SpeciesSubstrAlphaFinder(params['species1'],params['species2'])
alphafinder.run()

