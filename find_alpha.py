print "Loading libraries..."
import sys, os, lxml.html, xml.sax, sqlite3, itertools, urllib2, subprocess, random, cookielib
import urllib, os.path, mmap, zlib, re, xml.etree.cElementTree, MySQLdb

import eventlet.green.urllib2
import eventlet.green.subprocess as gsub

class AlphaFinder():
	def __init__(self):
		pass
	def run(self):
		self.conn = MySQLdb.Connection('localhost')
		self.conn.text_factory = str
		self.cursor = self.conn.cursor()
		self.cursor.execute('USE a;')
		if not os.path.isdir(self.run_name):
			os.mkdir(self.run_name)
		self.update_roundup_species_list()
		assert len(self.l_taxon_ids())
		assert len(self.r_taxon_ids())
		self.fetch_roundup_orthology_assignment()
		self.download_full_genomes()
		self.load_genomes_into_db()
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
			self.cursor.execute('INSERT INTO roundup_species (ncbi_taxon, name, num_genes) VALUES (%s,%s,%s)', (ncbi_taxon, name, gene_count))
		self.cursor.execute('COMMIT')
		print " Done."
	def download_full_genomes(self):
		failed_taxons = []
		sys.stdout.write("Downloading full genomes...\n")
		sys.stdout.flush()
		for ncbi_taxon in itertools.chain(self.l_taxon_ids(), self.r_taxon_ids()):
			cachefile_name = '%s/full_proteomes/%i.xml'%(self.run_name,ncbi_taxon)
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
		#self.cursor.execute('SELECT COUNT(name) FROM sqlite_master WHERE type="table" AND name=%s',[tname])
		self.cursor.execute('SELECT COUNT(*) FROM information_schema.tables WHERE table_schema="a" AND table_name=%s',[tname])
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
			return None
		# Fetch the XML file of results
		req = urllib2.Request('http://roundup.hms.harvard.edu'+xmllinks[0])
		response = opener.open(req).read()
		outfile = open('%s/roundup_tmp.xml'%(self.run_name),'w')
		outfile.write(response)
		outfile.close()
		os.rename('%s/roundup_tmp.xml'%(self.run_name), '%s/roundup.xml'%(self.run_name))
		print " Done."
		return True
	def load_genomes_into_db(self):
		msg = "Loading genomes into db...\n"
		sys.stdout.write(msg)
		sys.stdout.flush()
		taxons = self.l_taxon_ids()
		taxons.extend(self.r_taxon_ids())
		if not self.table_exists('uniprot_accession'):
			self.cursor.execute('BEGIN')
			self.cursor.execute("""
				CREATE TABLE uniprot_accession (
					accession VARCHAR(10),
					shortname VARCHAR(255),
					UNIQUE (accession,shortname));""")
			self.cursor.execute("CREATE INDEX uniprot_accession_idx0 ON uniprot_accession (accession);")
			self.cursor.execute("CREATE INDEX uniprot_accession_idx1 ON uniprot_accession (accession,shortname);")
			self.cursor.execute('COMMIT')
		if not self.table_exists('uniprot_dbReference'):
			self.cursor.execute('BEGIN')
			self.cursor.execute("""
				CREATE TABLE uniprot_dbReference (
					shortname VARCHAR(255),
					type VARCHAR(255),
					id VARCHAR(255),
					aa_seq_id TEXT DEFAULT NULL,
					nuc_seq_id TEXT DEFAULT NULL,
					UNIQUE (shortname, type, id)
				);""")
			self.cursor.execute("CREATE INDEX uniprot_dbReference_idx0 ON uniprot_dbReference (shortname);")
			self.cursor.execute("CREATE INDEX uniprot_dbReference_idx1 ON uniprot_dbReference (type);")
			self.cursor.execute("CREATE INDEX uniprot_dbReference_idx2 ON uniprot_dbReference (shortname,type,id);")
			self.cursor.execute('COMMIT')
		if not self.table_exists('uniprot'):
			self.cursor.execute('BEGIN')
			self.cursor.execute("""
				CREATE TABLE uniprot (
					shortname VARCHAR(255) UNIQUE,
					zippedxml BLOB,
					aa TEXT DEFAULT NULL,
					nuc TEXT DEFAULT NULL,
					aa_len INTEGER DEFAULT 0,
					mass REAL DEFAULT 0
				);""")
			self.cursor.execute("CREATE INDEX uniprot_idx0 ON uniprot (shortname);")
			self.cursor.execute('COMMIT')
		if not self.table_exists('uniprot_seq'):
			self.cursor.execute('BEGIN')
			self.cursor.execute("""
				CREATE TABLE uniprot_seq (
					shortname VARCHAR(255) UNIQUE,
					aa TEXT,
					aa_valid INTEGER DEFAULT 0,
					nuc TEXT DEFAULT NULL,
					nuc_valid INTEGER DEFAULT 0,
					aa_len INTEGER,
					mass INTEGER
				);""")
			self.cursor.execute("CREATE INDEX uniprot_seq_idx0 ON uniprot_seq (shortname);")
			self.cursor.execute("CREATE INDEX uniprot_seq_idx1 ON uniprot_seq (aa_valid);")
			self.cursor.execute("CREATE INDEX uniprot_seq_idx2 ON uniprot_seq (nuc_valid);")
			self.cursor.execute('COMMIT')
		for ncbi_taxon in taxons:
			self.cursor.execute('SELECT proteome_in_db FROM roundup_species WHERE ncbi_taxon=%s',[ncbi_taxon])
			if self.cursor.fetchone()[0] == 1:
				continue
			file_name = '%s/full_proteomes/%i.xml'%(self.run_name,ncbi_taxon)
			with open(file_name,'r+b') as f:
				print "\tProcessing: %s"%file_name
				fmap = mmap.mmap(f.fileno(), 0)
				for matchobj in re.finditer(r"(<entry.{100,}?</entry>)", fmap, flags=re.MULTILINE+re.DOTALL):
					self.cursor.execute('BEGIN')
					entry_str = matchobj.group(0)
					strz = zlib.compress(entry_str)
					tree = xml.etree.cElementTree.fromstring(entry_str)
					shortname = tree.findall('name')[0].text
					self.cursor.execute("REPLACE INTO uniprot (shortname, zippedxml) VALUES (%s,%s)",(shortname,strz))
					for accession in tree.findall('accession'):
						self.cursor.execute("""
							REPLACE INTO uniprot_accession (accession, shortname)
							VALUES (%s,%s)""",\
							(accession.text,shortname))
					for dbref in tree.findall('dbReference'):
						db_type, db_id, aa_seq_id, nuc_seq_id = dbref.get('type'), dbref.get('id'), None, None
						for prop in dbref.findall('property'):
							if prop.get('type') == 'protein sequence ID':
								aa_seq_id = prop.get('value')
							elif prop.get('type') == 'nucleotide sequence ID':
								nuc_seq_id = prop.get('value')
						self.cursor.execute("""
							REPLACE INTO uniprot_dbReference (shortname, type, id, aa_seq_id, nuc_seq_id)
							VALUES (%s,%s,%s,%s,%s);""",\
							(shortname, db_type, db_id, aa_seq_id, nuc_seq_id))
					for seq in tree.findall('sequence'):
						aa = seq.text
						aa_len = int(seq.get('length'))
						mass = float(seq.get('mass'))
						self.cursor.execute("UPDATE uniprot SET aa=%s, aa_len=%s, mass=%s WHERE shortname=%s",[aa,aa_len,mass,shortname])
					self.cursor.execute('COMMIT')
			#self.cursor.execute('BEGIN')
			self.cursor.execute('UPDATE roundup_species SET proteome_in_db=1 WHERE ncbi_taxon=%s',[ncbi_taxon])
			self.cursor.execute('COMMIT')

class SpeciesSubstrAlphaFinder(AlphaFinder):
	def __init__(self, l_species, r_species):
		self.run_name = l_species.split(' ')[0] + '-' + r_species.split(' ')[0]
		self.l_species = l_species
		self.r_species = r_species
	def l_taxon_names(self):
		searchstr = "%%%s%%"%self.l_species
		self.cursor.execute('SELECT name FROM roundup_species WHERE name LIKE %s AND blacklisted=0', [searchstr])
		return [row[0] for row in self.cursor.fetchall()]
	def l_taxon_ids(self):
		searchstr = "%%%s%%"%self.l_species
		self.cursor.execute('SELECT ncbi_taxon FROM roundup_species WHERE name LIKE %s AND blacklisted=0', [searchstr])
		return [row[0] for row in self.cursor.fetchall()]
	def r_taxon_names(self):
		searchstr = "%%%s%%"%self.r_species
		self.cursor.execute('SELECT name FROM roundup_species WHERE name LIKE %s AND blacklisted=0', [searchstr])
		return [row[0] for row in self.cursor.fetchall()]
	def r_taxon_ids(self):
		searchstr = "%%%s%%"%self.r_species
		self.cursor.execute('SELECT ncbi_taxon FROM roundup_species WHERE name LIKE %s AND blacklisted=0', [searchstr])
		return [row[0] for row in self.cursor.fetchall()]
alphafinder = SpeciesSubstrAlphaFinder('Escherichia coli', 'Salmonella enterica subsp. enterica serovar Typhimurium')
alphafinder.run()

