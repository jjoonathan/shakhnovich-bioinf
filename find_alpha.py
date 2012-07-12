print "Loading libraries..."
import sys, os, lxml.html, xml.sax, sqlite3, itertools, urllib2, subprocess, random

import eventlet.green.urllib2
import eventlet.green.subprocess as gsub

class AlphaFinder():
	def __init__(self):
		pass
	def run(self):
		self.conn = sqlite3.Connection('data.db', isolation_level=None)
		self.cursor = self.conn.cursor()
		self.update_roundup_species()
		self.download_full_genomes()
	def update_roundup_species(self, force=False):
		sys.stdout.write("Updating roundup species...")
		sys.stdout.flush()
		if self.table_exists('roundup_species') and force != True:
			print " Done. Used cache."
			return
		self.cursor.execute('BEGIN')
		self.cursor.execute('CREATE TABLE roundup_species (ncbi_taxon INTEGER PRIMARY KEY, name TEXT, num_genes INTEGER, proteome_download_complete INTEGER DEFAULT 0);')
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
	def download_full_genomes(self, taxons):
		failed = False
		sys.stdout.write("Downloading full genomes...")
		sys.stdout.flush()
		for ncbi_taxon in taxons:
			if len(self.cursor.execute('SELECT * FROM roundup_species WHERE proteome_download_complete=1 AND ncbi_taxon=?',[ncbi_taxon]).fetchall()) == 1:
				continue
			print "--------------------------- Fetching Complete Proteome of %i -------------"%ncbi_taxon
			url = 'http://www.uniprot.org/uniprot/?query=organism%%3a%i+keyword%%3a181&format=xml'%(ncbi_taxon,)
			args = ('curl', '--create-dirs', '-C', '-', '-o', 'full_proteomes/%i.xml'%ncbi_taxon, url)
			print "\n"+' '.join(args)+"\n"+str(ncbi_taxon)+"\n"
			subproc = subprocess.Popen(args, stdout=sys.stdout, stderr=sys.stderr)
			if subproc.wait() == 0:
				self.cursor.execute('BEGIN')
				self.cursor.execute('UPDATE roundup_species SET proteome_download_complete=1 WHERE ncbi_taxon=?',[ncbi_taxon])
				self.cursor.execute('COMMIT')
			else:
				failed = True
		print " Done."
	def table_exists(self, tname):
		self.cursor.execute('SELECT COUNT(name) FROM sqlite_master WHERE type="table" AND name=?',[tname])
		return bool(self.cursor.fetchone()[0])
class SpeciesSubstrAlphaFinder(AlphaFinder):
	def __init__(self, l_species, r_species):
		self.run_name = l_species.split(' ')[0] + '-' + r_species.split(' ')[0]
		self.l_species = l_species
		self.r_species = r_species
	def is_species_1(species_name):
		return re.search(self.l_species, species_name) != None
	def is_species_2(species_name):
		return re.search(self.r_species, species_name) != None
	def bridges(species1_names, species2_names):
		""" Takes two sets of strings and returns one or more (s1, s2) tuples """
		k12 = filter(lambda s: re.search('K-12',s)!=None, species1_names)[0]
		r = random.Random()
		r.seed(1)
		return [(r.choice(species1_names), r.choice(species2_names)) for i in range(3)]
