print "Loading libraries..."
import sys, os, lxml.html, xml.sax, sqlite3, itertools, urllib2

import eventlet.green.urllib2
import eventlet.green.subprocess as gsub

class AlphaFinder():
	def __init__(self):
		self.conn = sqlite3.Connection('data.db', isolation_level=None)
		self.cursor = self.conn.cursor()
		self.update_roundup_species()
	def update_roundup_species(self, force=False):
		sys.stdout.write("Updating roundup species...")
		sys.stdout.flush()
		if self.table_exists('roundup_species') and force != True:
			print " Done. Used cache."
			return
		self.cursor.execute('BEGIN')
		self.cursor.execute('CREATE TABLE roundup_species (ncbi_taxon INTEGER PRIMARY KEY, name TEXT, num_genes INTEGER);')
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
	def table_exists(self, tname):
		self.cursor.execute('SELECT COUNT(name) FROM sqlite_master WHERE type="table" AND name=?',[tname])
		return bool(self.cursor.fetchone()[0])
print "Running main..."		
AlphaFinder() 
