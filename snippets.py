############################## Load genemaps.json
import itertools
import cjson
gm=cjson.decode(open('genemaps.json').read())
geneFamilies=gm['geneFamilies']
geneToSpecies=gm['geneToSpecies']
geneToOrthologs=gm['geneToOrthologs']

############################## Upload genemaps.json to sqlite
import sqlite3
conn=sqlite3.connect('example.db')
c=conn.cursor()
def entryFromGene(g):
    species = geneToSpecies[g]
    orthlgs = ' '.join(geneToOrthologs[g])
    ontl = None
    return (g, species, orthlgs, ontl)
c.executemany('INSERT INTO Genes VALUES (?,?,?,?)', itertools.imap(entryFromGene,geneToSpecies.iterkeys()))
conn.commit()

############################## Add nucleotide_seq_fasta and uniprot_xml TEXT columns
import sqlite3
conn=sqlite3.connect('example.db')
c=conn.cursor()
c.execute('ALTER TABLE Genes ADD COLUMN nucleotide_seq_fasta TEXT')
c.execute('ALTER TABLE Genes ADD COLUMN uniprot_xml TEXT')
conn.commit()

############################## Add fasta sequences
import itertools, cjson, sqlite3
gs=cjson.decode(open('gene_sequences.json').read())
conn=sqlite3.connect('example.db')
c=conn.cursor()
c.executemany('UPDATE Genes SET nucleotide_seq_fasta==? WHERE uniprot_id==?', itertools.imap(lambda (g,seq): (seq,g), gs.iteritems()))
conn.commit()


############################## Time to get set of all genes
#Result: 4.7s
time python <<HERE_DOC
import sqlite3
conn=sqlite3.connect('example.db')
c=conn.cursor()
c.execute('SELECT uniprot_id FROM Genes')
lst = c.fetchall()
print len(lst)
HERE_DOC

#Result: 
time python <<HERE_DOC
import cjson
gm=cjson.decode(open('genemaps.json').read())
geneFamilies=gm['geneFamilies']
geneToSpecies=gm['geneToSpecies']
geneToOrthologs=gm['geneToOrthologs']
print len(gm.keys())
HERE_DOC
