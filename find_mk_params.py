import sys, code, re, string, itertools, cPickle, urllib, lxml
import cookielib, time, os, signal, random, traceback, subprocess
import select, httplib, sqlite3

import lxml.html                    # easy_install lxml
import lxml.etree

import cjson                        # easy_install python-cjson

import eventlet                     # easy_install eventlet
from eventlet.green import urllib2
import eventlet.green.subprocess as gsub

from PyQt4.QtCore import *          # yum install PyQt4
from PyQt4.QtGui import *
from PyQt4.QtWebKit import QWebPage

appInstance = QApplication(sys.argv)
run_name = 'EColi-Salmonella'
email = 'jjoonathan@gmail.com'
def is_species_1(species_name):
    return re.search('Escherichia', species_name) != None
def is_species_2(species_name):
    return re.search('Salmonella', species_name) != None
def bridges(species1_names, species2_names):
    """ Takes two sets of strings and returns one or more (s1, s2) tuples """
    k12 = filter(lambda s: re.search('K-12',s)!=None, species1_names)[0]
    return [(k12, species2_names[0]), (k12, species2_names[1]), (k12, species2_names[2])]


# Folders
#   ecoli_ortholog.py
#   <run_name>/species_names.json    // keys: allSpecies, species1Names, species2Names
#   <run_name>/roundup/subspeciesA---subspeciesB.xml    //raw cache of roundup response
#   <run_name>/genemaps.json         // keys: geneToSpecies, geneToOrthologs, geneFamilies 
#   <run_name>/gene_sequences.json   // flat: gene_name -> gene seq in FASTA format+'\n'
#   <run_name>/gene_families.json   // flat: gene_name -> gene seq in FASTA format+'\n'
#   <run_name>/clustalin/92347894    // unaligned gene family input
#   <run_name>/clustalout/09833900   // aligned gene family output

class Crawler:
    geneToOrthologs = {}
    geneToSpecies = {}
    geneSequences = {}
    geneFamilies = None  # A list of sets containing the proteins in that family
    species1Names = None
    species2Names = None
    speciesPairs = []
    malformedXMLFiles = []
    def main(self):
        self.ensure_output_dirs_exist()
        self.connect_to_output_db()
        
        self.find_species_supported_by_roundup()
        self.determine_target_species()
        self.build_list_of_species_pairs_to_scan_for_orthologs()

        self.fetch_uncached_orthologs()
        self.load_gene_list()
        self.find_gene_families()
        # self.output_gene_families()
        self.fetch_gene_sequences()
        self.align_families()
        self.mktest_families()
        exit(0)

    def ensure_output_dirs_exist(self):
        if not os.path.isdir(run_name):
            os.mkdir(run_name)
	if not os.path.isdir(run_name+'/clustalin'):
            os.mkdir(run_name+'/clustalin')
	if not os.path.isdir(run_name+'/clustalout'):
            os.mkdir(run_name+'/clustalout')
	if not os.path.isdir(run_name+'/roundup'):
            os.mkdir(run_name+'/roundup')
	if not os.path.isdir(run_name+'/mktest_out'):
            os.mkdir(run_name+'/mktest_out')

    def connect_to_output_db(self):
        dbfile = run_name+'/genes.db'
        self.db_connection = sqlite3.connect(dbfile)
        self.db = self.db_connection.cursor()

############################################# find species-species pairs to lookup in roundup #############################################
    def find_species_supported_by_roundup(self):
        try:
            print "Loading list of species supported by roundup from db..."
            self.db.execute('SELECT name FROM Species')
            self.allSpecies = map(lambda tuple: tuple[0], self.db.fetchall())
            return
        except sqlite3.OperationalError:
            pass
        print "\x1B[1FFetching list of species supported by roundup..."
        self.db.execute('CREATE TABLE Species (\n'+\
                            'id INTEGER PRIMARY KEY,\n'+\
                            'name TEXT)')
        while True:
            try:
                result = urllib2.urlopen('http://roundup.hms.harvard.edu/retrieve/').read()
                self.allSpecies = re.findall('{\s*"name":\s*"([^"]+)"\s*,',result)
                allSpeciesAsSQLTuples = map(lambda i: (i,self.allSpecies[i],), xrange(len(self.allSpecies)))
                self.db.executemany('INSERT INTO Species (id,name) VALUES (?)', allSpeciesAsSQLTuples)
                self.db_connection.commit()
                break
            except Exception as e:
                print "Error: could not load gene list because: %s"%e
                time.sleep(5)
    
    def determine_target_species(self):
        self.species1Names = filter(is_species_1, self.allSpecies)
        self.species2Names = filter(is_species_2, self.allSpecies)
        self.bridgePairs = bridges(self.species1Names, self.species2Names)

############################################# fetch_uncached_orthologs #############################################
    def build_list_of_species_pairs_to_scan_for_orthologs(self):
        combs1 = len(self.species1Names)*(len(self.species1Names)-1)/2
        combs2 = len(self.species2Names)*(len(self.species2Names)-1)/2
        self.speciesPairs.extend(self.bridgePairs)
        self.speciesPairs.extend(itertools.combinations(self.species1Names,2))
        self.speciesPairs.extend(itertools.combinations(self.species2Names,2))
        print "That's %i combinations of species1, %i of species2, %i bridges."%(combs1,combs2,len(self.bridgePairs))
        
    def fetch_uncached_orthologs(self):
        self.downloader_pool = eventlet.greenpool.GreenPool(size=5)
        self.pairs_to_download = []
        numPairs = len(self.speciesPairs)
        for i in xrange(numPairs):
            l,r = self.speciesPairs[i]
            if i%20 == 0:
                print "%i%% (%i/%i)\x1B[1F"%(int(i*100.0/numPairs),i,numPairs)
            if not os.path.isfile('%s/roundup/%s.xml'%(run_name,self.cache_name(l,r))):
                self.pairs_to_download.append((l,r))
        num_to_dl = len(self.pairs_to_download)
        print "Fetching %i uncached combinations of species..."%num_to_dl
        pdp = self.downloader_pool.imap(self.fetch_pair, self.pairs_to_download)
        i=0
        for response in pdp:
            i+=1
            cachename = self.cache_name(*response)
            print "%i%% (%i/%i): %s\x1B[1F"%(int(i*100.0/num_to_dl), i, num_to_dl, cachename)

    def cache_name(self, lSpecies, rSpecies):
        name = lSpecies+'---'+rSpecies
        valid_chrs = '-_.() %s%s'%(string.ascii_letters, string.digits)
        filename = ''.join(c for c in name if c in valid_chrs)
        return filename

    def fetch_pair(self, (lSpecies, rSpecies)):
        while True:
            try:
                self.attempt_fetch_pair((lSpecies,rSpecies))
                break
            except urllib2.URLError as e:
                print "Error fetching (%s,%s): %s"%(lSpecies,rSpecies,e)
                continue
        return (lSpecies,rSpecies)

    def attempt_fetch_pair(self, (lSpecies, rSpecies)):
        # First grab the CSRF, cookies with a GET request
        cjar = cookielib.CookieJar()
        opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cjar))
        request = urllib2.Request('http://roundup.hms.harvard.edu/retrieve/')
        response = opener.open(request).read()
        rtree = lxml.html.document_fromstring(response)
        key = rtree.xpath('//input[@name="csrfmiddlewaretoken"]/@value')[0]
        # Now use the session ID to submit the form
        params = []
        params.append(('csrfmiddlewaretoken',key))
        params.append(('genomes_filter','E'))
        params.append(('genomes_filter','B'))
        params.append(('genomes_filter','A'))
        params.append(('genomes_filter','V'))
        params.append(('genome_choices','Escherichia coli 536'))
        params.append(('genomes', '%s\r\n%s\r\n'%(lSpecies,rSpecies)))
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
        rtree = lxml.html.document_fromstring(response_str)
        links = rtree.xpath('//a/@href')
        xmllinks = filter(lambda s: s.find('tt=xml')!=-1, links)
        if len(xmllinks)!=1:
            print "ERROR: bad number of XML links (%s)"%xmllinks
            return None
        # Fetch the XML file of results
        req = urllib2.Request('http://roundup.hms.harvard.edu'+xmllinks[0])
        response = opener.open(req).read()
        open('%s/roundup/%s.xml'%(run_name,self.cache_name(lSpecies,rSpecies)),'w').write(response)
        return self.cache_name(lSpecies,rSpecies)

############################################# load_gene_list #############################################
    def load_gene_list(self):
        if os.path.isfile('%s/genemaps.json'%run_name):
            print "Loading cached gene list..."
            f = open('%s/genemaps.json'%run_name)
            stored_list = f.read()
            print "Decoding gene list..."
            save_dict = cjson.decode(stored_list)
            self.geneToOrthologs = save_dict['geneToOrthologs']
            self.geneToSpecies = save_dict['geneToSpecies']
            self.geneFamilies = save_dict.get('geneFamilies',None)
            if self.geneFamilies != None:
                for i in xrange(len(self.geneFamilies)):
                    self.geneFamilies[i] = set(self.geneFamilies[i])
        else:
            self.construct_gene_list()
            self.save_gene_list()

    def construct_gene_list(self):
        print "Creating global gene list..."
        lastpercent = -1
        total = len(self.speciesPairs)
        for i in xrange(len(self.speciesPairs)):
            l, r = self.speciesPairs[i]
            self.add_genes_to_list(l,r)
            if True:
                print "%i%% (%i/%i)\x1B[1F"%(int(((i+1)*100)//total), i+1, total)
        if self.malformedXMLFiles:
            for f in self.malformedXMLFiles:
                os.unlink(f)
            print "Malformed XML files purged. Please run script again."
            exit(0)
                

    def add_genes_to_list(self, lSpecies, rSpecies):
        fname = '%s/roundup/%s.xml'%(run_name,self.cache_name(lSpecies,rSpecies))
        fcontents = open(fname,'r').read()
        try:
            tree = lxml.etree.fromstring(fcontents)
        except lxml.etree.XMLSyntaxError:
            print "Problem with file '%s'."%fname
            self.malformedXMLFiles.append(fname)
            return
        idnumToProt = {}
        for speciesElmt in tree.xpath('//*[local-name()="species"]'):
            species_name = speciesElmt.get('name')
            for geneElmt in speciesElmt.xpath('.//*[local-name()="gene"]'):
                prtId = geneElmt.get('protId')
                if prtId not in self.geneToOrthologs:
                    self.geneToOrthologs[prtId] = {}
                self.geneToSpecies[prtId] = species_name
                idnumToProt[geneElmt.get('id')] = prtId
        for edgeElmt in tree.xpath('//*[local-name()="orthologGroup"]'):
            dist = float(edgeElmt.xpath('./*[local-name()="score"]/@value')[0])
            ids = edgeElmt.xpath('*[local-name()="geneRef"]/@id')
            proteins = map(lambda idno: idnumToProt[idno], ids)
            for p1 in proteins:
                for p2 in proteins:
                    if p1 == p2:
                        continue
                    self.geneToOrthologs[p1][p2] = dist
                    self.geneToOrthologs[p2][p1] = dist

    def save_gene_list(self):
        save_dict = {'geneToOrthologs' : self.geneToOrthologs}
        save_dict['geneToSpecies'] = self.geneToSpecies
        if self.geneFamilies != None:
            save_dict['geneFamilies'] = map(list,self.geneFamilies)
        open('%s/genemaps.json'%run_name,'w').write(cjson.encode(save_dict))


############################################# find_gene_families #############################################
    def find_gene_families(self):
        if self.geneFamilies != None:
            return
        print "Determining gene families..."
        self.geneFamilies = []
        proteins = set(self.geneToOrthologs.iterkeys())
        while proteins:
            # Build up the protein family
            visitedMembers = set()
            unvisitedMembers = set((proteins.pop(),))
            while unvisitedMembers:
                mbr = unvisitedMembers.pop()
                visitedMembers.add(mbr)
                for child in self.geneToOrthologs[mbr].iterkeys():
                    if child not in visitedMembers:
                        unvisitedMembers.add(child)
                        proteins.discard(child)
            family = visitedMembers
            self.geneFamilies.append(family)
        self.save_gene_list()

############################################# output_gene_families #############################################
    def output_gene_families(self):
        if os.path.isfile('%s/gene_families.txt'%run_name):
            return
        print "Writing output..."
        species_name_to_col_no = {}
        fout = open('%s/gene_families.txt'%run_name,'w')
        fout.write('Family#\t')
        # Print header
        for i in xrange(len(self.allSpecies)):
            species_name = self.allSpecies[i]
            species_name_to_col_no[species_name] = i
            fout.write(species_name+'\t')
        # Loop through each family
        num_families = len(self.geneFamilies)
        num_species = len(self.allSpecies)
        for familyNum in xrange(num_families):
            sys.stdout.write("%i%% (%i/%i)\x1B[1G"%(int(familyNum*100.0/num_families), familyNum, num_families))
            sys.stdout.flush()
            fout.write('\n%s\t'%str(familyNum))
            family = self.geneFamilies[familyNum]
            genes = []
            # Sort genes into (index, species) cells
            for g in self.geneFamilies:
                genes.append(set())
            for g in family:
                species_name = self.geneToSpecies[g]
                genes[species_name_to_col_no[species_name]].add(g)
            genes = map(list, genes)
            max_genes_in_species = 0
            for gset in genes:
                max_genes_in_species = max(max_genes_in_species, len(gset))
            # Print out the cells
            for j in xrange(max_genes_in_species):
                if j!=0:
                    fout.write('\t')
                for i in xrange(num_species):
                    try:
                        fout.write(genes[i][j]+'\t')
                    except IndexError:
                        fout.write('\t')
                fout.write('\n')
            fout.write('\n')

############################################# fetch_gene_sequences #############################################
    def fetch_gene_sequences(self):
        print "Fetching family FASTA files..."
        try:
            self.geneSequences = cjson.decode(open('%s/gene_sequences.json'%run_name).read())
        except:
            self.geneSequences = {}
        genes_to_fetch = list(set(self.geneToSpecies.iterkeys())-set(self.geneSequences.iterkeys()))
        if len(genes_to_fetch) == 0:
            print "All genes already fetched."
            return
        self.genes_fetched=0
        self.total_genes_to_fetch = len(genes_to_fetch)
        while len(genes_to_fetch) > 0:
            # Chop off a chunk of 1000 genes, fetch them, write them to the output file
            current_chunk = genes_to_fetch[:1000]
            genes_to_fetch = genes_to_fetch[1000:]
            for g in downloader_pool.imap(fetch_gene, current_chunk):
                i += 1
                total_genes_to_fetch = len(genes_to_fetch) + len(current_chunk)                
            fname='%s/gene_sequences.json'%run_name
            if os.path.isfile(fname):
                os.rename(fname,fname+'_old')
                open(fname,'w').write(cjson.encode(self.geneSequences))
                os.unlink(fname+'_old')
            else:
                open(fname,'w').write(cjson.encode(self.geneSequences))
        print "Done."


    def fetch_genes(genes):
        downloader_pool = eventlet.greenpool.GreenPool(size=8)
        while len(genes)>0:
            pool = downloader_pool.imap(fetch_gene, list(genes))
            genes = []
            for success, gene in pool:
                if not success:
                    genes.append(gene)
                else:
                    self.genes_fetched += 1
                i, tot = self.genes_fetched, self.total_genes_to_fetch
                print "%i%% (%i/%i)\x1B[1F"%(i*100.0/tot, i, tot)

    def fetch_gene(g):
        fetchurl = 'http://www.uniprot.org/uniprot/%s.xml'%g
        try:
            response = urllib2.urlopen(fetchurl).read()
            xmldoc = lxml.etree.fromstring(response)
            xmlns = {'up':'http://uniprot.org/uniprot'}
            seq_refs = xmldoc.xpath('//up:dbReference[@type="EMBL"]/up:property[@type="protein sequence ID"]/@value',namespaces=xmlns)
            if len(seq_refs)==0:
                print "Error: %s had no associated sequence references."%g
                return
        except Exception as e:
            print "Error '%s'->%s"%(fetchurl,e)
            return (False, g)
        try:
            url = "http://www.ebi.ac.uk/ena/data/view/%s&display=fasta"%(seq_refs[0],)
            self.geneSequences[g] = urllib2.urlopen(url).read()
        except Exception as e:
            print "Error2 %s %s"%(url, e)
            return (False, g)
        return (True, g)

############################################# clustal #############################################               
    def align_families(self):
        print "Aligning families..."
        num_families = len(self.geneFamilies)
        self.n_values = {}  # gene family id -> number of genes in family of species1 origin
        stderr = open('%s/clustal_stderr.txt'%run_name,'w')
        if not os.path.isdir(run_name+'/clustalin'):
            os.mkdir(run_name+'/clustalin')
        if not os.path.isdir(run_name+'/clustalout'):
            os.mkdir(run_name+'/clustalout')
        aligner_pool = eventlet.greenpool.GreenPool(size=16)
        i = 0
        for result in aligner_pool.imap(self.align_family, xrange(num_families)):
            i += 1
            sys.stdout.write("%i%% (%i/%i)\x1B[1G"%(int(i*100.0/num_families), i, num_families))
            sys.stdout.flush()            
        print "All family alignments processed."
        self.n_values = cjson.decode(open(run_name+'/n_values.json').read())
    
    def align_family(self, familynum):
        i = familynum
        family = list(self.geneFamilies[i])
        family.sort()
        familyid = repr(family).__hash__()
        ifname = '%s/clustalin/%s.fasta'%(run_name,familyid)
        ofname = '%s/clustalout/%s.fasta'%(run_name,familyid)
        if os.path.isfile(ofname):
            return
        s1,s2 = [],[]
        for g in family:  # Sort the genes into species1, species2 buckets
            s = self.geneToSpecies[g]
            if s in self.species1Names:
                s1.append(g)
            else:
                s2.append(g)
        if len(s1)<10 or len(s2)<10:
            return
        # print '+%i:%i family'%(len(s1),len(s2))
        clustal_inpt_f = open(ifname,'w')
        for g in itertools.chain(s1,s2):
            clustal_inpt_f.write(self.geneSequences.get(g,''))
        clustal_inpt_f.close()
        args = ('clustalw2','-INFILE='+ifname,'-OUTFILE='+ofname,'-OUTPUT=FASTA')
        try:
            pass
            #PIPE = gsub.PIPE
            #proc = gsub.Popen(args, stdin=None, stdout=None, stderr=stderr)
            #proc.wait()
        except OSError as e:
            print "clustalw2 terminated: %s"%e
        self.n_values[str(familyid)] = len(s1)
        if i%10 == 0 or i == len(self.geneFamilies):
            fname = run_name+'/n_values.json'
            if os.path.isfile(fname):
                os.rename(fname,fname+'_old')
                open(run_name+'/n_values.json','w').write(cjson.encode(self.n_values))
                os.unlink(fname+'_old')
            else:
                open(run_name+'/n_values.json','w').write(cjson.encode(self.n_values))

    def mktest_families(self):
        print "MKtesting families..."
        num_families = len(self.geneFamilies)
        mktest_pool = eventlet.greenpool.GreenPool(size=4)
        i = 0
        for result in mktest_pool.imap(self.mktest_family, xrange(num_families)):
            i += 1
            sys.stdout.write("%i%% (%i/%i)\x1B[1G"%(int(i*100.0/num_families), i, num_families))
            sys.stdout.flush()

    def mktest_family(self, i):
        family = list(self.geneFamilies[i])
        family.sort()
        familyid=repr(family).__hash__()
        ifname = '%s/clustalout/%s.fasta'%(run_name,familyid)
        ofname = '%s/mktest_out/%s.txt'%(run_name,familyid)
        if os.path.isfile(ofname):
            return
        try:
            in_num = self.n_values[str(familyid)]
        except KeyError as e:
            if os.path.isfile(ifname):
                print "Failed to mktest family %s (no n)."%str(familyid)
            return
        args = ('mktest','-i',ifname,'-n',str(in_num))
        ofile = open(ofname,'w')
        try:
            proc = gsub.Popen(args,stdout=ofile,stderr=ofile)
            proc.wait()
        except OSError:
            pass
        ofile.close()
        
            
def ctrlc(a,b):
    traceback.print_stack()
    exit(0)
if __name__ == '__main__':
	signal.signal(signal.SIGINT, ctrlc)
	c = Crawler()
	c.main()
