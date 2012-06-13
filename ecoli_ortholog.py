import sys, code, re, string, itertools, cPickle, urllib, lxml
import cookielib, time, os, signal
import lxml.html
import lxml.etree
import cjson                        # easy_install python-cjson
import eventlet                     # easy_install eventlet
from eventlet.green import urllib2
from PyQt4.QtCore import *          # yum install PyQt4
from PyQt4.QtGui import *
from PyQt4.QtWebKit import QWebPage
appInstance = QApplication(sys.argv)

class Crawler:
    geneToOrthologs = {}
    geneToSpecies = {}
    geneFamilies = None  # A list of sets containing the proteins in that family
    targetSpecies = None
    def main(self):
        try:
            os.stat('cache/ecoliNames')
            self.ecoliNames = cjson.decode(open('cache/ecoliNames','r').read())
            print "Found %i cached E. Coli names"%len(self.ecoliNames)
            self.fetch_uncached_orthologs()
        except OSError:
            self.fetch_organism_list()
            while self.targetSpecies == None:
                time.sleep(.05)
                appInstance.processEvents()
        if not self.load_gene_list():
            self.construct_gene_list()
            self.save_gene_list()
        if self.geneFamilies == None:
            self.find_gene_families()
            self.save_gene_list()
        self.output_gene_families()
        exit(0)

    def fetch_organism_list(self):
        self.webpage = QWebPage()
        self.webpage.loadFinished.connect(self.process_organism_list)
        self.webpage.loadStarted.connect(self.fetch_organism_list_started)
        self.webpage.mainFrame().load(QUrl('http://roundup.hms.harvard.edu/retrieve/'))

    def fetch_organism_list_started(self):
        print "Fetching organism list..."

    def process_organism_list(self, bool):
        print "Processing organism list..."
        organisms_query = 'select#id_genome_choices'
        organisms_element = self.webpage.mainFrame().findAllElements(organisms_query).at(0)
        elmt = organisms_element.firstChild()
        self.organismNames = []
        while True:
            if elmt == organisms_element.lastChild():
                break
            self.organismNames.append(str(elmt.attribute('value')))
            elmt = elmt.nextSibling()
        print "Found %i organisms."%len(self.organismNames)
        ecoliP = lambda s: re.search('Escherichia',s) != None
        self.targetSpecies = filter(ecoliP, self.organismNames)
        print "Found %i species of E. Coli."%len(self.targetSpecies)
        open('cache/targetSpecies','w').write(cjson.encode(self.targetSpecies))
        self.fetch_uncached_orthologs()

    def fetch_uncached_orthologs(self):
        self.downloader_pool = eventlet.greenpool.GreenPool(size=5)
        pairs_to_download = []
        combs = len(self.targetSpecies)*(len(self.targetSpecies)-1)/2
        print "That's %i species-species combinations."%combs
        for (l,r) in itertools.combinations(self.targetSpecies,2):
            try:
                os.stat('cache/%s.xml'%self.cache_name(l,r))
            except OSError:
                pairs_to_download.append((l,r))
        print "Fetching %i uncached combinations of species..."%len(pairs_to_download)
        pdp = self.downloader_pool.imap(self.fetch_pair, pairs_to_download)
        for response in pdp:
            print "Received: %s"%response
        
    def cache_name(self, lSpecies, rSpecies):
        name = lSpecies+'---'+rSpecies
        valid_chrs = '-_.() %s%s'%(string.ascii_letters, string.digits)
        filename = ''.join(c for c in name if c in valid_chrs)
        return filename

    def fetch_pair(self, (lSpecies, rSpecies)):
        print 'Requesting: %s'%self.cache_name(lSpecies,rSpecies)
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
        params.append(('evalue','1e-20'))
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
            return
        # Fetch the XML file of results
        req = urllib2.Request('http://roundup.hms.harvard.edu'+xmllinks[0])
        response = opener.open(req).read()
        open('cache/%s.xml'%self.cache_name(lSpecies,rSpecies),'w').write(response)
        return self.cache_name(lSpecies,rSpecies)

    def construct_gene_list(self):
        print "Creating global gene list..."
        allpairs = list(itertools.combinations(self.targetSpecies,2))
        lastpercent = -1
        total = len(allpairs)
        for i in xrange(len(allpairs)):
            l, r = allpairs[i]
            self.add_genes_to_list(l,r)
            thispercent = int(((i+1)*100)//total)
            if True:
                print "%i%% (%i/%i)"%(thispercent, i+1, total)
                lastpercent = thispercent

    def add_genes_to_list(self, lSpecies, rSpecies):
        fname = 'cache/%s.xml'%self.cache_name(lSpecies,rSpecies)
        fcontents = open(fname,'r').read()
        tree = lxml.etree.fromstring(fcontents)
        idnumToProt = {}
        for speciesElmt in tree.xpath('//*[local-name()="species"]'):
            species_name = speciesElmt.get('name')
            for geneElmt in speciesElmt.xpath('//*[local-name()="gene"]'):
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
        open('genemaps.json','w').write(cjson.encode(save_dict))

    def load_gene_list(self):
        try:
            f = open('genemaps.json')
            print "Loading cached gene list..."
            stored_list = f.read()
        except IOError:
            return False
        print "Decoding gene list..."
        save_dict = cjson.decode(stored_list)
        self.geneToOrthologs = save_dict['geneToOrthologs']
        self.geneToSpecies = save_dict['geneToSpecies']
        self.geneFamilies = save_dict.get('geneFamilies',None)
        if self.geneFamilies != None:
            for i in xrange(len(self.geneFamilies)):
                self.geneFamilies[i] = set(self.geneFamilies[i])
        return True

    def find_gene_families(self):
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
        # print "Checking consistency of gene families..."
        # for i in xrange(len(self.geneFamilies)):
        #     for j in xrange(len(self.geneFamilies)):
        #         gf1, gf2 = self.geneFamilies[i], self.geneFamilies[j]
        #         if gf1==gf2:
        #             continue
        #         for g in gf1:
        #             if g in gf2:
        #                 print "%s in %s and %s"%(g, gf1, gf2)
        #                 code.interact(None,None,locals())
        # print "Gene families determined to be consistent."
        

    def output_gene_families(self):
        print "Writing output..."
        species_name_to_col_no = {}
        fout = open('gene_families.txt','w')
        fout.write('Family#\t')
        # Print header
        for i in xrange(len(self.targetSpecies)):
            species_name = self.targetSpecies[i]
            species_name_to_col_no[species_name] = i
            fout.write(species_name+'\t')
        # Loop through each family
        lastPercentDone = -1
        num_families = len(self.geneFamilies)
        num_species = len(self.targetSpecies)
        for familyNum in xrange(num_families):
            currentPercentDone = int(familyNum*100/num_families)
            if currentPercentDone != lastPercentDone:
                lastPercentDone = currentPercentDone
                print "%i%% (%i/%i)"%(currentPercentDone, familyNum, num_families)
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
                
        

signal.signal(signal.SIGINT, lambda a,b: exit(0))
c = Crawler()
c.main()
