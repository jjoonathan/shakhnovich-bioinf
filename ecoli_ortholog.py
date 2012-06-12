import sys, code, re, string, itertools, cPickle, urllib, lxml, os, signal
import cookielib
import lxml.html
import lxml.etree
import cjson                        # easy_install python-cjson
import eventlet                     # easy_install eventlet
from eventlet.green import urllib2
from PyQt4.QtCore import *          # yum install PyQt4
from PyQt4.QtGui import *
from PyQt4.QtWebKit import QWebPage
appInstance = QApplication(sys.argv)

class Protein:
    instances={}
    @classmethod
    def withID(self, protID):
        try:
            return Protein.instances[protID]
        except KeyError:
            return Protein(protID)
    def __init__(self, protID):
        self.name = protID
        self.orthologs = {}  # Protein -> distance
        Protein.instances[protID]=self
    def __cmp__(self, other):
        return self.name.__cmp__(other)
    def __str__(self):
        return '<Protein %s>'%self.name
    def __hash__(self):
        return self.name.__hash__()
    def __repr__(self):
        return 'Protein.withId(%s)'%self.name

class Crawler:
    def main(self):
        try:
            os.stat('cache/ecoliNames')
            self.ecoliNames = cjson.decode(open('cache/ecoliNames','r').read())
            print "Found %i cached E. Coli names"%len(self.ecoliNames)
            self.fetch_uncached_orthologs()
        except OSError:
            self.fetch_organism_list()
        if not self.load_gene_list():
            self.construct_gene_list()
            self.save_gene_list()
        self.print_gene_families()

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
        self.ecoliNames = filter(ecoliP, self.organismNames)
        print "Found %i species of E. Coli."%len(self.ecoliNames)
        open('cache/ecoliNames','w').write(cjson.encode(self.ecoliNames))
        self.fetch_uncached_orthologs()

    def fetch_uncached_orthologs(self):
        self.downloader_pool = eventlet.greenpool.GreenPool(size=5)
        pairs_to_download = []
        combs = len(self.ecoliNames)*(len(self.ecoliNames)-1)
        print "That's %i species-species combinations."%combs
        for (l,r) in itertools.combinations(self.ecoliNames,2):
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
            return
        # Fetch the XML file of results
        req = urllib2.Request('http://roundup.hms.harvard.edu'+xmllinks[0])
        response = opener.open(req).read()
        open('cache/%s.xml'%self.cache_name(lSpecies,rSpecies),'w').write(response)
        return self.cache_name(lSpecies,rSpecies)

    def construct_gene_list(self):
        print "Creating global gene list..."
        allpairs = list(itertools.combinations(self.ecoliNames,2))
        lastpercent = -1
        total = len(allpairs)
        for i in xrange(len(allpairs)):
            l, r = allpairs[i]
            self.add_genes_to_list(l,r)
            thispercent = int(((i+1)*100)//total)
            if True:
                print "%i%% (%i/%i)"%(thispercent, i, total)
                lastpercent = thispercent

    def add_genes_to_list(self, lSpecies, rSpecies):
        fname = 'cache/%s.xml'%self.cache_name(lSpecies,rSpecies)
        fcontents = open(fname,'r').read()
        tree = lxml.etree.fromstring(fcontents)
        idnumToProt = {}
        for geneElmt in tree.xpath('//*[local-name()="gene"]'):
            prt = Protein.withID(geneElmt.get('protId'))
            idnumToProt[geneElmt.get('id')] = prt
        for edgeElmt in tree.xpath('//*[local-name()="orthologGroup"]'):
            dist = float(edgeElmt.xpath('./*[local-name()="score"]/@value')[0])
            ids = edgeElmt.xpath('*[local-name()="geneRef"]/@id')
            proteins = map(lambda idno: idnumToProt[idno], ids)
            for p1 in proteins:
                for p2 in proteins:
                    p1.orthologs[p2] = dist

    def save_gene_list(self):
        storedict = {}
        for name, prot in Protein.instances.iteritems():
            orthologs = {}
            for o,dist in prot.orthologs.iteritems():
                orthologs[o.name] = dist
            storedict[name] = orthologs
        open('cache/genelist.json','w').write(cjson.encode(storedict))

    def load_gene_list(self):
        try:
            f = open('cache/genelist.json')
            print "Loading cached gene list..."
            stored_list = f.read()
        except IOError:
            return false
        print "Decoding gene list..."
        raw_genelist = cjson.decode(stored_list)
        print "Initializing gene data structures..."
        del stored_list
        for name in raw_genelist.iterkeys():
            Protein.withID(name)  # initalize all the objects
        for name, edges in raw_genelist.iteritems():
            prt = Protein.withID(name)
            for ortholog, distance in edges.iteritems():
                prt.orthologs[Protein.withID(ortholog)] = distance
        print "Done initializing gene data structures."
                              
            

        

    def print_gene_families(self):
        print "Printing gene families"
        proteins = set(Protein.instances.iterkeys())
        while proteins:
            # Build up the protein family
            visitedMembers = set()
            unvisitedMembers = set((proteins.pop()))
            while unvisitedMembers:
                mbr = unvisitedMembers.pop()
                visitedMembers.add(mbr)
                for child in mbr.orthologs.iterkeys():
                    if child not in visitedMembers:
                        unvisitedMembers.add(child)
            family = visitedMembers
            print 'Family: %s'%family
        

signal.signal(signal.SIGINT, lambda a,b: exit(0))
c = Crawler()
c.main()
appInstance.exec_()
# code.interact('', None, locals())
