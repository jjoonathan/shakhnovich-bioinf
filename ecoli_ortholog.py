import sys, code, re, string, itertools, cPickle, urllib, lxml, os
from lxml import etree
import eventlet                     # easy_install eventlet
from eventlet.green import urllib2
from PyQt4.QtCore import *          # yum install PyQt4
from PyQt4.QtGui import *
from PyQt4.QtWebKit import QWebPage
appInstance = QApplication(sys.argv)

class crawler:
    def main(self):
        self.fetch_organism_list()

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
        self.fetch_pairwise_orthologs()

    def fetch_pairwise_orthologs(self):
        self.downloader_pool = eventlet.greenpool.GreenPool(size=5)
        pairs_to_download = []
        for (l,r) in itertools.combinations(self.ecoliNames,2):
            try:
                os.stat('cache/%s'%self.cache_name(l,r))
            except OSError:
                pairs_to_download.append((l,r))
        print "That's %i uncached combinations of species"%len(pairs_to_download)
        pdp = self.downloader_pool.imap(self.fetch_pair, pairs_to_download)
        for response in pdp:
            print response
        
    def cache_name(self, lSpecies, rSpecies):
        name = lSpecies+'---'+rSpecies
        valid_chrs = '-_.() %s%s'%(string.ascii_letters, string.digits)
        filename = ''.join(c for c in name if c in valid_chrs)

    def fetch_pair(self, (lSpecies, rSpecies)):
        print '-------------- %s ------------'%self.cache_name(lSpecies,rSpecies)
        params = []
        params.append(('genomes_filter','E'))
        params.append(('genomes_filter','B'))
        params.append(('genomes_filter','A'))
        params.append(('genomes_filter','V'))
        params.append(('genomes', '%s\n%s'%(lSpecies,rSpecies)))
        params.append(('divergence','.8'))
        params.append(('evalue','1e-5'))
        params.append(('distance_lower_limit',''))
        params.append(('distance_upper_limit',''))
        loc = 'http://roundup.hms.harvard.edu/retrieve/'
        dat = urllib.urlencode(params)
        # First grab the session ID with a GET request
        # We have to use Qt again because etree barfs on 0x2c
        # UPDATE: no run loop integration for Qt, green trees :(
        # wp = QWebPage()
        # callback = lambda success: self.fetch_pair_2(wp, success, lSpecies, rSpecies)
        # wp.loadFinished.connect(callback)
        # wp.mainFrame().load(QUrl('http://roundup.hms.harvard.edu/retrieve/'))
        response = urllib2.urlopen('http://roundup.hms.harvard.edu/retrieve/').read()
        parser = etree.XMLParser(recover=True)
        rtree = etree.fromstring(response, parser)
        key = rtree.find('input[@id=csrfmiddlewaretoken]').attrib['value']
        print "(%s,%s): %s"%(lSpecies,rSpecies,key)
        
        


    def fetch_pair_2(self, success, webpage, lSpecies, rSpecies):
        key_query = 'input#csrfmiddlewaretoken'
        key_element = webpage.mainFrame().findAllElements(key_query).at(0)
        key = key_element.attribute('value')
        print '(%s,%s): %s'%(lSpecies,rSpecies,key)
        # Then 
        #response = urllib2.urlopen(loc, dat).read()

c = crawler()
c.main()
appInstance.exec_()
# code.interact('', None, locals())
