import sys, code
from PyQt4.QtCore import *
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
        print self.organismNames
        exit(0)

c = crawler()
c.main()
appInstance.exec_()
# code.interact('', None, locals())
