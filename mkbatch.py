#!/usr/bin/env python
import os, sys, glob, cjson, re
if len(sys.argv) != 2:
    print "Usage:"
    print "../mkbatch.py file_number_%i.ext > ../fnames"
    print "../mkbatch.py -restore < ../fnames"
    exit(0)
if sys.argv[1] == '-restore':
    reverse_map = cjson.decode(sys.stdin.read())
    for fn,newn in reverse_map.iteritems():
        os.rename(newn,fn)
else:
    fnames = os.listdir('.')
    reverse_map = {}
    i=0
    for fn in fnames:
        i += 1
        newn = re.sub('%i',str(i),sys.argv[1])
        reverse_map[fn]=newn
        os.rename(fn,newn)
    print cjson.encode(reverse_map)

