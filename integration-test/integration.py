#!/usr/bin/env python

import os
import sys
from getopt import getopt
import unittest 

class TestSuite(): 
    def __init__(self, tests):
        self.tests = tests

    def suite(self): 
        alltests = unittest.TestSuite()  
        for module in map(__import__, self.tests):  
            alltests.addTest(unittest.findTestCases(module))  
        return alltests  


if __name__ == "__main__":
    opts, args = getopt(sys.argv[1:], "", ['exe='])
    hashOpts = {}
    for kv in opts: hashOpts[kv[0]] = kv[1]
    os.environ["ITEST_EXECUTABLE"] = hashOpts['--exe']

    if len(args) == 0:
        print >> sys.stderr, "No test specified!"
        sys.exit(1)

    testSuite = TestSuite(args)
    print "RV", unittest.main(defaultTest="testSuite.suite", argv=[sys.argv[0]])
