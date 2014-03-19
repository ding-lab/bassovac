#!/usr/bin/env python

import os
import os.path
import shutil
import sys
import tempfile
import unittest

class TestExpectedOutput(unittest.TestCase):
    def setUp(self):
        self.bvacExecutable = os.getenv("ITEST_EXECUTABLE")
        if self.bvacExecutable == None:
            raise ValueError("ITEST_EXECUTABLE unset")
        if os.access(self.bvacExecutable, os.X_OK) == False:
            raise ValueError("bassovac executable %s not executable" %(self.bvacExecutable))

        self.tmpDir = tempfile.mkdtemp()

        self.testDir = os.path.dirname(sys.argv[0])
        self.dataDir = os.path.join(self.testDir, "data")
        self.normalBam = os.path.join(self.dataDir, "n.bam")
        self.tumorBam  = os.path.join(self.dataDir, "t.bam")
        self.refSeq = "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa"

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def runBassovac(self, args):
        cmdline = "%s %s" %(self.bvacExecutable, " ".join(args))
        return os.system(cmdline) 

    def test_expected1(self):
        outputFile = os.path.join(self.tmpDir, "actual.out")
        expectedFile = os.path.join(self.dataDir, "expected.out")
        rv = self.runBassovac(
            (
                "-f", self.refSeq,
                "-n", self.normalBam,
                "-t", self.tumorBam,
                "-o", outputFile,
                "-x -q 1",
                "--normal-purity 1",
                "--tumor-purity 0.76"
            )
        )
        self.assertEqual(0, rv)
        expectedData = open(expectedFile).read()
        actualData   = open(outputFile).read()

        self.assertEqual(expectedData, actualData)

if __name__ == "__main__":
    unittest.main()
