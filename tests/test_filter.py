"""
Test classifier
"""

import logging
import os

import filecmp
import sys

from classifier import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)


class TestFilter(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        main(['filter'] + [str(a) for a in arguments])

    log_info = 'classifier filter {}'

    thisdatadir = os.path.join(datadir, 'TestFilter')

    def test01(self):
        """
        Minimal inputs.
        """

        thisdatadir = self.thisdatadir

        blast = os.path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        filter_out = os.path.join(outdir, 'blast.csv.bz2')

        args = [
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--out', filter_out,
            blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(blast, filter_out))

    def test02(self):
        """
        min-pident 99, max-pident 100
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        blast = os.path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        filter_out = os.path.join(outdir, 'blast.csv.bz2')

        filter_ref = os.path.join(
            thisdatadir, this_test, 'blast.csv.bz2')

        args = [
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--max-pident', '100',
            '--min-pident', '99',
            '--out', filter_out,
            blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(filter_ref, filter_out))

    def test03(self):
        """
        min qcovs
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        blast = os.path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        filter_out = os.path.join(outdir, 'blast.csv.bz2')

        filter_ref = os.path.join(
            thisdatadir, this_test, 'blast.csv.bz2')

        args = [
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--min-qcovs', '99',
            '--out', filter_out,
            blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(filter_ref, filter_out))

    def test04(self):
        """
        all together
        """

        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        outdir = self.mkoutdir()

        blast = os.path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        filter_out = os.path.join(outdir, 'blast.csv.bz2')

        filter_ref = os.path.join(
            thisdatadir, this_test, 'blast.csv.bz2')

        args = [
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--min-pident', '98',
            '--max-pident', '99',
            '--min-qcov', '99',
            '--out', filter_out,
            blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(filter_ref, filter_out))
