"""
Test classifier
"""

import bz2
import csv
import filecmp
import logging
import os
import sys

from classifier import main

from tests import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)


class TestClassify(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        main(['classify'] + [str(a) for a in arguments])

    log_info = 'classifier classify {}'

    copy_numbers = os.path.join(datadir, 'rrnDB_16S_copy_num.csv.bz2')

    thisdatadir = os.path.join(datadir, 'TestClassify')

    def test01(self):
        """
        Minimal inputs.
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = os.path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = os.path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--out', classify_out,
            '--details-out', details_out,
            '--seq-info', seq_info,
            blast,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test16(self):
        """
        test no blast result
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = os.path.join(thisdatadir, 'blast.csv.bz2')
        specimen_map = os.path.join(thisdatadir, 'map_single.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = os.path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--specimen-map', specimen_map,
            '--out', classify_out,
            '--details-out', details_out,
            '--seq-info', seq_info,
            blast,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test02(self):
        """
        Include weights.
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        weights = os.path.join(thisdatadir, 'weights.csv.bz2')
        taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = os.path.join(thisdatadir, 'blast.csv.bz2')
        specimen_map = os.path.join(thisdatadir, 'map_single.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = os.path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--specimen-map', specimen_map,
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--weights', weights,
            '--out', classify_out,
            '--details-out', details_out,
            '--seq-info', seq_info,
            blast,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test03(self):
        """
        Include specimen-map.
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        specimen_map = os.path.join(thisdatadir, 'map.csv.bz2')
        taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = os.path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = os.path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--specimen-map', specimen_map,
            '--out', classify_out,
            '--details-out', details_out,
            '--seq-info', seq_info,
            blast,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test06(self):
        """
        All together
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        weights = os.path.join(thisdatadir, 'weights.csv.bz2')
        specimen_map = os.path.join(thisdatadir, 'map.csv.bz2')
        taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = os.path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = os.path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--specimen-map', specimen_map,
            '--weights', weights,
            '--copy-numbers', self.copy_numbers,
            '--out', classify_out,
            '--details-out', details_out,
            '--seq-info', seq_info,
            blast,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test07(self):
        """
        Test validation of type strains
        """

        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        blast = os.path.join(thisdatadir, this_test, 'blast.csv.bz2')
        taxonomy = os.path.join(thisdatadir, this_test, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, this_test, 'seq_info.csv.bz2')
        specimen_map = os.path.join(thisdatadir, this_test, 'map.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = os.path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--specimen-map', specimen_map,
            '--seq-info', seq_info,
            '--details-out', details_out,
            '--out', classify_out,
            blast, taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test08(self):
        """
        Test empty blast_results file
        """

        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv')
        details_out = os.path.join(outdir, 'details.csv')

        # create blank blast.csv file
        blast = os.path.join(outdir, 'blast.csv')
        open(blast, 'w').close()
        taxonomy = os.path.join(thisdatadir, this_test, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, this_test, 'seq_info.csv.bz2')

        args = [
            '--seq-info', seq_info,
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--details-out', details_out,
            '--out', classify_out,
            blast, taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertFalse(os.path.isfile(classify_out))
        self.assertFalse(os.path.isfile(details_out))

    def test09(self):
        """
        Test all [no blast results] classification with hits
        """

        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        blast = os.path.join(thisdatadir, this_test, 'blast.csv.bz2')
        taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')
        weights = os.path.join(thisdatadir, 'weights.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = os.path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--seq-info', seq_info,
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--weights', weights,
            '--details-out', details_out,
            '--out', classify_out,
            blast, taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test10(self):
        """
        Parse non-default blast result files that have headers
        """

        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv')

        blast = os.path.join(thisdatadir, 'blast_extrafields.csv.bz2')
        specimen_map = os.path.join(thisdatadir, 'map_single.csv.bz2')
        taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')

        args = [
            '--seq-info', seq_info,
            '--specimen-map', specimen_map,
            '--out', classify_out,
            blast,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))

    def test11(self):
        """
        Test dynamic thresholding

        github issue #32
        """

        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        blast = os.path.join(thisdatadir, this_test, 'blast.csv.bz2')
        taxonomy = os.path.join(thisdatadir, this_test, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, this_test, 'seq_info.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = os.path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--seq-info', seq_info,
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--details-out', details_out,
            '--out', classify_out,
            blast, taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test13(self):
        """
        Test ordering of assignment_id

        Note: This test was inspired by a case when upgrading
        from numpy 1.9.2 to 1.10.1 where, with sort_index we were using
        an "unstable" sort on a field (specimen) with only one value, and
        our classifications became unsorted.  Data here is derived from that
        test, previously would have failed this test, and subsequently passes
        """
        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        blast = os.path.join(thisdatadir, this_test, 'blast.csv.bz2')
        taxonomy = os.path.join(datadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, this_test, 'seq_info.csv.bz2')
        weights = os.path.join(thisdatadir, this_test, 'weights.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')

        args = [
            '--seq-info', seq_info,
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--specimen', 'specimen',
            '--weights', weights,
            '--out', classify_out,
            blast, taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))

    def test14(self):
        """
        Test --hits-below-threshold
        """
        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        blast = os.path.join(thisdatadir, 'blast.csv.bz2')
        taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')
        specimen_map = os.path.join(thisdatadir, 'map_single.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = os.path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--seq-info', seq_info,
            '--specimen-map', specimen_map,
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--hits-below-threshold',
            '--details-out', details_out,
            '--out', classify_out,
            blast, taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test15(self):
        """
        Test --best-n-hits
        """
        thisdatadir = self.thisdatadir

        # Blast results contain:
        # 2 strep mutans, 0 mismatch
        # 1 strep troglodytae, 2 mismatch
        # 1 strep infantarius, 3 mismatch
        # All (artificially) > 99%
        blast = os.path.join(thisdatadir, 'blast_3strepto.csv')
        taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        args = [
            '--seq-info', seq_info,
            '--best-n-hits', 3,
            '--details-out', details_out,
            '--out', classify_out,
            blast, taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        names = set()
        for line in csv.DictReader(bz2.open(details_out, mode='rt')):
            names.add(line['tax_name'])

        """
        Normally we would expect 4 details rows spanning 3 tax_names, lending
        to the classification "Streptococcus infantarius/mutans*/troglodytae
        With --best-n-hits, we expect 3 details rows
        lending to the classification
        """
        self.assertTrue(len(names) == 2)

    def test17(self):
            """
            min-pident 99, max-pident 100
            """

            this_test = sys._getframe().f_code.co_name

            thisdatadir = self.thisdatadir

            taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
            seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')
            blast = os.path.join(thisdatadir, 'blast.csv.bz2')

            outdir = self.mkoutdir()

            classify_out = os.path.join(outdir, 'classifications.csv.bz2')
            details_out = os.path.join(outdir, 'details.csv.bz2')

            classify_ref = os.path.join(
                thisdatadir, this_test, 'classifications.csv.bz2')
            details_ref = os.path.join(
                thisdatadir, this_test, 'details.csv.bz2')

            args = [
                '--seq-info', seq_info,
                '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
                '--max-pident', '100',
                '--min-pident', '99',
                '--out', classify_out,
                '--details-out', details_out,
                blast,
                taxonomy]

            log.info(self.log_info.format(' '.join(map(str, args))))

            self.main(args)

            self.assertTrue(filecmp.cmp(classify_ref, classify_out))
            self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test18(self):
        """
        min qcovs
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        taxonomy = os.path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = os.path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = os.path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = os.path.join(outdir, 'classifications.csv.bz2')
        details_out = os.path.join(outdir, 'details.csv.bz2')

        classify_ref = os.path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = os.path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--seq-info', seq_info,
            '--columns', 'qseqid,sseqid,pident,qstart,qend,qlen,qcovs',
            '--min-qcovs', '99',
            '--out', classify_out,
            '--details-out', details_out,
            blast,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))
