# This file is part of classifier
#
#    classifier is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    classifier is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with classifier.  If not, see <http://www.gnu.org/licenses/>.

"""
Filter blast results by qcovs or pident given minimum thresholds
"""

import sys
import logging
import pandas

from classifier import utils, sequtils

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'blast', help='tabular blast file of query and subject hits')
    parser.add_argument(
        '--min-qcovs',
        type=float,
        help=('miminum coverage in blast results'))
    parser.add_argument(
        '--min-pident',
        type=float,
        help=('miminum coverage in blast results'))
    parser.add_argument(
        '--max-pident',
        type=float,
        help=('miminum coverage in blast results'))
    parser.add_argument(
        '--limit',
        help=('take head number of lines'))

    header_parser = parser.add_argument_group(
        title='alignment input header-less options',
        description=('assumed comma-seperated with header if not specified'))
    columns_parser = header_parser.add_mutually_exclusive_group(required=False)
    columns_parser.add_argument(
        '--blast6in', '-b',
        action='store_true',
        help=('header-less blast-like tab-separated input'))
    columns_parser.add_argument(
        '--columns', '-c',
        help=('specify columns for header-less comma-seperated values'))

    outs_parser = parser.add_argument_group('output options')
    outs_parser.add_argument(
        '-o', '--out',
        default=sys.stdout,
        metavar='FILE',
        help="classification results [default: stdout]")


def raw_filtering(blast_results, min_coverage=None,
                  max_identity=None, min_identity=None):
    """run raw hi, low and coverage filters and output log information
    """

    blast_results_len = len(blast_results)

    if min_coverage:
        # run raw hi, low and coverage filters
        blast_results = blast_results[
            blast_results['qcovs'] >= min_coverage]

        blast_results_post_len = len(blast_results)

        len_diff = blast_results_len - blast_results_post_len
        if len_diff:
            log.warn('dropping {} sequences below '
                     'coverage threshold'.format(len_diff))

        blast_results_len = blast_results_post_len

    if max_identity:
        blast_results = blast_results[
            blast_results['pident'] <= max_identity]

        blast_results_post_len = len(blast_results)

        len_diff = blast_results_len - blast_results_post_len
        if len_diff:
            log.warn('dropping {} sequences above max_identity'.format(
                len_diff))

        blast_results_len = blast_results_post_len

    if min_identity:
        blast_results = blast_results[
            blast_results['pident'] >= min_identity]

        blast_results_post_len = len(blast_results)

        len_diff = blast_results_len - blast_results_post_len
        if len_diff:
            log.warn('dropping {} sequences below min_identity'.format(
                len_diff))

        blast_results_len = blast_results_post_len

    return blast_results


def action(args):
    log.info('loading alignment file')
    if args.blast6in:
        blast = pandas.read_csv(
            args.blast,
            sep='\t',
            names=['qseqid,sseqid,pident,length,mismatch,gapopen,'
                   'qstart,qend,sstart,send,evalue,bitscore'],
            dtype=sequtils.DTYPES)
    else:
        blast = pandas.read_csv(
            args.blast,
            dtype=sequtils.DTYPES,
            names=args.columns.split(',') if args.columns else None)

    if args.min_qcovs:
        blast_len = len(blast)
        blast = blast[blast['qcovs'] >= args.min_qcovs]
        msg = 'removed {} hits under --min_qcovs {}'.format(
            blast_len - len(blast), args.min_qcovs)
        log.info(msg)

    if args.min_pident:
        blast_len = len(blast)
        blast = blast[blast['pident'] >= args.min_pident]
        msg = 'removed {} hits under --min_pident {}'.format(
            blast_len - len(blast), args.min_pident)
        log.info(msg)

    if args.max_pident:
        blast_len = len(blast)
        blast = blast[blast['pident'] <= args.max_pident]
        msg = 'removed {} hits over --max_pident {}'.format(
            blast_len - len(blast), args.max_pident)
        log.info(msg)

    blast.to_csv(
        args.out,
        compression=utils.get_compression(args.out),
        sep='\t' if args.blast6in else ',',
        header=not (args.blast6in or bool(args.columns)),
        float_format='%.2f',
        index=False)
