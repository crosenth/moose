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

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'blast', help='tabular blast file of query and subject hits')
    parser.add_argument(
        '--min-qcovs',
        help=('miminum coverage in blast results'))
    parser.add_argument(
        '--min-pident',
        help=('miminum coverage in blast results'))
    parser.add_argument(
        '--max-pident',
        help=('miminum coverage in blast results'))
    parser.add_argument(
        '--limit',
        help=('take head number of lines'))

    blast_parser = parser.add_argument_group('blast input options')
    blast_parser.add_argument(
        '--columns',
        help=('if no header specify columns names'))
    blast_parser.add_argument(
        '--csv', action='store_true', help='default: tabular')

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
    blast = pandas.read_csv(
        args.blast,
        dtype={'qseqid': str,
               'sseqid': str,
               'pident': float,
               'qcovs': float,
               'mismatch': float,
               'qstart': int,
               'qlen': int,
               'qend': int,
               'qcovs': float},
        names=args.columns,
        sep=',' if args.csv else '\t',
        limit=args.limit)

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
        header=not bool(args.columns),
        index=False,
        sep=',' if args.csv else '\t')
