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
Calculate a global qcovs given columns qstart, qend and qlen
"""

import sys
import logging
import pandas

from classifier import utils

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'blast', help='tabular blast file of query and subject hits')

    blast_parser = parser.add_argument_group('blast input options')
    blast_parser.add_argument(
        '--columns',
        help=('if no header specify columns names'))
    blast_parser.add_argument(
        '--tsv', action='store_true', help='default: csv')

    outs_parser = parser.add_argument_group('output options')
    outs_parser.add_argument(
        '-o', '--out',
        default=sys.stdout,
        metavar='FILE',
        help="classification results [default: stdout]")


def action(args):
    blast = pandas.read_csv(
        args.blast,
        dtype=utils.ALIGNMENT_DTYPES,
        names=args.columns.split(',') if args.columns else None,
        sep='\t' if args.tsv else ',')

    if not all(q in blast.columns for q in ['qstart', 'qend', 'qlen']):
        sys.exit('blast input must have columns [qstart qend qlen]')

    blast.loc[:, 'qcovs'] = (
        (blast['qend'] - blast['qstart'] + 1) / blast['qlen'] * 100)

    blast.to_csv(
        args.out,
        compression=utils.get_compression(args.out),
        header=not bool(args.columns),
        index=False,
        float_format='%.2f',
        sep='\t' if args.tsv else ',')
