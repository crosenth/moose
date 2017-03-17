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
Provides taxonomic lineages for each classifier assignment per assignment_id
"""

import sys
import logging
import pandas

from classifier import utils

log = logging.getLogger(__name__)


def build_parser(parser):
    """
    build the parser
    """

    parser.add_argument(
        'assignments',
        help='output from classify')
    parser.add_argument(
        'details',
        help='with assignment_ids and assignment_tax_ids')
    parser.add_argument(
        'taxtable',
        help='with lineages')
    parser.add_argument(
        '--tax_ids',
        action='store_false',
        dest='names',
        help='')
    parser.add_argument(
        '-o', '--out',
        default=sys.stdout,
        metavar='FILE',
        help="classification results [default: stdout]")


#  Need this to identify rank columns in the taxtable
TAX_TABLE_NON_RANKS = ['tax_id', 'parent_id', 'rank', 'tax_name']


def action(args):
    assignments = pandas.read_csv(
        args.assignments,
        usecols=['specimen', 'assignment', 'assignment_id'],
        dtype=utils.ALIGNMENT_DTYPES)

    details = pandas.read_csv(
        args.details,
        usecols=['specimen', 'assignment_id', 'assignment_tax_id'],
        dtype=utils.ALIGNMENT_DTYPES)
    details = details.drop_duplicates()

    taxtable = pandas.read_csv(
        args.taxtable,
        dtype=str)
    taxtable = taxtable.set_index('tax_id')

    lineages = assignments.merge(
        details, on=['specimen', 'assignment_id'], how='left')
    lineages = lineages.merge(
        taxtable, left_on='assignment_tax_id', right_index=True, how='left')
    lineages = lineages.dropna(axis=1, how='all')

    if args.names:
        ranks = [c for c in taxtable if c not in TAX_TABLE_NON_RANKS]
        lineage_ranks = [c for c in lineages.columns if c in ranks]
        for col in lineage_ranks:
            lineages = lineages.merge(
                taxtable[['tax_name']],
                left_on=col,
                right_index=True,
                how='left',
                suffixes=('', '_'))
            lineages = lineages.drop(col, axis=1)
            lineages = lineages.rename(columns={'tax_name_': col})

    lineages.to_csv(args.out, index=False)
