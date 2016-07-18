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
Pivot classifier output on read counts with
assignments as rows and specimens as columns
"""

import sys
import logging
import pandas

from classifier import utils

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'assignments', help='output from classify')
    parser.add_argument(
        '--details',
        help=('for tax_id column in output'))

    outs_parser = parser.add_argument_group('output options')
    outs_parser.add_argument(
        '-o', '--out',
        default=sys.stdout,
        metavar='FILE',
        help="classification results [default: stdout]")


def action(args):
    assignments = pandas.read_csv(
        args.assignments,
        dtype=utils.ALIGNMENT_DTYPES)

    rows = ['assignment', 'best_rank']

    if args.details:
        details = pandas.read_csv(
            args.details,
            dtype=utils.ALIGNMENT_DTYPES)

        assignments = assignments.merge(details, how='left')

        by = ['specimen', 'assignment_id']
        for (specimen, assignment_id), df in assignments.groupby(by=by):
            tax_ids = ','.join(df['condensed_id'].drop_duplicates().tolist())
            assignments.loc[df.index, 'tax_id'] = tax_ids

        rows.insert(1, 'tax_id')

    pivot = pandas.pivot_table(
        assignments,
        index=rows,
        columns='specimen',
        values='reads')
    pivot.index = pivot.index.rename('rank', level='best_rank')

    pivot.to_csv(args.out, float_format='%.0f')
