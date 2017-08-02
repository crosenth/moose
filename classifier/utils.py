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
common utilities
"""

import bz2
import gzip
import itertools
import logging
import os
import sys


log = logging.getLogger(__name__)

ALIGNMENT_DTYPES = {
    'qseqid': str,
    'sseqid': str,
    'pident': float,
    'mismatch': float,
    'qstart': int,
    'qend': int,
    'qlen': int,
    'gapopen': int,
    'sstart': int,
    'send': int,
    'evalue': float,
    'bitscore': float,
    'length': int,
    'qcovs': float,
    'specimen': str,
    'assignment_id': str,
    'tax_name': str,
    'rank': str,
    'assignment_tax_name': str,
    'assignment_rank': str,
    'pident': float,
    'tax_id': str,
    'assignment_tax_id': str,
    'condensed_id': str,
    'accession': str,
    'qseqid': str}


def get_compression(io):
    if io is sys.stdout:
        compression = None
    else:
        compress_ops = {'.gz': 'gzip', '.bz2': 'bz2'}
        ext = os.path.splitext(io)[-1]
        compression = compress_ops.get(ext, None)
    return compression


def opener(mode='rt'):
    """Factory for creating file objects

    Keyword Arguments:
        - mode -- A string indicating how the file is to be opened. Accepts the
            same values as the builtin open() function.
        - bufsize -- The file's desired buffer size. Accepts the same values as
            the builtin open() function.
    """

    def open_file(f):
        if f is sys.stdout or f is sys.stdin:
            return f
        elif f == '-':
            return sys.stdin if 'r' in mode else sys.stdout
        elif f.endswith('.bz2'):
            return bz2.BZ2File(f, mode)
        elif f.endswith('.gz'):
            return gzip.open(f, mode)
        else:
            return open(f, mode)

    return open_file


def groupbyl(li, key=None, as_dict=False):
    groups = sorted(li, key=key)
    groups = itertools.groupby(groups, key=key)
    groups = ((g, list(l)) for g, l in groups)
    if as_dict:
        return(dict(groups))
    else:
        return groups
