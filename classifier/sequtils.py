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

import itertools
import logging
import utils

from operator import itemgetter

log = logging.getLogger(__name__)

RANKS = ['root', 'superkingdom', 'phylum', 'class',
         'order', 'family', 'genus', 'species']


DTYPES = {
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
    'bitscore': int,
    'length': int,
    'qcovs': float}


def condense_ids(assignments,
                 taxonomy,
                 ranks=RANKS,
                 floor_rank=None,
                 ceiling_rank=None,
                 max_size=3,
                 rank_thresholds={}):
    """
    assignments = [tax_ids...]
    taxonomy = {taxid:taxonomy}

    Functionality: Group items into taxonomic groups given max rank sizes.
    """

    floor_rank = floor_rank or ranks[-1]
    ceiling_rank = ceiling_rank or ranks[0]

    if not taxonomy:
        raise TypeError('taxonomy must not be empty or NoneType')

    if floor_rank not in ranks:
        msg = '{} not in ranks: {}'.format(floor_rank, ranks)
        raise TypeError(msg)

    if ceiling_rank not in ranks:
        msg = '{} not in ranks: {}'.format(ceiling_rank, ranks)
        raise TypeError(msg)

    if ranks.index(floor_rank) < ranks.index(ceiling_rank):
        msg = '{} cannot be lower rank than {}'.format(
            ceiling_rank, floor_rank)
        raise TypeError(msg)

    # set rank to ceiling
    try:
        assignments = {a: taxonomy[a][ceiling_rank] for a in assignments}
    except KeyError:
        print assignments
        error = ('Assignment id not found in taxonomy.')
        raise KeyError(error)

    def condense(groups, ceiling_rank=ceiling_rank, max_size=max_size):
        new_groups = {}

        for a, r in groups.items():
            new_groups[a] = taxonomy[a][ceiling_rank] or r

        num_groups = len(set(new_groups.values()))

        if rank_thresholds.get(ceiling_rank, max_size) < num_groups:
            return groups

        groups = new_groups

        # return if we hit the floor
        if ceiling_rank == floor_rank:
            return groups

        # else move down a rank
        ceiling_rank = ranks[ranks.index(ceiling_rank) + 1]

        # recurse each branch down the tax tree
        for _, g in utils.groupbyl(groups.items(), itemgetter(1)):
            g = condense(dict(g),
                         ceiling_rank,
                         max_size - num_groups + 1)
            groups.update(g)

        return groups

    return condense(assignments)


def compound_assignment(assignments, taxonomy):
    """
    Create taxonomic names based on 'assignmnets', which are a set of
    two-tuples: {(tax_id, is_starred), ...} where each tax_id is a key
    into taxdict, and is_starred is a boolean indicating whether at
    least one reference sequence had a parirwise alignment identity
    score meeting some thresholed. 'taxdict' is a dictionary keyed by
    tax_id and returning a dict of taxonomic data corresponding to a
    row from the taxonomy file. If 'include_stars' is False, ignore
    the second element of each tuple in 'assignments' and do not
    include asterisks in the output names.
    assignments = [(tax_id, is_starred),...]
    taxonomy = {taxid:taxonomy}
    Functionality: see format_taxonomy
    """

    if not taxonomy:
        raise TypeError('taxonomy must not be empty or NoneType')

    assignments = ((taxonomy[i]['tax_name'], a) for i, a in assignments)
    assignments = zip(*assignments)

    return format_taxonomy(*assignments, asterisk='*')


def format_taxonomy(names, selectors, asterisk='*'):
    """
    Create a friendly formatted string of taxonomy names. Names will
    have an asterisk value appended *only* if the cooresponding
    element in the selectors evaluates to True.
    """

    names = itertools.izip_longest(names, selectors)
    names = ((n, asterisk if s else '')
             for n, s in names)  # add asterisk to selected names
    names = set(names)
    names = sorted(names)  # sort by the name plus asterisk
    names = itertools.groupby(names, key=itemgetter(0))  # group by just the names
    # prefer asterisk names which will be at the bottom
    names = (list(g)[-1] for _, g in names)
    names = (n + a for n, a in names)  # combine names with asterisks
    # assume species names have exactly two words

    def is_species(s):
        return len(s.split()) == 2

    names = sorted(names, key=is_species)
    names = itertools.groupby(names, key=is_species)

    tax = []

    for species, assigns in names:
        if species:
            # take the species word and combine them with a '/'
            assigns = (a.split() for a in assigns)
            # group by genus name
            assigns = itertools.groupby(assigns, key=itemgetter(0))
            assigns = ((k, map(itemgetter(1), g))
                       for k, g in assigns)  # get a list of the species names
            assigns = ('{} {}'.format(k, '/'.join(g))
                       for k, g in assigns)  # combine species names with '/'

        tax.extend(assigns)

    return ';'.join(sorted(tax))
