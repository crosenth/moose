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

import logging
import utils

from operator import itemgetter

log = logging.getLogger(__name__)

RANKS = ['root', 'superkingdom', 'phylum', 'class',
         'order', 'family', 'genus', 'species']


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
