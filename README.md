Groups pairwise alignments by taxonomy and alignment scores.  Works safely 
with large data sets using the Python Data Analysis Library.

Note - This package requires a fully constructed taxonomy table generated
by Taxtastic (https://github.com/fhcrc/taxtastic).  The taxtable must reflect
the taxonomy of the refset being aligned against.

```
positional arguments:
  alignments            alignment file with query and subject sequence hits
                        and optional header
  seq_info              File mapping reference seq name to tax_id
  taxonomy              Table defining the taxonomy for each tax_id

optional arguments:
  -h, --help            show this help message and exit

alignment input header-less options:
  assumed comma-seperated with header if not specified

  --blast6in, -b        header-less blast-like tab-separated input
  --columns COLUMNS, -c COLUMNS
                        specify columns for header-less comma-seperated values

filtering options:
  --best-n-hits BEST_N_HITS
                        For each qseqid sequence, filter out all but the best
                        N hits. Used in conjunction with alignment "mismatch"
                        column.
  --max-pident MAX_PIDENT
                        miminum coverage of aligments
  --min-cluster-size INTEGER
                        minimum cluster size to include in classification
                        output [1]
  --min-pident MIN_PIDENT
                        miminum coverage of alignments
  --min-qcovs MIN_QCOVS
                        miminum coverage of alignments [None]

assignment options:
  --starred PERCENT     Names of organisms for which at least one reference
                        sequence has pairwise identity with a query sequence
                        of at least PERCENT will be marked with an asterisk
                        [100.0]
  --max-group-size INTEGER
                        group multiple target-rank assignments that excede a
                        threshold to a higher rank [3]
  --split-condensed-assignments
                        Group final assignment classifications before
                        assigning condensed taxonomic ids

other input options:
  --copy-numbers CSV    Estimated 16s rRNA gene copy number for each tax_ids
                        (CSV file with columns: tax_id, median)
  --rank-thresholds CSV
                        Columns [tax_id,ranks...]
  --specimen LABEL      Single group label for reads
  --specimen-map CSV    CSV file with columns (name, specimen) assigning
                        sequences to groups.
  -w CSV, --weights CSV
                        Optional headless csv file with columns 'seqname',
                        'count' providing weights for each query sequence
                        described in the alignment input (used, for example,
                        to describe cluster sizes for corresponding cluster
                        OTUs).

output options:
  --details-full        do not limit out_details to largest cluster per
                        assignment [False]
  --include-ref-rank INCLUDE_REF_RANK
                        Given a single rank (species,genus,etc), include each
                        reference sequence's tax_id as $\{rank\}_id and its
                        taxonomic name as $\{rank\}_name in details output
  --hits-below-threshold
                        Hits that were below the best-rank threshold will be
                        included in the details
  --details-out FILE    Optional details of taxonomic assignments.
  -o FILE, --out FILE   classification results [default: stdout]
```
