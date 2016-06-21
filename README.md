# TODO

1. Document features

```
positional arguments:
  blast                 tabular blast file of query and subject hits
  seq_info              File mapping reference seq name to tax_id
  taxonomy              Table defining the taxonomy for each tax_id

optional arguments:
  -h, --help            show this help message and exit
  --threads NUM         Number of threads (CPUs). Can also specify with
                        environment variable THREADS_ALLOC. [28]

blast input options:
  --has-header          if blast data has a header
  --columns COLUMNS     column specifiers. global query coverage can be
                        calculated by including the qstart, qend and qlen
                        specifiers. column "mismatch will filter out all but
                        the best N hits based on the number of mismatches.
  --tab                 default: csv

blast filtering options:
  --limit LIMIT         limit number of blast results
  --min-pident PERCENT  minimum identity threshold for accepting matches
  --max-pident PERCENT  maximum identity threshold for accepting matches
  --min-cluster-size INTEGER
                        minimum cluster size to include in classification
                        output [1]
  --min-qcovs PERCENT   percent of alignment coverage of blast result
  --best-n-hits BEST_N_HITS
                        For each query sequence, filter out all but the best N
                        hits. Used in conjunction with blast "mismatch"
                        column.

assignment options:
  --starred PERCENT     Names of organisms for which at least one reference
                        sequence has pairwise identity with a query sequence
                        of at least PERCENT will be marked with an asterisk
                        [100.0]
  --max-group-size INTEGER
                        group multiple target-rank assignments that excede a
                        threshold to a higher rank [3]
  --split-condensed-assignments
                        Do not combine condensed identical assignments

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
                        described in the blast input (used, for example, to
                        describe cluster sizes for corresponding cluster
                        centroids).

output options:
  --details-full        do not limit out_details to only larget cluster per
                        assignment
  --include-ref-rank INCLUDE_REF_RANK
                        Given a single rank (species,genus,etc), include each
                        reference sequence's tax_id as $\{rank\}_id and its
                        taxonomic name as $\{rank\}_name in details output
  --hits-below-threshold
                        Hits that were below the best-rank threshold will be
                        included in the details
  -O FILE, --details-out FILE
                        Optional details of taxonomic assignments.
  -o FILE, --out FILE   classification results [default: stdout]
```

1. Decide what features are needed
    * calculate query coverage
    * filter Blast/(aligner) hits (everythin gunder blast filtering options)
    * expand results using weights file
    * required columns [qseqid,sseqid,pident,qstart,qend,qlen]
        * change to [qseqid,sseqid,pident]
        * optional qcovs
    *  pplacer outputs
1. Standalone organization
    * subcommands
    * executables
    * command arguments
