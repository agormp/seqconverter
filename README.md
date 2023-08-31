# seqconverter

![](https://img.shields.io/badge/version-2.13.1
-blue)
[![PyPI downloads](https://static.pepy.tech/personalized-badge/seqconverter?period=total&units=none&left_color=black&right_color=blue&left_text=PyPI%20downloads&service=github)](https://pepy.tech/project/seqconverter)

The command-line program `seqconverter` can read and write text files containing aligned or unaligned DNA or protein sequences. The program understands most standard and some non-standard formats (fasta, Nexus, Phylip, Clustal, Stockholm, tab, raw, Genbank, How). The tool can be used to convert between sequence file formats, and is also able to perform various manipulations and analyses of sequences.

## Availability

The seqconverter source code is available on GitHub: https://github.com/agormp/seqconverter. The executable can be installed from PyPI: https://pypi.org/project/seqconverter/

## Installation

```
python3 -m pip install seqconverter
```

Upgrading to latest version:

```
python3 -m pip install --upgrade seqconverter
```

## Dependencies

seqconverter relies on the [sequencelib library](https://github.com/agormp/sequencelib) and the [NumPy package](https://numpy.org), which are automatically included when using pip to install.

## Highlights

* Can be used to convert between sequence file formats but also able to perform many other manipulations and analyses of sequences.
* Read and write **aligned** sequences in the following formats:
	* fasta
	* Nexus
	* Phylip
	* Clustal
	* Stockholm (so far only read)
	* tab
	* raw
* Read and write **unaligned** sequences in the following formats:
	* fasta
	* tab
	* raw
	* Genbank
	* How
* Writes to stdout, so output can be used in pipes or redirected to file
* Also accepts input on stdin
* Options to select or discard sequences based on one of several criteria: name matches regular expression, name in NAMEFILE, sequence contains  specific residues on specific positions, duplicate (identical) sequences, duplicate names, sequence has many gaps at ends (<=> is shorter than other sequences), random sample of given size, ...
* Options to select or remove columns from alignment based on one of several criteria: all gaps, some gaps, more than fraction gaps, more than fration endgaps, conserved, specified indices, ...
* Extract all overlapping windows of specified size
* Options to rename one or more sequences based on various criteria
* Options to merge identically named sequences from multiple sequence files (end-to-end or discarding automatically discovered overlaps)
* Options to automatically create Nexus charset commands based on merging multiple individual files (e.g., one charset/partition per gene).
* Can automatically write MrBayes block with template for commands to run partitioned analysis, also based on merging multiple separate sequence alignments.
* Can translate and find reverse complement for DNA sequences
* Options to obtain summary information about sequences and alignments: number of seqs, names, lengths, composition (overall or per sequence), nucleotide diversity (pi), site summary (how many columns are variable, contain multiple residues, contain gaps, or contain IUPAC ambiguity symbols)
* More...
* Underlying library has been optimized for high speed and low memory consumption
* Really has too many options, but does useful stuff (and has been created based on what I needed for own projects)

## Quick start usage examples

These examples highlight some of the options available. For the full list use option -h to get help.

### Get help:

```
seqconverter -h
```

----

### Convert aligned sequences in fasta format to nexus, 70 characters per line

Note: output is written to the terminal so you need to use redirection to store in a file.

```
seqconverter -I fasta -O nexus --width 70 myalignment.fasta > myalignment.nexus
```

----

### Select all sequences whose name match the regular expression "seq_1[0-9]+"

```
seqconverter -I fasta -O fasta --select "seq_1[0-9]+" myseqs.fasta > subset.fasta
```

----

### Discard all sequences whose name match the regular expression "seq_1[0-9]+":

```
seqconverter -I fasta -O fasta --discard "seq_1[0-9]+" myseqs.fasta > subset.fasta
```

----

### Select random subset of 50 sequences from input file

```
seqconverter -I fasta -O fasta --subsample 50 myseqs.fasta > subset.fasta
```

----

### Select all sequences containing a Lysine at position 484 and a Tyrosine at position 501 

```
seqconverter -I clustal -O fasta --filterpos 484K,501Y myalignment.aln > voc.fasta
```

----

### Select columns 50-150 from ClustalW formatted alignment file, write output in fasta

This command extracts columns 50-150 (inclusive, with numbering starting at 1) from alignment in Clustal format, and writes output in fasta format to file (using redirection):

```
seqconverter -I clustal -O fasta --subseq 50,150 myalignment.aln > aligment_50_150.fasta
```

----

### Remove columns, where one or more residues are gaps, from alignment:

This command will remove alignment columns if any sequence has a gap in that position.

```
seqconverter -I fasta -O fasta --remgapcols myalignment.fasta > gapfree.fasta
```

----

### Remove columns, where more than 75% have endgaps, from alignment:

This command will remove alignment columns if more than 75% of sequences have endgaps in that position. An endgap is defined as a contiguous gappy region at either the beginning or end of a sequence, and are often a result of missing data (the gaps then do not represent insertion or deletion events).

```
seqconverter -I fasta -O fasta --remendgapcols 0.75 myalignment.fasta > endgapfree.fasta
```

----

### Concatenate identically named sequences from separate input files:

Sequences are pasted end to end in the same order as the input files. All input files must contain same number of sequences, and sequences in different files must have same name (for instance each file could contain an alignment of the sequences for a specific gene from a number of different species, and each sequence could then have the name of the species). The order of sequences in different files does not matter.

When used with the --charset (and possibly --mbpartblock) option this can be used to set up a partitioned analysis in MrBayes or BEAST (see below).

```
seqconverter -I fasta -O fasta --paste gene1.fasta gene2.fasta gene3.fasta > concat.fasta
```

----

### Concatenate sequences from multiple files, create partitioned Nexus file containing charset command

This command concatenates identically named sequences from separate input alignments, creating a partitioned Nexus file with `charset` specification. Start and stop indices for different charsets are automatically derived from lengths of sub-alignments. Charsets are named based on the names of included files.

This can be used for phylogenetic analyses in BEAST or MrBayes where different genomic regions (e.g., genes) have different substitution models. Note: sequences in each file need to have identical names (e.g. name of species).

```
seqconverter -I fasta -O nexus --paste --charset \
                               gene1.fasta gene2.fasta gene3.fasta > partitioned.nexus
```

----

### Concatenate sequences from multiple files, create partitioned Nexus file with commands to run MrBayes or BEAST analysis

This command does the same as the example above, and additionally adds a MrBayes block containing commands to run a partitioned analysis. The commands have sensible default values (e.g., setting DNA substution models to "nst=mixed" and unlinking most parameters across partitions). Optimally the commands should be tweaked according to the concrete data set. Importing the Nexus file in BEAUTI should result in setting most corresponding options for a BEAST run (but check, and remember to set priors etc.)

```
seqconverter -I fasta -O nexus --paste --charset --mbpartblock \
                               gene1.fasta gene2.fasta gene3.fasta > partitioned.nexus
```


## Usage

```
usage: seqconverter [-h] [-I FORMAT] [-O FORMAT] [--width WIDTH] [--subsample N]
                    [--select "REGEXP"] [--discard "REGEXP"] [--subset NAMEFILE]
                    [--remseqs NAMEFILE] [--filterpos VARIANT[,VARIANT,...]] [--filterdupseq]
                    [--filterdupname] [--remendgapseqs N] [--subseq START,STOP]
                    [--subseqrename] [--windows WSIZE] [--degap] [--remcols INDEX-LIST]
                    [--remambigcols] [--remgapcols] [--remallgapcols] [--remfracgapcols FRAC]
                    [--remendgapcols FRAC] [--remconscols] [--remhmminsertcols]
                    [--rename OLD,NEW] [--renamenumber BASENAME] [--appendnumber]
                    [--renameregexp "REGEXP"] [--regdupfix] [--savenames FILE]
                    [--restorenames FILE] [--gbname FIELD1[,FIELD2,FIELD3,...]] [--paste]
                    [--overlap] [--minoverlap N] [--multifile] [--charset] [--mbpartblock]
                    [--revcomp] [--translate] [--num] [--len] [--com] [--seqcom]
                    [--ignoregaps] [--nam] [--div] [--sit] [--debug]
                    [SEQFILE ...]

positional arguments:
  SEQFILE               One or more sequence files

options:
  -h, --help            show this help message and exit
  --debug               Print longer error messages

File formats:
  -I FORMAT             Input format: auto, fasta, nexus, phylip, clustal, stockholm,
                        genbank, tab, raw, how [default: auto]
  -O FORMAT             Output format: fasta, nexus, nexusgap, phylip, clustal, tab, raw, how
                        [default: fasta]
  --width WIDTH         Print sequences with WIDTH characters per line [default: 60]

Retrieve subset of sequences:
  --subsample N         Randomly extract N sequences from sequence set
  --select "REGEXP"     Select sequences where substring of name matches regular expression
  --discard "REGEXP"    Discard sequences where substring of name matches regular expression
  --subset NAMEFILE     Retrieve sequences listed in NAMEFILE
  --remseqs NAMEFILE    Discard sequences listed in NAMEFILE
  --filterpos VARIANT[,VARIANT,...]
                        Retrieve sequences containing specific residues on specific
                        positions. Syntax is: <POS><RESIDUE>, possibly in a comma-separated
                        list. Example: 484K,501Y
  --filterdupseq        Remove duplicate sequences (keeping one of each); print names of
                        removed sequences on stderr.
  --filterdupname       Remove sequences with duplicate names (keeping one of each). If this
                        option is not set (default): stop execution on duplicate names.
  --remendgapseqs N     Discard sequences with endgaps longer than N positions

Extracting or removing parts of sequences:
  --subseq START,STOP   Extract subsequence, positions START to STOP, from alignment
  --subseqrename        When extracting sub-sequences: add '_START_STOP' to seqnames
  --windows WSIZE       Extract all overlapping sequence windows of size WSIZE
  --degap               Remove all gap characters from sequences
  --remcols INDEX-LIST  Remove listed columns from alignment. Columns can be indicated as
                        comma-separated list of indices, and as ranges. Example:
                        --remcols=10,15,22-40,57
  --remambigcols        Remove columns where one or more residues are ambiguity symbols
                        (e.g., N for nucleotides)
  --remgapcols          Remove columns where one or more residues are gaps
  --remallgapcols       Remove columns that are all-gaps
  --remfracgapcols FRAC
                        Remove columns that contain > FRAC fraction gaps
  --remendgapcols FRAC  Remove columns where > FRAC fraction have endgaps
  --remconscols         Remove conserved columns from alignment
  --remhmminsertcols    When reading Stockholm format file from HMMer's hmmalign: remove
                        columns corresponding to insert states

Renaming sequences:
  --rename OLD,NEW      Rename single sequence from OLD to NEW
  --renamenumber BASENAME
                        Rename all sequences to this form: BASENAME_001, ...
  --appendnumber        Append numbering at end of existing sequence names (SeqA_001,
                        SeqXYZ_002, ...
  --renameregexp "REGEXP"
                        Rename sequences by deleting parts of names matching regular
                        expression in REGEXP
  --regdupfix           Fix duplicate names, created by regexp, by appending numbers to
                        duplicates (seqA, seqA_2, ...)
  --savenames FILE      Save renaming information in FILE for later use
  --restorenames FILE   Restore original names using info previously saved in FILE
  --gbname FIELD1[,FIELD2,FIELD3,...]
                        For Genbank input: construct sequence names from the list of named
                        fields, in the specified order

Combining multiple sequence files:
  --paste               Concatenate identically named sequences from separate input files.
                        Sequences are pasted end to end in the same order as the input files.
                        All input files must contain same number of sequences, and sequences
                        in different files must have same name.(To see partitions choose
                        nexus output, or output to multiple partition files).
  --overlap             Similar to --paste, but for input alignments that overlap partly.
                        Overlap is discovered automatically and partition boundaries are then
                        set such that each partition is covered by a unique set of genes. (To
                        see partitions choose nexus output, or output to multiple partition
                        files).
  --minoverlap N        Minimum overlap required for merging input alignments (default: set
                        automatically based on seq lengths)
  --multifile           Outputs to multiple files (one per partition) instead of stdout.
                        Partitions are generated automatically based on other options.
  --charset             Appends Nexus form charset block listing partitions in data (forces
                        output in Nexus format). Charsets and partitions are generated
                        automatically based on other options.
  --mbpartblock         Appends MrBayes block with commands for running partitioned analysis
                        (forces output in Nexus format). Charsets and partitions are
                        generated automatically based on other options.

DNA manipulations:
  --revcomp             Return reverse complement of sequence(s). Requires sequences to be
                        DNA.
  --translate           Translate DNA into amino acid sequences (requires sequences to be
                        DNA, in frame, and length multiple of 3)

Summaries:
  No sequences are printed when these options are used

  --num                 Print number of sequences
  --len                 Print summary of sequence lengths
  --com                 Print overall sequence composition
  --seqcom              Print composition for each individual sequence. Output is one line
                        per residue-type per sequence: seqname, residue-type, freq, count,
                        seqlength
  --ignoregaps          When reporting composition: do not count gap symbols
  --nam                 Print names of sequences
  --div                 (For alignments) Print nucleotide diversity (=average pairwise
                        sequence difference): mean, std, min, max
  --sit                 (For alignments) Print site summary: how many columns are variable,
                        contain multiple residues, contain gaps, or contain IUPAC ambiguity
                        symbols. Also keeps track of overlaps between these categories.
```


