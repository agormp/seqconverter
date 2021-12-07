# seqconverter: command line program for reading, writing, and manipulating sequence files

![](https://img.shields.io/badge/version-2.1.1-blue)

The command-line program `seqconverter` can read and write text files containing aligned or unaligned DNA or protein sequences. The program understands most standard and some non-standard formats (fasta, Nexus, Phylip, Clustal, tab, raw, Genbank, How). The program can perform various manipulations on the sequences.

## Availability

The seqconverter source code is available on GitHub: https://github.com/agormp/seqconverter and can be installed from PyPI: https://pypi.org/project/seqconverter/

## Installation

```
python3 -m pip install seqconverter
```

## Highlights

* Can be used to convert between sequence file formats but also does other things
* Read and write **aligned** sequences in the following formats:
	* fasta
	* Nexus
	* Phylip
	* Clustal
	* tab
	* raw
* Read and write **unaligned** sequences in the following formats:
	* fasta
	* tab
	* raw
	* Genbank
	* How
* Writes to stdout, so output can be used in pipes or redirected to file
* Extract subsequence (specified columns) from alignment
* Extract all overlapping windows of specified size
* Extract named sequences from set of sequences
* Randomly sample from set of sequence
* Remove columns from alignment based on one of several criteria (all gaps, some gaps, more than fraction gaps, conserved, specified indices)
* Rename sequences automatically or using file with pairs of "oldname newname"
* Generate partitioned Nexus file with `charset` specification automatically from separate files containing identically named sequences (sequences are concatenated end to end in same order as files).
* More...
* Underlying library has been optimized for high speed and low memory consumption
* Really has too many options, but does useful stuff (and has been created based on what I needed for own projects)

## Usage examples

Get help:
```
seqconverter -h
```

----------------------------------------------------------------

Convert aligned sequences in fasta format to nexus (Note: output is written to the terminal so you need to use redirection to store in a file):
```
seqconverter -I fasta -O nexus myalignment.fasta > myalignment.nexus
```

----------------------------------------------------------------

Extract columns 50-150 (inclusive, with numbering starting at 1) from alignment in Clustal format, write output in fasta format to file (using redirection):
```
seqconverter -I clustal -O fasta --subseq 50,150 myalignment.aln > aligment_50_150.fasta
```

----------------------------------------------------------------

Extract all sequences containing a Lysine at position 484 and a Tyrosine at position 501 from set of amino acid sequences:
```
seqconverter -I clustal -O fasta --filterpos 484K,501Y myalignment.aln > voc.fasta
```

----------------------------------------------------------------

Remove columns where one or more residues are gaps from alignment:
```
seqconverter -I fasta -O fasta --remgapcols myalignment.fasta > gapfree.fasta
```

----------------------------------------------------------------

Concatenate identically named sequences from separate input files:
```
seqconverter -I fasta -O fasta --paste alignm1.fasta alignm2.fasta alignm3.fasta > concat.fasta
```

----------------------------------------------------------------

Concatenate identically named sequences from separate input files, creating partitioned Nexus file with `charset` specification. This can be used for phylogenetic analyses in BEAST or MrBayes where different genomic regions (e.g., genes) have different substitution models. Note: sequences in each file need to have identical names (e.g. name of species) and sequences in each file needs to be already aligned.
```
seqconverter -I fasta -O nexus --paste --charset gene1.fasta gene2.fasta gene3.fasta > partitioned.nexus
```


## Usage

```
usage: seqconverter [-h] [-I FORMAT] [-O FORMAT] [--nocomments] [--subseq START,STOP]
                    [--subseqrename] [--windows WSIZE] [--subsample N] [--subset NAMEFILE]
                    [--remseqs NAMEFILE] [--filterpos VARIANT[,VARIANT,...]] [--gff FILE]
                    [--gffsymbol ANNOT] [--degap] [--remcols INDEX LIST] [--remambigcols]
                    [--remgapcols] [--remallgapcols] [--remfracgapcols FRAC] [--remconscols]
                    [--paste] [--overlap] [--minoverlap N] [--multifile] [--charset]
                    [--mbpartblock] [--bestblock] [--filterdupseq] [--filterdupname]
                    [--gbname FIELD1[,FIELD2,FIELD3,...]] [--appendnumber] [--rename OLD,NEW]
                    [--renamenumber BASENAME] [--renameregexp "REGEXP"] [--regdupfix]
                    [--savenames FILE] [--restorenames FILE] [--revcomp] [--translate]
                    [--summary] [--names] [--debug]
                    SEQFILE [SEQFILE ...]

positional arguments:
  SEQFILE               One or more sequence files

optional arguments:
  -h, --help            show this help message and exit
  -I FORMAT             Input format: auto, fasta, tab, raw, genbank, how, clustal, phylip,
                        nexus [default: auto]
  -O FORMAT             Output format: fasta, tab, raw, how, clustal, phylip, nexus,
                        nexusgap, nexusmulti [default: fasta]
  --nocomments          Do not print comments in output (only print seqnames)
  --subseq START,STOP   Extract subsequence, positions START to STOP, from alignment
  --subseqrename        When extracting sub-sequences: add '_START_STOP' to seqnames
  --windows WSIZE       Extract all overlapping sequence windows of size WSIZE
  --subsample N         Randomly extract N sequences from sequence set
  --subset NAMEFILE     Extract sequences listed in NAMEFILE
  --remseqs NAMEFILE    Remove sequences listed in NAMEFILE
  --filterpos VARIANT[,VARIANT,...]
                        Extract sequences containing specific residues on specific positions.
                        Syntax is: <POS><RESIDUE>, possibly in a comma-separated list.
                        Example: 484K,501Y
  --gff FILE            Get annotation info for sequences from GFF-formatted FILE
  --gffsymbol ANNOT     Use the character ANNOT in the sequence annotation strings derived
                        from GFF file
  --degap               Remove all gap characters from sequences
  --remcols INDEX LIST  Remove listed columns from alignment. Columns can be indicated as
                        comma-separated list of indices, and as ranges. Example:
                        --remcols=10,15,22-40,57
  --remambigcols        Remove columns where one or more residues are ambiguity symbols
                        (e.g., N for nucleotides)
  --remgapcols          Remove columns where one or more residues are gaps
  --remallgapcols       Remove columns that are all-gaps
  --remfracgapcols FRAC
                        Remove columns that contain > FRAC fraction gaps
  --remconscols         Remove conserved columns from alignment
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
  --bestblock           Appends MrBayes block with commands for running BEST analysis to
                        output file (forces output in Nexus format). Charsets and partitions
                        are generated automatically so they correspond to input files (one
                        partition for each file).
  --filterdupseq        Remove duplicate sequences (keeping one of each); print names of
                        removed sequences on stderr.
  --filterdupname       Remove sequences with duplicate names (keeping one of each). If this
                        option is not set (default): stop execution on duplicate names.
  --gbname FIELD1[,FIELD2,FIELD3,...]
                        For Genbank input: construct sequence names from the list of named
                        fields, in the specified order
  --appendnumber        Append numbering at end of existing sequence names (SeqA_001,
                        SeqXYZ_002, ...
  --rename OLD,NEW      Rename sequence from OLD to NEW
  --renamenumber BASENAME
                        Rename sequences to this form: BASENAME_001, ...
  --renameregexp "REGEXP"
                        Rename sequences by deleting parts of names matching regular
                        expression in REGEXP
  --regdupfix           Fix duplicate names, created by regexp, by appending numbers to
                        duplicates (seqA, seqA_2, ...)
  --savenames FILE      Save renaming information in FILE for later use
  --restorenames FILE   Restore original names using info previously saved in FILE
  --revcomp             Return reverse complement of sequence(s).
  --translate           Translate DNA into amino acid sequences (requires sequences to be
                        DNA, in frame, and length multiple of 3)
  --summary             Print summary of data set (names, number, lengths, composition,
                        etc.). No sequences are output.
  --names               Print names of sequences in data set.
  --debug               Print longer error messages
Command me, Master>
```


