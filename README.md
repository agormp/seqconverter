# seqconverter

![](https://img.shields.io/badge/version-3.2.5-blue)
[![PyPI downloads](https://static.pepy.tech/personalized-badge/seqconverter?period=total&units=none&left_color=black&right_color=blue&left_text=PyPI%20downloads&service=github)](https://pepy.tech/project/seqconverter)
[![DOI](https://zenodo.org/badge/433106305.svg)](https://zenodo.org/doi/10.5281/zenodo.10411474)

The command-line program `seqconverter` can read and write text files containing aligned or unaligned DNA or protein sequences. The program understands most standard and some non-standard formats (fasta, Nexus, Phylip, Clustal, Stockholm, tab, raw, Genbank, How). The tool can be used to convert between sequence file formats, and is also able to perform various manipulations and analyses of sequences.

## Availability

The seqconverter source code is available on GitHub: https://github.com/agormp/seqconverter. The executable can be installed from PyPI: https://pypi.org/project/seqconverter/

## Version 3

Version 3 has recently been released, and contains a number of changes to the user-interface compared to version 2.x.x. For a full overview see notes in the latest release.

## Installation

```
python3 -m pip install seqconverter
```

Upgrading to latest version:

```
python3 -m pip install --upgrade seqconverter
```
## Citation

To cite seqconverter: use the link in the right sidebar under About --> Cite this repository.

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
* Options to select or remove columns from alignment based on one of several criteria: some gaps, more than fraction gaps, more than fration endgaps, conserved, specified indices, random sample of columns, ...
* Extract all overlapping windows of specified size
* Options to rename one or more sequences based on various criteria
* Options to concatenate identically named sequences from multiple sequence files (end-to-end or discarding automatically discovered overlaps)
* Options to automatically create Nexus charset commands based on merging multiple individual files (e.g., one charset/partition per gene).
* Can automatically write MrBayes block with template for commands to run partitioned analysis, also based on merging multiple separate sequence alignments.
* Can translate and find reverse complement for DNA sequences
* Options to obtain summary information about sequences and alignments: number of seqs, names, lengths, composition (overall or per sequence), nucleotide diversity (pi), site summary (how many columns are variable, contain multiple residues, contain gaps, or contain IUPAC ambiguity symbols, how many unique site patterns)
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

```
seqconverter --informat fasta --outformat nexus \
             --width 70 -i myalignment.fasta > myalignment.nexus
```

Note 1: output is written to the terminal so you need to use redirection to store in a file.
Note 2: input format will be automatically detected if not specified with --informat (this works well for standard file types)


----

### Select all sequences whose name match the regular expression "seq_1[0-9]+"

```
seqconverter --informat fasta --outformat fasta \
             --keepreg "seq_1[0-9]+" -i myseqs.fasta > subset.fasta
```
Note: default output format is fasta, so you do not need to specify --outformat fasta

----

### Discard all sequences whose name match the regular expression "seq_1[0-9]+":

```
seqconverter --informat fasta --outformat fasta \
             --remreg "seq_1[0-9]+" -i myseqs.fasta > subset.fasta
```

----

### Select random subset of 50 sequences from input file

```
seqconverter --informat fasta --outformat fasta \
             --sampleseq 50 -i myseqs.fasta > subset.fasta
```

----

### Select all sequence variants containing a Lysine at position 484 and a Tyrosine at position 501 

```
seqconverter --informat clustal --outformat fasta \
             --keepvar 484K 501Y -i myalignment.aln > voc.fasta
```

----

### Select columns 50-150 from ClustalW formatted alignment file, write output in fasta

```
seqconverter --informat clustal --outformat fasta \
             --keepcols 50-150 -i myalignment.aln > aligment_50_150.fasta
```

----

### Remove columns, where one or more residues are gaps, from alignment:

```
seqconverter --informat fasta --outformat fasta \
             --remgapcols -i myalignment.fasta > gapfree.fasta
```

----

### Remove columns, where >= 75% are gaps, from alignment:

```
seqconverter --informat fasta --outformat fasta \
             --remgapcols 0.75 -i myalignment.fasta > fewergaps.fasta
```

----

### Remove columns, where more than 75% have endgaps, from alignment:

This command will remove alignment columns if more than 75% of sequences have endgaps in that position. An endgap is defined as a contiguous gappy region at either the beginning or end of a sequence, and are often a result of missing data (the gaps then do not represent insertion or deletion events).

```
seqconverter --informat fasta --outformat fasta \
             --remendgapcols 0.75 -i myalignment.fasta > fewer_endgaps.fasta
```

----

### Concatenate identically named sequences from separate input files:

Sequences are pasted end to end in the same order as the order of the input files. All input files must contain the same number of sequences, and sequences in different files must have same name (for instance each file could contain an alignment of the sequences for a specific gene from a number of different species, and each sequence could then have the name of the species). The order of sequences in different files does not matter.

When used with the --charset (and possibly --mb) option this can be used to set up a partitioned analysis in MrBayes or BEAST (see below).

```
seqconverter --informat fasta --outformat fasta \
             --paste -i gene1.fasta -i gene2.fasta -i gene3.fasta > concat.fasta
```

----

### Concatenate sequences from multiple files, create partitioned Nexus file containing charset command

This command concatenates identically named sequences from separate input alignments, creating a partitioned Nexus file with `charset` specification. Start and stop indices for different charsets are automatically derived from lengths of sub-alignments. Charsets are named based on the names of included files.

This can be used for phylogenetic analyses in BEAST or MrBayes where different genomic regions (e.g., genes) have different substitution models. Note: sequences in each file need to have identical names (e.g. name of species).

```
seqconverter --outformat nexus --paste \
             --charset -i gene1.fasta -i gene2.fasta -i gene3.fasta > partitioned.nexus
```

----

### Concatenate sequences from multiple files, create partitioned Nexus file with commands to run MrBayes or BEAST analysis

This command does the same as the example above, and additionally adds a MrBayes block containing commands to run a partitioned analysis. The commands have sensible default values (e.g., setting DNA substution models to "nst=mixed" and unlinking most parameters across partitions). Optimally the commands should be tweaked according to the concrete data set. Importing the Nexus file in BEAUTI should result in setting most corresponding options for a BEAST run (but check, and remember to set priors etc.)

```
seqconverter --outformat nexus --paste \
             --charset --mb -i gene1.fasta -i gene2.fasta -i gene3.fasta > partitioned.nexus
```


## Usage

```
usage: seqconverter [-h] [-i SEQFILE] [--informat FORMAT] [--outformat FORMAT]
                    [--width WIDTH] [--sampleseq N] [--keepreg "REGEXP"]
                    [--remreg "REGEXP"] [--keepname NAMEFILE] [--remname NAMEFILE]
                    [--keepvar VARIANT [VARIANT ...]] [--remdupseq] [--remdupname]
                    [--remendgapseqs MIN] [--samplecols N]
                    [--keepcols INDEX_OR_RANGE [INDEX_OR_RANGE ...]]
                    [--remcols INDEX_OR_RANGE [INDEX_OR_RANGE ...]]
                    [--remgapcols [FRAC]] [--remambigcols [FRAC]]
                    [--remendgapcols [FRAC]] [--remconscols] [--windows WSIZE] [--degap]
                    [--rename OLD NEW] [--renamenum BASENAME]
                    [--renameregex "OLD_REGEX" "NEW_STRING"] [--saverename NAMEFILE]
                    [--renamefile NAMEFILE] [--gbname FIELD1[,FIELD2,FIELD3,...]]
                    [--paste] [--overlap [MIN]] [--multifile] [--charset] [--mb]
                    [--revcomp] [--translate READING_FRAME] [--nam] [--num] [--len]
                    [--sit] [--com] [--comseq] [--div] [--ignoregaps] [--debug]

options:
  -h, --help            show this help message and exit
  --debug               Print longer error messages

Input/Output:
  -i SEQFILE            One or more sequence files (repeat -i SEQFILE option for each
                        input file). If -i SEQFILE is not given: take input from stdin
                        (typically from a UNIX pipe).
  --informat FORMAT     Input format: auto, fasta, nexus, phylip, clustal, stockholm,
                        genbank, tab, raw, how [default: auto]
  --outformat FORMAT    Output format: fasta, nexus, phylip, clustal, tab, raw, how
                        [default: fasta]
  --width WIDTH         Print sequences with WIDTH characters per line [default: 60]

Select subset of sequences:
  --sampleseq N         Randomly sample N sequences from sequence set
  --keepreg "REGEXP"    Select sequences where substring of name matches regular
                        expression
  --remreg "REGEXP"     Discard sequences where substring of name matches regular
                        expression
  --keepname NAMEFILE   Select sequences listed in NAMEFILE
  --remname NAMEFILE    Discard sequences listed in NAMEFILE
  --keepvar VARIANT [VARIANT ...]
                        Select sequences containing specific variants, i.e., specific
                        residues on specific positions. Syntax for specifying VARIANT
                        is: <POS><RESIDUE> (e.g., 484K). Multiple variants can be
                        specifyed simultaneously separated by blanks. Example: --keepvar
                        484K 501Y
  --remdupseq           Remove duplicate sequences (keeping one of each, randomly
                        selected).
  --remdupname          Remove sequences with duplicate names (keeping one of each,
                        randomly selected). If this option is not set (default): stop
                        execution on duplicate names.
  --remendgapseqs MIN   Discard sequences with endgaps >= MIN positions. Endgaps are
                        defined as contiguous block of gap symbols at either end of
                        sequence.

Select subset of positions in sequences:
  --samplecols N        Randomly sample N columns from alignment
  --keepcols INDEX_OR_RANGE [INDEX_OR_RANGE ...]
                        Keep alignment columns indicated by one or more INDEX_OR_RANGE
                        values. INDEX_OR_RANGE values are either a single position
                        (e.g., 15) or a range (e.g., 20-37). Multiple values shold be
                        separated by blanks. Example: --keepcols 10 15 22-40 57
  --remcols INDEX_OR_RANGE [INDEX_OR_RANGE ...]
                        Remove alignment columns indicated by one or more INDEX_OR_RANGE
                        values. INDEX_OR_RANGE values are either a single position
                        (e.g., 15) or a range (e.g., 20-37). Multiple values shold be
                        separated by blanks. Example: --remcols 10 15 22-40 57
  --remgapcols [FRAC]   Remove columns that contain any gaps. If FRAC (number between
                        0-1) given: Remove columns where the fraction of gaps >= FRAC.
  --remambigcols [FRAC]
                        Remove columns where one or more residues are ambiguity symbols
                        (e.g., N for nucleotides). If FRAC (number between 0-1) given:
                        Remove columns where the fraction of ambiguity symbols >= FRAC.
  --remendgapcols [FRAC]
                        Remove columns where one or more sequences have endgaps. If FRAC
                        (number between 0-1) given: Remove columns where the fraction of
                        sequences having endgaps is >= FRAC. Endgaps are defined as
                        contiguous block of gap symbols at either end of sequence
  --remconscols         Remove conserved columns from alignment
  --windows WSIZE       For each sequence in input: extract all overlapping sequence
                        windows of size WSIZE
  --degap               Remove all gap characters from sequences

Renaming sequences:
  --rename OLD NEW      Rename single sequence from OLD to NEW
  --renamenum BASENAME  Rename all sequences to this form: BASENAME_001, ...
  --renameregex "OLD_REGEX" "NEW_STRING"
                        Rename sequences: Replace occurrences of regular expression
                        OLD_REGEX with NEW_STRING
  --saverename NAMEFILE
                        Save renaming information in NAMEFILE for later use
  --renamefile NAMEFILE
                        Replace names in sequence file using OLDNAME NEWNAME pairs in
                        NAMEFILE. Not all names need to be listed. Note: can be used to
                        restore names saved with --saverename during previous renaming.
  --gbname FIELD1[,FIELD2,FIELD3,...]
                        For Genbank input: construct sequence names from the list of
                        named fields, in the specified order

Combining multiple sequence files:
  --paste               Concatenate identically named sequences from separate input
                        files. Sequences are pasted end to end in the same order as the
                        order of the input files. All input files must contain same
                        number of sequences, and sequences in different files must have
                        same name. (Order of sequences in individual file is not
                        important).To see partitions choose nexus output, or output to
                        multiple partition files.
  --overlap [MIN]       Similar to --paste, but for input alignments that overlap partly
                        at their ends. End-overlaps are discovered automatically and
                        partition boundaries are then set such that each partition is
                        covered by a unique set of genes. To see partitions choose nexus
                        output, or output to multiple partition files. MIN: (optional,
                        integer) minimum number of overlapping residues required for
                        merging input alignments (default: set automatically based on
                        seq lengths)
  --multifile           Outputs to multiple files (one per partition) instead of stdout.
                        Partitions are generated automatically based on other options.
  --charset             Appends Nexus form charset block listing partitions in data
                        (forces output in Nexus format). Charsets and partitions are
                        generated automatically based on other options.
  --mb                  Appends MrBayes block with commands for running partitioned
                        analysis (forces output in Nexus format). Charsets and
                        partitions are generated automatically based on other options.

DNA manipulations:
  --revcomp             Return reverse complement of sequence(s). Requires sequences to
                        be DNA.
  --translate READING_FRAME
                        Translate input DNA sequences into amino acid sequences.
                        READING_FRAME: either 1, 2, or 3, where 1 means start
                        translation from first nucleotide in sequences. Translation
                        includes as many full-length codons as possible, given
                        READING_FRAME.

Summaries:
  No sequences are printed when these options are used

  --nam                 Print names of sequences
  --num                 Print number of sequences
  --len                 Print summary of sequence lengths
  --sit                 (For alignments) Print site summary: how many columns are
                        variable, contain multiple residues, contain gaps, or contain
                        IUPAC ambiguity symbols. Also keeps track of overlaps between
                        these categories, and the number of unique site patterns
                        (columns)
  --com                 Print overall sequence composition
  --comseq              Print composition for each individual sequence. Output is one
                        line per residue-type per sequence: seqname, residue-type, freq,
                        count, seqlength
  --div                 (For alignments) Print nucleotide diversity (=average pairwise
                        sequence difference): mean, std, min, max
  --ignoregaps          When computing composition or diversity: do not count gap
                        symbols
```


