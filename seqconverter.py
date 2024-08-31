#!/usr/bin/env python3

# Anders Gorm Pedersen, agpe@dtu.dk, 2012-2022
# Section for Bioinformatics, DTU Health Technology, Technical University of Denmark
# Converts between different sequence file formats. Performs various manipulations on sequences

# Has grown organically according to what was needed in my own projects.
# Badly needs refactoring and re-thinking of what options should be available
# Nevertheless does do lots of useful stuff with sequence files

##########################################################################################
##########################################################################################

import sys
import math
import os.path
import sequencelib as sq
import argparse
import re
from pathlib import Path

################################################################################################

def main():
    parser = build_parser()
    args = parser.parse_args()
    try:
        args = check_commandline(args)
        seqs,args = read_seqs(args)
        check_args_alignment_issues(args)

        if args.dupseqfilter:
            seqs.removedupseqs()

        if args.keepvar:
            seqs = positionfilter(seqs, args)

        seqs = change_seqs(seqs, args)

        if args.multifile:
            write_partitions(seqs, args)
        else:
            if any([args.s_nam, args.s_len, args.s_num, args.s_div, args.s_com, args.s_comseq, args.s_sit]):
                print_summary(seqs, args)
            else:
                print_seqs(seqs, args)

                if args.charset:
                    print("")
                    print(seqs.charsetblock())

                if args.mbpartblock:
                    print("")
                    print(seqs.mbpartblock())

    except sq.SeqError as exc:
        if args.debug:
            import traceback
            traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write("Error: {}\n".format(exc.errormessage))

################################################################################################

def build_parser():
    parser = argparse.ArgumentParser()

    #########################################################################################

    fileg = parser.add_argument_group("Input/Output")

    fileg.add_argument("-i", action="append", dest="filelist", metavar='SEQFILE',
                       default=None, type=Path,
                       help="One or more sequence files (repeat -i SEQFILE option for " +
                        "each input file). If -i SEQFILE is not given: take input from stdin "
                        +"(typically from a UNIX pipe).")

    fileg.add_argument("--informat", action='store', metavar="FORMAT", default="auto",
                      choices=["auto", "fasta", "nexus", "phylip", "clustal", "stockholm",
                               "genbank", "tab", "raw", "how"],
                      help="Input format: %(choices)s [default: %(default)s]")

    fileg.add_argument("--outformat", action='store', metavar="FORMAT", default="fasta",
                      choices=["fasta", "nexus", "phylip", "clustal", "tab", "raw", "how"],
                      help="Output format:  %(choices)s [default: %(default)s]")

    fileg.add_argument("--width", action="store", type=int, dest="width", metavar="WIDTH", default=60,
                        help="Print sequences with WIDTH characters per line [default: %(default)s] " 
                        + "Use the special value -1 (--width -1) to print each sequence in its "
                        + "entirety on a single line, regardless of its length.")

    #########################################################################################

    subsetg = parser.add_argument_group("Selecting subset of sequences")

    subsetg.add_argument("--sampleseq", action="store", type=int, dest="samplesize", metavar="N",
                        help="Randomly sample N sequences from sequence set")

    subsetg.add_argument("--keepreg", action="store", metavar='"REGEXP"',
                          help="Select sequences where substring of name matches regular expression")

    subsetg.add_argument("--remreg", action="store", metavar='"REGEXP"',
                          help="Discard sequences where substring of name matches regular expression")

    subsetg.add_argument("--keepname", action="store", dest="namefile", metavar="NAMEFILE",
                        help="Select sequences listed in NAMEFILE")

    subsetg.add_argument("--remname", action="store", dest="remfile", metavar="NAMEFILE",
                        help="Discard sequences listed in NAMEFILE")

    subsetg.add_argument("--keepvar", nargs='+', type=str, metavar="VARIANT",
                          help="Select sequences containing specific variants, i.e., specific "
                              + "residues on specific positions. "
                              + "Syntax for specifying VARIANT is: <POS><RESIDUE> (e.g., 484K). "
                              + "Multiple variants can be specifyed simultaneously separated "
                              + "by blanks. Example: --keepvar 484K 501Y")

    subsetg.add_argument("--remdupseq", action="store_true", dest="dupseqfilter",
                          help="Remove duplicate sequences (keeping one of each, randomly selected).")

    subsetg.add_argument("--remdupname", action="store_true", dest="dupnamefilter",
                          help="Remove sequences with duplicate names (keeping one of each, randomly selected). If this option is not set (default): stop execution on duplicate names.")

    subsetg.add_argument("--remendgapseqs", action="store", type=int, metavar="MIN",
                        help="Discard sequences with endgaps >= MIN positions. "
                            + "Endgaps are defined as contiguous block of gap symbols at either end of sequence.")

    #########################################################################################

    seqpartg = parser.add_argument_group("Selecting subset of positions in sequences")

    seqpartg.add_argument("--samplecols", action="store", type=int, metavar="N",
                        help="Randomly sample N columns from alignment")

    seqpartg.add_argument("--keepcols",  nargs='+', type=str, metavar="INDEX_OR_RANGE",
                          help="Keep alignment columns indicated by one or more INDEX_OR_RANGE values. "
                             + "INDEX_OR_RANGE values are either a single position (e.g., 15) "
                             + "or a range (e.g., 20-37). Multiple values shold be separated by "
                             + "blanks. Example: --keepcols 10 15 22-40 57")

    seqpartg.add_argument("--remcols", nargs='+', type=str, metavar="INDEX_OR_RANGE",
                          help="Remove alignment columns indicated by one or more INDEX_OR_RANGE values. "
                             + "INDEX_OR_RANGE values are either a single position (e.g., 15) "
                             + "or a range (e.g., 20-37). Multiple values shold be separated by "
                             + "blanks. Example: --remcols 10 15 22-40 57")

    seqpartg.add_argument("--remgapcols", nargs='?', const=0.0, default=None,
                          metavar="FRAC", action="store", type=float,
                          help="Remove columns that contain any gaps. "
                             + "If FRAC (number between 0-1) given: Remove columns "
                             + "where the fraction of gaps >= FRAC.")

    seqpartg.add_argument("--remambigcols", nargs='?', const=0.0, default=None,
                          metavar="FRAC", action="store", type=float,
                          help="Remove columns where one or more residues are ambiguity symbols "
                             + "(e.g., N for nucleotides). "
                             + "If FRAC (number between 0-1) given: Remove columns "
                             + "where the fraction of ambiguity symbols >= FRAC.")

    seqpartg.add_argument("--remendgapcols", nargs='?', const=0.0, default=None,
                          metavar="FRAC", action="store", type=float, dest="endgapfrac",
                          help="Remove columns where one or more sequences have endgaps. "
                             + "If FRAC (number between 0-1) given: Remove columns "
                             + "where the fraction of sequences having endgaps is >= FRAC. "
                             + "Endgaps are defined as contiguous block of gap symbols at either end of sequence")

    seqpartg.add_argument("--remconscols", action="store_true", dest="remconscols",
                          help="Remove conserved columns from alignment")

    seqpartg.add_argument("--remhmminsertcols", action="store_true", dest="remhmminsertcols",
                          help=argparse.SUPPRESS)  # Secret option
                          # help="When reading Stockholm format file from HMMer's hmmalign: remove columns corresponding to insert states")

    seqpartg.add_argument("--windows", action="store", type=int, dest="wsize", metavar="WSIZE",
                          help="For each sequence in input: extract all overlapping sequence windows "
                             + "of size WSIZE")

    seqpartg.add_argument("--degap", action="store_true", dest="degap",
                          help="Remove all gap characters from sequences")

    #########################################################################################

    renameg = parser.add_argument_group("Renaming sequences")

    renameg.add_argument("--rename", nargs=2, metavar=('OLD', 'NEW'),
                          help="Rename single sequence from OLD to NEW")

    renameg.add_argument("--renamenum", action="store", dest="renamenum", metavar="BASENAME",
                          help="Rename all sequences to this form: BASENAME_001, ...")

    renameg.add_argument("--renamereg",  nargs=2, metavar=('"OLD_REGEX"', '"NEW_STRING"'),
                          help="Rename sequences: Replace occurrences of regular expression "
                             + "OLD_REGEX with NEW_STRING")

    renameg.add_argument("--saverename", action="store", dest="savenamefile", metavar="NAMEFILE",
                          help="Save renaming information in NAMEFILE for later use")

    renameg.add_argument("--renamefile", action="store", metavar="NAMEFILE",
                          help="Replace names in sequence file using OLDNAME NEWNAME pairs "
                             + "in NAMEFILE. Not all names need to be listed. "
                             + "Note: can be used to restore names saved "
                             + "with --saverename during previous renaming.")

    renameg.add_argument("--gbname", action="store", dest="gbname", metavar="FIELD1[,FIELD2,FIELD3,...]",
                        help="For Genbank input: construct sequence names from the list of named fields, in the specified order")

    #########################################################################################

    multifileg = parser.add_argument_group("Combining multiple sequence files")

    multifileg.add_argument("--paste", action="store_true", dest="paste",
                          help="Concatenate identically named sequences from separate input files. " +
                          "Sequences are pasted end to end in the same order as the order of the input files. " +
                          "All input files must contain same number of sequences, and sequences " +
                          "in different files must have same name. (Order of sequences in individual " +
                          "file is not important)." +
                          "To see partitions choose nexus output, or output to multiple partition files.")

    multifileg.add_argument("--overlap", nargs='?', const=0, type=int, default=None, metavar="MIN",
                          help="Similar to --paste, but for input alignments that overlap partly at their " +
                          "ends. " +
                          "End-overlaps are discovered automatically and partition boundaries are then set " +
                          "such that each partition is covered by a unique set of genes. " +
                          "To see partitions choose nexus output, or output to multiple partition files. " +
                          "MIN: (optional, integer) minimum number of overlapping residues required for " +
                          "merging input alignments " +
                          "(default: set automatically based on seq lengths)")

    multifileg.add_argument("--multifile", action="store_true", dest="multifile",
                        help="Outputs to multiple files (one per partition) instead of stdout. " +
                          "Partitions are generated automatically based on other options.")

    multifileg.add_argument("--charset", action="store_true", dest="charset",
                        help="Appends Nexus form charset block listing partitions in data (forces output in Nexus format). " +
                          "Charsets and partitions are generated automatically based on other options.")

    multifileg.add_argument("--mb", action="store_true", dest="mbpartblock",
                        help="Appends MrBayes block with commands for running partitioned analysis (forces output in Nexus format). " +
                          "Charsets and partitions are generated automatically based on other options.")

    #########################################################################################

    dnag = parser.add_argument_group("DNA manipulations")

    dnag.add_argument("--revcomp", action="store_true", dest="revcomp",
                      help="Return reverse complement of sequence(s). Requires sequences to be DNA.")

    dnag.add_argument("--translate", type=int, default=None, metavar="READING_FRAME",
                      choices = [1,2,3], dest="reading_frame",
                      help="Translate input DNA sequences into amino acid sequences. "
                         + "READING_FRAME: either 1, 2, or 3, where 1 means start translation "
                         + "from first nucleotide in sequences. Translation includes as "
                         + "many full-length codons as possible, given READING_FRAME.")


    #########################################################################################

    summaryg = parser.add_argument_group("Summaries", description="No sequences are printed when these options are used")

    summaryg.add_argument("--nam", action="store_true", dest="s_nam",
                      help="Print names of sequences")

    summaryg.add_argument("--num", action="store_true", dest="s_num",
                      help="Print number of sequences")

    summaryg.add_argument("--len", action="store_true", dest="s_len",
                      help="Print summary of sequence lengths")

    summaryg.add_argument("--sit", action="store_true", dest="s_sit",
                      help="""(For alignments) Print site summary: how many columns are variable, contain multiple residues,
                              contain gaps, or contain IUPAC ambiguity symbols. Also keeps track of overlaps between
                              these categories, and the number of unique site patterns (columns)""")

    summaryg.add_argument("--com", action="store_true", dest="s_com",
                      help="Print overall sequence composition")

    summaryg.add_argument("--comseq", action="store_true", dest="s_comseq",
                      help="Print composition for each individual sequence. Output is one line per residue-type per sequence: seqname, residue-type, freq, count, seqlength")

    summaryg.add_argument("--div", action="store_true", dest="s_div",
                      help="(For alignments) Print nucleotide diversity (=average pairwise sequence difference):"
                      + " mean, std, min, max")

    summaryg.add_argument("--ignoregaps", action="store_true", dest="s_ignoregaps",
                      help="When computing composition or diversity: do not count gap symbols")

    #########################################################################################

    debug = parser.add_argument_group("Debugging")

    parser.add_argument("--debug", action="store_true", dest="debug",
                      help="Print longer error messages")

    #########################################################################################

    return parser

################################################################################################

def check_commandline(args):

    if args.filelist is None:
        args.filelist = ["-"]

    if args.informat == "auto":
        args.informat = "autodetect"         # Long name required by sq, short better for user interface...

    # Set flags indicating whether input and/or output is aligned, and whether file should be read as alignment
    # Note: args that require manipulation of columns, and --paste also implies aligned sequences
    # (reading as alignment results in check of equal sequence lengths)
    # Note 2: This is somewhat messy. Should really be handled by argparse mutually exclusive arguments
    alignin = alignout = args.aligned = False
    if any([args.informat=="nexus", args.informat=="phylip", args.informat=="clustal", args.informat=="stockholm"]):
        alignin = True
    if any([args.outformat=="nexus", args.outformat=="phylip", args.outformat=="clustal"]):
        alignout = True
    if any([alignin, alignout, args.remcols, args.remambigcols, args.remgapcols,
            args.remconscols, args.keepcols, args.paste, args.overlap]):
        args.aligned = True

    # If option --gbname is set: force input format to "genbank"
    if args.gbname:
        args.informat = "genbank"

    # Perform sanity checks on combination of args
    # (Note: if autodetection of filetype is requested, then many checks are impossible - defer error checking to sq)

    # Sanity check #1: option --degap cannot be used with any option that requires aligned sequences
    if args.degap and any([args.keepcols, args.remambigcols, args.remgapcols, args.remconscols,
                           args.mbpartblock, args.charset, args.multifile,
                           args.overlap, args.s_div, args.s_sit, args.paste]):
        raise sq.SeqError("Option --degap cannot be used with any option that requires aligned sequences".format(args.outformat))

    # Sanity check #2: option --degap cannot be used if output is an alignment
    if args.degap and alignout:
        raise sq.SeqError("Removal of all gap characters (--degap) not possible for output in {}".format(args.outformat))

    # Sanity check #3: option --keepcols cannot be used together with args --remambigcols, --remgapcols, or --remconscols
    if args.keepcols and args.remambigcols:
        raise sq.SeqError("Options --keepcols and --remambigcols cannot be used simultaneously")
    if args.keepcols and args.remgapcols:
        raise sq.SeqError("Options --keepcols and --remgapcols cannot be used simultaneously")
    if args.keepcols and args.remconscols:
        raise sq.SeqError("Options --keepcols and --remconscols cannot be used simultaneously")

    # Sanity check #5: option --saverename requires option --renamenum or --renamereg
    if args.savenamefile and not (args.renamenum or args.renamereg):
        raise sq.SeqError("Option --saverename requires option --renamenum or --renamereg")

    # Sanity check #6: args --renamenum and --renamefile are incompatible
    if args.renamenum and  args.renamefile:
        raise sq.SeqError("Option --renamenum and --renamefile cannot be used simultaneously")

    # Sanity check #7: args --renamereg and --renamefile are incompatible
    if args.renamereg and  args.renamefile:
        raise sq.SeqError("Option --renamereg and --renamefile cannot be used simultaneously")

    # Sanity check #8: args --renamereg and --paste are incompatible
    if args.renamereg and  args.paste:
        raise sq.SeqError("Option --renamereg and --paste cannot be used simultaneously: first rename, then paste resulting files")

    # Sanity check #9c: option --overlap is not meaningful unless several input files are provided
    if args.overlap and len(args.filelist) == 1:
        raise sq.SeqError("Option --overlap requires multiple input files")

    return (args)

################################################################################################

def read_seqs(args):

    for filename in args.filelist:
        if filename != "-" and not os.path.isfile(filename):
            raise sq.SeqError("File %s not found." % filename)

    # Read sequences from all files
    seqlist = []          # List of either Seq_set or Seq_alignment objects
    for filename in args.filelist:
        if args.gbname:
            seqfile = sq.Genbankfilehandle(filename, namefromfields=args.gbname)
        else:
            seqfile = sq.Seqfile(filename, args.informat)

        # If filetype was autodetected: check if this should be read as an alignment:
        if args.informat == "autodetect" and isinstance(seqfile, sq.Alignfile_reader):
            args.aligned = True

        if args.aligned:
            seqlist.append(seqfile.read_alignment(args.dupnamefilter))
        else:
            seqlist.append(seqfile.read_seqs(args.dupnamefilter))

    # Post processing to check if this is alignment, despite other flags not indicating it
    # This can happen if alignment is in e.g. fasta format, and none of the alignment status-triggering flags are used
    # Assume sequences are aligned if all sequences have same length:
    if not args.aligned:
        baselength = len(seqlist[0][0])      # Length of first sequence in first sequence object
        change_to_alignment = True
        for seqset in seqlist:
            for seq in seqset:
                if len(seq) != baselength:
                    change_to_alignment = False

        if change_to_alignment:
            newseqlist = []
            for seqset in seqlist:
                newseqset = sq.Seq_alignment()
                newseqset.addseqset(seqset)
                newseqlist.append(newseqset)
            seqlist = newseqlist
            args.aligned = True

    # If automatic overlap detection has been requested: find overlaps and concatenate or merge sub-alignments
    if args.overlap:

        # Create consensus sequence from each seq in seqset
        consensus_list = []
        alignment_dict = {}
        for alignment in seqlist:
            consensus_list.append( alignment.consensus() )
            alignment_dict[alignment.name] = alignment

        # Find overlap between consensus sequences, get list of contigs
        ra = sq.Read_assembler( consensus_list )
        contiglist = ra.assemble(minoverlap=args.overlap)

        # Build concatenated alignment from non-overlapping sub-alignments
        # Order and coordinates for each subalignment is based on info about regions in each contig from contiglist
        seqs = None
        for contig in contiglist:
            regionlist = contig.regions()
            for region in regionlist:
                sub_name, sub_start, sub_stop = region.get_align_tuple()     # Format is: (alignment_name, subseq_start, subseq_stop)
                subalignment = alignment_dict[sub_name].subseq(sub_start, sub_stop, rename=False)
                new_name = ".".join(region.name_list)       # Set name of partition to collection of all seqs in this overlap (possibly only one)
                subalignment.name = new_name
                if seqs is None:
                    seqs = sq.Seq_alignment(new_name)
                    seqs.addseqset(subalignment)
                else:
                    seqs = seqs.appendalignment(subalignment)

    # For other options than overlap: Collect sequences in one sequence collection object
    else:
        seqs = seqlist[0]
        for seqset in seqlist[1:]:
            if args.paste:
                seqs = seqs.appendalignment(seqset)
            else:
                seqs.addseqset(seqset, args.dupnamefilter)

    return seqs, args

################################################################################################

def check_args_alignment_issues(args):
    # After read_seqs we know whether sequences correspond to alignment.
    # Check if any requested option clashes with alignment status
    if not args.aligned:
        bad_options_list = []
        if args.keepcols:
            bad_options_list.append("--keepcols")
        if args.remendgapseqs:
            bad_options_list.append("--remendgapseqs")
        if args.remambigcols:
            bad_options_list.append("--remambigcols")
        if args.remgapcols:
            bad_options_list.append("--remgapcols")
        if args.endgapfrac:
            bad_options_list.append("--remendgapcols")
        if args.remconscols:
            bad_options_list.append("--remconscols")
        if args.mbpartblock:
            bad_options_list.append("--mbpartblock")
        if args.charset:
            bad_options_list.append("--charset")
        if args.multifile:
            bad_options_list.append("--multifile")
        if args.overlap:
            bad_options_list.append("--overlap")
        if args.s_div:
            bad_options_list.append("--div")
        if args.s_sit:
            bad_options_list.append("--sit")
        if args.paste:
            bad_options_list.append("--paste")
        if bad_options_list:
            bad_option_string = ", ".join(bad_options_list)
            raise sq.SeqError("The option(s) below require aligned sequences:\n\t{}".format(bad_option_string))

################################################################################################

def parse_indexlist(seqs, indexlist):
    """Parse INDEX-LIST string, create list of all individual positions implied by syntax
    Note: input assumed to be 1-indexed, output is 0-indexed ("slicesyntax")
    Example: input: 1-5,10,15-18, output: [0,1,2,3,4,9,14,15,16,17]"""

    poslist = []
    for item in indexlist:
       if "-" in item:
           (start,stop) = item.split("-")
           start = int(start) - 1       # Slice syntax: Numbering starts at 0
           stop = int(stop)             # Slice syntax: stop is one more than last index (=> do not subtract 1 from non-slice syntax)
           poslist.extend(list(range(start,stop)))
       else:
           poslist.append(int(item) - 1)

    return poslist

################################################################################################

def positionfilter(seqs, args):

    # Extract positions and residues from input of the form: ["484K", "501Y", "987W"].
    # Save as list of (pos, residue) tuples
    varlist = []
    for item in args.keepvar:
        residue = item[-1]
        pos = int(item[:-1]) - 1     # Convert to python array numbering from natural numbering
        varlist.append( (pos, residue) )

    # Iterate over sequences, and select those that contain all requested variants
    selseqs = seqs.__class__()       # Create new object of same class as seqs (don't know which subclass they are)
    for seq in seqs:
        select = True
        for (pos,res) in varlist:
            if seq[pos] != res:
                select = False
        if select:
            selseqs.addseq(seq)

    return selseqs

################################################################################################

def change_seqs(seqs, args):

    # Extract random subsample of sequences
    if args.samplesize:
        seqs = seqs.subsample(args.samplesize)

    # Extract sequences whose names match regexp in "keepreg"
    if args.keepreg:
        newseqs = sq.Seq_set()
        regexp = args.keepreg
        for seq in seqs:
            if re.search(regexp, seq.name):
                newseqs.addseq(seq)
        seqs = newseqs

    # Discard sequences whose names match regexp in "remreg"
    if args.remreg:
        newseqs = sq.Seq_set()
        regexp = args.remreg
        for seq in seqs:
            if not re.search(regexp, seq.name):
                newseqs.addseq(seq)
        seqs = newseqs

    # Extract sequences named in options.namefile
    if args.namefile:
        namelist = []
        namef = open(args.namefile)
        for line in namef:
            namelist.append(line.strip())
        namef.close()
        seqs = seqs.subset(namelist)

    # Remove sequences named in args.remfile
    if args.remfile:
        namelist = []
        namef = open(args.remfile)
        for line in namef:
            namelist.append(line.strip())
        namef.close()
        seqs.remseqs(namelist)

    # Remove sequences with endgaps longer than args.remendgapseqs positions
    if args.remendgapseqs is not None:
        seqs.remendgapseqs(args.remendgapseqs)

    # Removal of gaps from individual sequences
    if args.degap:
        # If sequences have been read as alignment: first convert to Seq_set:
        if args.aligned:
            args.aligned = False
            newseqs = sq.Seq_set()
            for seq in seqs:
                newseqs.addseq(seq)
            seqs = newseqs
        seqs.remgaps()

    # Random sampling of columns
    if args.samplecols:
        seqs.samplecols(args.samplecols)

    # Removal of listed columns
    if args.remcols:
        remlist = parse_indexlist(seqs, args.remcols)
        seqs.remcols(remlist)

    # Extraction of listed columns
    if args.keepcols:
        keeplist = parse_indexlist(seqs, args.keepcols)
        seqs.indexfilter(keeplist)

    # Removal of ambiguity columns
    if args.remambigcols is not None:
        if args.remambigcols == 0.0:
            seqs.remambigcol()
        else:
            seqs.remfracambigcol(args.remambigcols)

    # Removal of gappy columns
    if args.remgapcols is not None:
        if args.remgapcols == 0.0:
            seqs.remgapcol()
        else:
            seqs.remfracgapcol(args.remgapcols)

    # Removal of columns with > FRAC fraction endgaps
    if args.endgapfrac is not None:
        seqs.remendgapcol(args.endgapfrac)

    # Removal of conserved columns
    if args.remconscols:
        seqs.remconscol()

    # Removal of insert state columns (in output from HMMer's hmmalign)
    if args.remhmminsertcols:
        seqs.rem_hmmalign_insertcol()

    # Extraction of sequence windows
    if args.wsize:
        newseqs = sq.Seq_set()
        for seq in seqs:
            for seqwin in seq.windows(args.wsize, rename=True):
                newseqs.addseq(seqwin)
        seqs = newseqs

    # Renaming
    if args.rename:
        (oldname, newname) = args.rename.split(",")
        seqs.changeseqname(oldname, newname)
    elif args.renamenum:
        seqs.rename_numbered(args.renamenum, args.savenamefile)
    elif args.renamereg:
        old_regex, new_string = args.renamereg
        seqs.rename_regexp(old_regex, new_string, args.savenamefile)

    # Restore names
    if args.renamefile:
        seqs.transname(args.renamefile)


    # Translation
    if args.reading_frame:
        seqs = seqs.translate(args.reading_frame)

    # Reverse complement
    if args.revcomp:
        seqs = seqs.revcomp()

    # Return altered sequences
    return seqs

################################################################################################

def print_summary(seqs, args):

    if args.s_nam:
        for name in sorted(seqs.seqnamelist):
            print(name)

    if args.s_num:
        print("Number of sequences: {:10d}\n".format(len(seqs)))

    if args.s_len:
        if not args.aligned:
            seqlenlist = []
            lensum = 0
            for seq in seqs:
                seqlen = len(seq)
                seqlenlist.append(seqlen)
                lensum += seqlen
            print("Minimum length:   {:13d}".format(min(seqlenlist)))
            print("Maximum length:   {:13d}".format(max(seqlenlist)))
            print("Average length:   {:13.2f}\n\n".format(lensum/len(seqs)))
        else:
            print("Alignment length:    {:10d}\n".format( seqs.alignlen() ))

    if args.s_div:
        avg, std, minpi, maxpi = seqs.sequence_diversity(ignoregaps=args.s_ignoregaps)
        print("Nucleotide diversity pi (average pairwise sequence distance)")
        print(f"    Mean:               {avg:.5f}")
        print(f"    Standard dev:       {std:.5f}")
        print(f"    Min:                {minpi:.5f}")
        print(f"    Max:                {maxpi:.5f}\n")

    if args.s_sit:
        numvar = 0
        nummulti = 0
        numgap = 0
        numambig = 0
        numvar_only = 0
        nummulti_gap = 0
        nummulti_ambig = 0
        nummulti_gap_ambig = 0
        numgap_ambig = 0

        ambigsymbols_set = set(seqs.ambigsymbols)
        sitepatterns = set()
        for col in seqs.columns():
            sitepatterns.add(tuple(col))
            columnset = set(col)
            if len(columnset) != 1:
                numvar += 1
                if "-" not in columnset and not (ambigsymbols_set & columnset):
                    numvar_only += 1
                if len(columnset - {"-", *ambigsymbols_set}) > 1:
                    nummulti += 1
                    if "-" in columnset:
                        nummulti_gap += 1
                        if ambigsymbols_set & columnset:
                            nummulti_gap_ambig += 1
                    if ambigsymbols_set & columnset:
                        nummulti_ambig += 1
            if "-" in columnset:
                numgap += 1
                if ambigsymbols_set & columnset:
                    numgap_ambig += 1
            if ambigsymbols_set & columnset:
                numambig += 1
        numconst = seqs.alignlen() - numvar

        align_length = seqs.alignlen()

        def percentage(num):
            return (num / align_length) * 100

        print("Site summary")
        print(f"    Constant sites:                      {numconst:>6,d}  {percentage(numconst):>4.1f}%")
        print(f"    Variable sites:                      {numvar:>6,d}  {percentage(numvar):>4.1f}%")
        print(f"    Multi-residue sites:                 {nummulti:>6,d}  {percentage(nummulti):>4.1f}%")
        print(f"    Gappy sites:                         {numgap:>6,d}  {percentage(numgap):>4.1f}%")
        print(f"    Ambiguous sites:                     {numambig:>6,d}  {percentage(numambig):>4.1f}%")
        print("")
        print(f"    Multi-residue-only sites:            {numvar_only:>6,d}  {percentage(numvar_only):>4.1f}%")
        print(f"    Multi-residue gappy sites:           {nummulti_gap:>6,d}  {percentage(nummulti_gap):>4.1f}%")
        print(f"    Multi-residue ambiguous sites:       {nummulti_ambig:>6,d}  {percentage(nummulti_ambig):>4.1f}%")
        print(f"    Multi-residue gappy-ambiguous sites: {nummulti_gap_ambig:>6,d}  {percentage(nummulti_gap_ambig):>4.1f}%")
        print(f"    Gappy-ambiguous sites:               {numgap_ambig:>6,d}  {percentage(numgap_ambig):>4.1f}%")
        print("")
        print(f"    Number of unique site patterns:      {len(sitepatterns):>6,d}\n")

    if args.s_com:
        compositiondict = seqs.composition(ignoregaps=args.s_ignoregaps)

        print("Composition:")
        print("    {:^7s}{:>8s}{:>10s}".format("Symbol", "Freq", "Count"))
        tot_count = 0
        for res,(count,freq) in compositiondict.items():
            tot_count += count

        nonambig = set(compositiondict.keys()) - seqs.ambigsymbols - set("-")
        ambig = set(compositiondict.keys()) - nonambig
        nonambig_sorted = sorted(list(nonambig))
        ambig_sorted = sorted(list(ambig))
        for res in nonambig_sorted:
            count,freq = compositiondict[res]
            print("    {:^7s}{:>8.3f}{:>10d} / {}".format(res, freq, count, tot_count))
        for res in ambig_sorted:
            count,freq = compositiondict[res]
            print("    {:^7s}{:>8.3f}{:>10d} / {}".format(res, freq, count, tot_count))
        print()

    if args.s_comseq:
        longestname = max(seqs.getnames(), key=len)
        maxlenname = max(len(longestname), len("name"))
        print("{name:<{nl}s}{res:^8s}{f:>8s}{n:>10s}{l:>10s}".format(name="name", nl=maxlenname+2,
                                                                    res="residue", f="freq",
                                                                    n="count", l="length"))
        for seq in seqs:
            length = len(seq)
            name = seq.name
            compositiondict = seq.composition(ignoregaps=args.s_ignoregaps)
            nonambig = set(compositiondict.keys()) - seq.ambigsymbols - set("-")
            ambig = set(compositiondict.keys()) - nonambig
            nonambig_sorted = sorted(list(nonambig))
            ambig_sorted = sorted(list(ambig))
            for res in nonambig_sorted:
                count,freq = compositiondict[res]
                print("{name:<{nl}s}{res:^8s}{f:>8.4f}{n:>10d}{l:>10d}".format(name=name, nl=maxlenname+2,
                                                                            res=res, f=freq, n=count, l=length))
            for res in ambig_sorted:
                count,freq = compositiondict[res]
                print("{name:<{nl}s}{res:^8s}{f:>8.4f}{n:>10d}{l:>10d}".format(name=name, nl=maxlenname+2,
                                                                            res=res, f=freq, n=count, l=length))

################################################################################################

def print_seqs(seqs, args, filehandle=sys.stdout):

    # NOTE: default is to write to stdout (print to terminal) but other filehandle can be set

    # Some options implies Nexus format output
    if args.charset or args.mbpartblock:
        args.outformat = "nexus"

    # Print sequences in requested format
    if args.outformat == "raw":
        filehandle.write(seqs.raw() + '\n')
    elif args.outformat == "tab":
        filehandle.write(seqs.tab() + '\n')
    elif args.outformat == "fasta":
        filehandle.write(seqs.fasta(width = args.width) + '\n')
    elif args.outformat == "how":
        filehandle.write(seqs.how(width = args.width) + '\n')
    elif args.outformat == "nexus":
        if any([args.paste, args.overlap, args.charset]) and (not args.multifile):
            parts = True
        else:
            parts = False
        filehandle.write(seqs.nexus(print_partitioned = parts, width = args.width) + '\n')
    elif args.outformat == "phylip":
        filehandle.write(seqs.phylip(width = args.width) + '\n')
    elif args.outformat == "clustal":
        filehandle.write(seqs.clustal(width = args.width) + '\n')

################################################################################################

def write_partitions(seqs, args):
    part_alignments = seqs.partitions_as_seqalignments()
    for alignment in part_alignments:
        outname = "partition_{alname}.{suffix}".format(alname=alignment.name, suffix=args.outformat)
        if os.path.isfile(outname):
            raise sq.SeqError("File with name {} already exists. Aborting to avoid overwriting.".format(outname))
        fh = open(outname, "w")
        print_seqs(alignment, args, filehandle=fh)
        fh.close()

################################################################################################

if __name__ == "__main__":
    main()
