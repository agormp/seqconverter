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
import os.path
import sequencelib as seqlib
import argparse

################################################################################################

def main():
    parser = build_parser()
    args = parser.parse_args()

    try:
        args = check_commandline(args)
        seqs,args = read_seqs(args)
        check_args_alignment_issues(args)

        if args.dupseqfilter:
            seqs = filterdupseqs(seqs)

        if args.filterpos:
            seqs = positionfilter(seqs, args)

        seqs = change_seqs(seqs, args)

        if args.multifile:
            write_partitions(seqs, args)
        else:
            if any([args.s_nam, args.s_len, args.s_num, args.s_div, args.s_com, args.s_seqcom, args.s_sit]):
                print_summary(seqs, args)
            else:
                print_seqs(seqs, args)

                if args.charset:
                    print("")
                    print(seqs.charsetblock())

                if args.mbpartblock:
                    print("")
                    print(seqs.mbpartblock())

    except seqlib.SeqError as exc:
        if args.debug:
            import traceback
            traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write("Error: {}\n".format(exc.errormessage))

################################################################################################

def build_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument("filelist", nargs="*", metavar='SEQFILE', default="-", help="One or more sequence files")

    #########################################################################################

    formatg = parser.add_argument_group("File formats")

    formatg.add_argument("-I", action='store', dest="informat", metavar="FORMAT", default="auto",
                      choices=["auto", "fasta", "nexus", "phylip", "clustal", "stockholm",
                               "genbank", "tab", "raw", "how"],
                      help="Input format: %(choices)s [default: %(default)s]")

    formatg.add_argument("-O", action='store', dest="outformat", metavar="FORMAT", default="fasta",
                      choices=["fasta", "nexus", "nexusgap", "phylip", "clustal", "tab", "raw", "how"],
                      help="Output format:  %(choices)s [default: %(default)s]")

    formatg.add_argument("--width", action="store", type=int, dest="width", metavar="WIDTH", default=60,
                        help="Print sequences with WIDTH characters per line [default: %(default)s]")

    #########################################################################################

    renameg = parser.add_argument_group("Renaming sequences")

    renameg.add_argument("--rename", action="store", dest="rename", metavar="OLD,NEW",
                          help="Rename single sequence from OLD to NEW")

    renameg.add_argument("--renamenumber", action="store", dest="renamenumber", metavar="BASENAME",
                          help="Rename all sequences to this form: BASENAME_001, ...")

    renameg.add_argument("--appendnumber", action="store_true", dest="appendnumber",
                          help="Append numbering at end of existing sequence names (SeqA_001, SeqXYZ_002, ...")

    renameg.add_argument("--renameregexp", action="store", dest="renameregexp", metavar='"REGEXP"',
                          help="Rename sequences by deleting parts of names matching regular expression in REGEXP")

    renameg.add_argument("--regdupfix", action="store_true", dest="fixdupnames",
                        help="Fix duplicate names, created by regexp, by appending numbers to duplicates (seqA, seqA_2, ...)")

    renameg.add_argument("--savenames", action="store", dest="savenamefile", metavar="FILE",
                          help="Save renaming information in FILE for later use")

    renameg.add_argument("--restorenames", action="store", dest="restorenames", metavar="FILE",
                          help="Restore original names using info previously saved in FILE")

    renameg.add_argument("--gbname", action="store", dest="gbname", metavar="FIELD1[,FIELD2,FIELD3,...]",
                        help="For Genbank input: construct sequence names from the list of named fields, in the specified order")

    #########################################################################################

    subsetg = parser.add_argument_group("Retrieve subset of sequences")

    subsetg.add_argument("--subsample", action="store", type=int, dest="samplesize", metavar="N",
                        help="Randomly extract N sequences from sequence set")

    subsetg.add_argument("--subset", action="store", dest="namefile", metavar="NAMEFILE",
                        help="Retrieve sequences listed in NAMEFILE")

    subsetg.add_argument("--remseqs", action="store", dest="remfile", metavar="NAMEFILE",
                        help="Discard sequences listed in NAMEFILE")

    subsetg.add_argument("--filterpos", action="store", dest="filterpos", metavar="VARIANT[,VARIANT,...]",
                          help="""Retrieve sequences containing specific residues on specific positions. Syntax is: <POS><RESIDUE>,
                          possibly in a comma-separated list. Example: 484K,501Y""")

    subsetg.add_argument("--filterdupseq", action="store_true", dest="dupseqfilter",
                          help="Remove duplicate sequences (keeping one of each); print names of removed sequences on stderr.")

    subsetg.add_argument("--filterdupname", action="store_true", dest="dupnamefilter",
                          help="Remove sequences with duplicate names (keeping one of each). If this option is not set (default): stop execution on duplicate names.")

    #########################################################################################

    seqpartg = parser.add_argument_group("Extracting or removing parts of sequences")

    seqpartg.add_argument("--subseq", action="store", dest="subseq", metavar="START,STOP",
                          help="Extract subsequence, positions START to STOP, from alignment")

    seqpartg.add_argument("--subseqrename", action="store_true", dest="subseqrename",
                          help="When extracting sub-sequences: add '_START_STOP' to seqnames")

    seqpartg.add_argument("--windows", action="store", type=int, dest="wsize", metavar="WSIZE",
                          help="Extract all overlapping sequence windows of size WSIZE")

    seqpartg.add_argument("--degap", action="store_true", dest="degap",
                          help="Remove all gap characters from sequences")

    seqpartg.add_argument("--remcols", action="store", dest="remcols", metavar="INDEX LIST",
                          help="Remove listed columns from alignment. Columns can be indicated as comma-separated list of indices, and as ranges. Example: --remcols=10,15,22-40,57")

    seqpartg.add_argument("--remambigcols", action="store_true", dest="remambigcols",
                          help="Remove columns where one or more residues are ambiguity symbols (e.g., N for nucleotides)")

    seqpartg.add_argument("--remgapcols", action="store_true", dest="remgapcols",
                          help="Remove columns where one or more residues are gaps")

    seqpartg.add_argument("--remallgapcols", action="store_true", dest="remallgapcols",
                        help="Remove columns that are all-gaps")

    seqpartg.add_argument("--remfracgapcols", action="store", type=float, dest="frac", metavar="FRAC",
                        help="Remove columns that contain > FRAC fraction gaps")

    seqpartg.add_argument("--remconscols", action="store_true", dest="remconscols",
                          help="Remove conserved columns from alignment")

    seqpartg.add_argument("--remhmminsertcols", action="store_true", dest="remhmminsertcols",
                          help="When reading Stockholm format file from HMMer's hmmalign: remove columns corresponding to insert states")

    #########################################################################################

    multifileg = parser.add_argument_group("Combining multiple sequence files")

    multifileg.add_argument("--paste", action="store_true", dest="paste",
                          help="Concatenate identically named sequences from separate input files. " +
                          "Sequences are pasted end to end in the same order as the input files. " +
                          "All input files must contain same number of sequences, and sequences " +
                          "in different files must have same name." +
                          "(To see partitions choose nexus output, or output to multiple partition files).")

    multifileg.add_argument("--overlap", action="store_true", dest="overlap",
                          help="Similar to --paste, but for input alignments that overlap partly. " +
                          "Overlap is discovered automatically and partition boundaries are then set " +
                          "such that each partition is covered by a unique set of genes. " +
                          "(To see partitions choose nexus output, or output to multiple partition files).")

    multifileg.add_argument("--minoverlap", action="store", type=int, dest="minoverlap", metavar="N", default=-1,
                        help="Minimum overlap required for merging input alignments (default: set automatically based on seq lengths)")

    multifileg.add_argument("--multifile", action="store_true", dest="multifile",
                        help="Outputs to multiple files (one per partition) instead of stdout. " +
                          "Partitions are generated automatically based on other options.")

    multifileg.add_argument("--charset", action="store_true", dest="charset",
                        help="Appends Nexus form charset block listing partitions in data (forces output in Nexus format). " +
                          "Charsets and partitions are generated automatically based on other options.")

    multifileg.add_argument("--mbpartblock", action="store_true", dest="mbpartblock",
                        help="Appends MrBayes block with commands for running partitioned analysis (forces output in Nexus format). " +
                          "Charsets and partitions are generated automatically based on other options.")

    #########################################################################################

    dnag = parser.add_argument_group("DNA manipulations")

    dnag.add_argument("--revcomp", action="store_true", dest="revcomp",
                      help="Return reverse complement of sequence(s). Requires sequences to be DNA.")

    dnag.add_argument("--translate", action="store_true", dest="translate",
                      help="Translate DNA into amino acid sequences (requires sequences to be DNA, in frame, and length multiple of 3)")

    #########################################################################################

    summaryg = parser.add_argument_group("Summaries", description="No sequences are printed when these options are used")

    summaryg.add_argument("--num", action="store_true", dest="s_num",
                      help="Print number of sequences")

    summaryg.add_argument("--len", action="store_true", dest="s_len",
                      help="Print summary of sequence lengths")

    summaryg.add_argument("--com", action="store_true", dest="s_com",
                      help="Print overall sequence composition")

    summaryg.add_argument("--seqcom", action="store_true", dest="s_seqcom",
                      help="Print composition for each individual sequence. Output is one line per residue-type per sequence: seqname, residue-type, freq, count, seqlength")

    summaryg.add_argument("--ignoregaps", action="store_true", dest="s_ignoregaps",
                      help="When reporting composition: do not count gap symbols")

    summaryg.add_argument("--nam", action="store_true", dest="s_nam",
                      help="Print names of sequences")

    summaryg.add_argument("--div", action="store_true", dest="s_div",
                      help="(For alignments) Print nucleotide diversity (=average pairwise sequence difference): mean and std")

    summaryg.add_argument("--sit", action="store_true", dest="s_sit",
                      help="""(For alignments) Print site summary: number of columns that are variable (not conserved),
                              number of columns that contain gaps, and number of columns that contain IUPAC ambiguity symbols""")


    #########################################################################################

    debug = parser.add_argument_group("Debugging")

    parser.add_argument("--debug", action="store_true", dest="debug",
                      help="Print longer error messages")

    #########################################################################################

    return parser

################################################################################################

def check_commandline(args):

    if args.informat == "auto":
        args.informat = "autodetect"         # Long name required by seqlib, short better for user interface...

    # Set flags indicating whether input and/or output is aligned, and whether file should be read as alignment
    # Note: args that require manipulation of columns, and --paste also implies aligned sequences
    # (reading as alignment results in check of equal sequence lengths)
    # Note 2: This is somewhat messy. Should really be handled by argparse mutually exclusive arguments
    alignin = alignout = args.aligned = False
    if any([args.informat=="nexus", args.informat=="phylip", args.informat=="clustal", args.informat=="stockholm"]):
        alignin = True
    if any([args.outformat=="nexus", args.outformat=="nexusgap", args.outformat=="phylip", args.outformat=="clustal"]):
        alignout = True
    if any([alignin, alignout, args.remcols, args.remambigcols, args.remgapcols, args.remallgapcols,
            args.remconscols, args.remhmminsertcols, args.subseq, args.paste, args.overlap]):
        args.aligned = True
    if args.frac is not None:        # Note: value may be zero, and bool(0) = False!
        args.aligned = True

    # If option --gbname is set: force input format to "genbank"
    if args.gbname:
        args.informat = "genbank"

    # Perform sanity checks on combination of args
    # (Note: if autodetection of filetype is requested, then many checks are impossible - defer error checking to seqlib)

    # Sanity check #1: option --degap cannot be used with any option that requires aligned sequences
    if args.degap and any([args.subseq, args.remambigcols, args.remgapcols, args.remallgapcols, args.remconscols,
                           args.frac, args.remhmminsertcols, args.mbpartblock, args.charset, args.multifile,
                           args.overlap, args.s_div, args.s_sit, args.paste]):
        raise seqlib.SeqError("Option --degap cannot be used with any option that requires aligned sequences".format(args.outformat))

    # Sanity check #2: option --degap cannot be used if output is an alignment
    if args.degap and alignout:
        raise seqlib.SeqError("Removal of all gap characters (--degap) not possible for output in {}".format(args.outformat))

    # Sanity check #3: option --subseq cannot be used together with args --remambigcols, --remgapcols, --remallgapcols or --remconscols
    if args.subseq and args.remambigcols:
        raise seqlib.SeqError("Options --subseq and --remambigcols cannot be used simultaneously")
    if args.subseq and args.remgapcols:
        raise seqlib.SeqError("Options --subseq and --remgapcols cannot be used simultaneously")
    if args.subseq and args.remallgapcols:
        raise seqlib.SeqError("Options --subseq and --remallgapcols cannot be used simultaneously")
    if args.subseq and args.remconscols:
        raise seqlib.SeqError("Options --subseq and --remconscols cannot be used simultaneously")

    # Sanity check #4: option --subseqrename can only be used in combination with --subseq
    if args.subseqrename and not args.subseq:
        raise seqlib.SeqError("Option --subseqrename (add '_START_STOP' to seqnames) requires option --subseq")

    # Sanity check #5: option --savenames requires option --renamenumber or --renameregexp
    if args.savenamefile and not (args.renamenumber or args.renameregexp):
        raise seqlib.SeqError("Option --savenames requires option --renamenumber or --renameregexp")

    # Sanity check #6: args --renamenumber and --restorenames are incompatible
    if args.renamenumber and  args.restorenames:
        raise seqlib.SeqError("Option --renamenumber and --restorenames cannot be used simultaneously")

    # Sanity check #7: args --renameregexp and --restorenames are incompatible
    if args.renameregexp and  args.restorenames:
        raise seqlib.SeqError("Option --renameregexp and --restorenames cannot be used simultaneously")

    # Sanity check #8: args --renameregexp and --paste are incompatible
    if args.renameregexp and  args.paste:
        raise seqlib.SeqError("Option --renameregexp and --paste cannot be used simultaneously: first rename, then paste resulting files")

    # Sanity check #9c: option --overlap is not meaningful unless several input files are provided
    if args.overlap and len(args.filelist) == 1:
        raise seqlib.SeqError("Option --overlap requires multiple input files")

    return (args)

################################################################################################

def read_seqs(args):

    for filename in args.filelist:
        if filename != "-" and not os.path.isfile(filename):
            raise seqlib.SeqError("File %s not found." % filename)

    # Read sequences from all files
    seqlist = []          # List of either Seq_set or Seq_alignment objects
    for filename in args.filelist:
        if args.gbname:
            seqfile = seqlib.Genbankfilehandle(filename, namefromfields=args.gbname)
        else:
            seqfile = seqlib.Seqfile(filename, args.informat)

        # If filetype was autodetected: check if this should be read as an alignment:
        if args.informat == "autodetect" and isinstance(seqfile, seqlib.Alignfile_reader):
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
                newseqset = seqlib.Seq_alignment()
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
        ra = seqlib.Read_assembler( consensus_list )
        contiglist = ra.assemble(minoverlap=args.minoverlap)

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
                    seqs = seqlib.Seq_alignment(new_name)
                    seqs.addseqset(subalignment)
                else:
                    seqs.appendalignment(subalignment)

    # For other options than overlap: Collect sequences in one sequence collection object
    else:
        seqs = seqlist[0]
        for seqset in seqlist[1:]:
            if args.paste:
                seqs.appendalignment(seqset)
            else:
                seqs.addseqset(seqset, args.dupnamefilter)

    return seqs, args

################################################################################################

def check_args_alignment_issues(args):
    # After read_seqs we know whether sequences correspond to alignment.
    # Check if any requested option clashes with alignment status
    if not args.aligned:
        bad_options_list = []
        if args.subseq:
            bad_options_list.append("--subseq")
        if args.remambigcols:
            bad_options_list.append("--remambigcols")
        if args.remgapcols:
            bad_options_list.append("--remgapcols")
        if args.remallgapcols:
            bad_options_list.append("--remallgapcols")
        if args.frac:
            bad_options_list.append("--remfracgapcols")
        if args.remconscols:
            bad_options_list.append("--remconscols")
        if args.remhmminsertcols:
            bad_options_list.append("--remhmminsertcols")
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
            raise SeqError("The option(s) below require aligned sequences:\n\t{}".format(bad_option_string))

################################################################################################

def filterdupseqs(seqs):
    simlist = seqs.removedupseqs()
    if simlist:
        sys.stderr.write("Lists of identical sequences:\n")
        for namelist in simlist:
            sys.stderr.write(" ".join(namelist))
            sys.stderr.write("\n")

    return(seqs)

################################################################################################

def positionfilter(seqs, args):

    # Extract positions and residues from input of the form: 484K,501Y,987W.
    # Save as list of (pos, residue) tuples
    varlist = []
    items = args.filterpos.split(",")
    for item in items:
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

    # Removal of gaps from individual sequences
    if args.degap:
        # If sequences have been read as alignment: first convert to Seq_set:
        if args.aligned:
            args.aligned = False
            newseqs = seqlib.Seq_set()
            for seq in seqs:
                newseqs.addseq(seq)
            seqs = newseqs
        seqs.remgaps()

    # Removal of listed columns
    if args.remcols:
        remlist = []
        items = args.remcols.split(",")
        for item in items:
            if "-" in item:
                (start,stop) = item.split("-")
                start = int(start) - 1          # Slice syntax: Numbering starts at 0
                stop = int(stop)                # Slice syntax: stop is one more than last index (=> do not need to subtract 1 from non-slice syntax)
                remlist.extend(list(range(start,stop)))
            else:
                remlist.append(int(item) - 1)
        seqs.remcols(remlist)

    # Removal of ambiguity columns
    if args.remambigcols:
        seqs.remambigcol()

    # Removal of gappy columns
    if args.remgapcols:
        seqs.remgapcol()

    # Removal of all-gap columns
    if args.remallgapcols:
        seqs.remallgapcol()

    # Removal of columns with > FRAC fraction gaps
    if args.frac is not None:
        seqs.remfracgapcol(args.frac)

    # Removal of conserved columns
    if args.remconscols:
        seqs.remconscol()

    # Removal of insert state columns (in output from HMMer's hmmalign)
    if args.remhmminsertcols:
        seqs.rem_hmmalign_insertcol()

    # Extraction of subsequence
    # Note: it is assumed that indexing starts at 1, and that stop is included ("slicesyntax=False")
    if args.subseq:
        (start, stop) = args.subseq.split(",")
        start = int(start)
        stop = int(stop)
        seqs = seqs.subseq(start, stop, slicesyntax=False, rename=args.subseqrename)

    # Extraction of sequence windows
    if args.wsize:
        newseqs = seqlib.Seq_set()
        for seq in seqs:
            for seqwin in seq.windows(args.wsize, rename=True):
                newseqs.addseq(seqwin)
        seqs = newseqs

    # Renaming
    if args.rename:
        (oldname, newname) = args.rename.split(",")
        seqs.changeseqname(oldname, newname)
    elif args.renamenumber:
        seqs.rename_numbered(args.renamenumber, args.savenamefile)
    elif args.renameregexp:
        seqs.rename_regexp(args.renameregexp, args.savenamefile, args.dupnamefilter, args.fixdupnames)
    # Note: option appendnumber is not under "else" clause: In principle it can be combined with previous renaming!
    if args.appendnumber:
        seqs.appendnumber(args.savenamefile)

    # Restore names
    if args.restorenames:
        seqs.transname(args.restorenames)


    # Translation
    if args.translate:
        seqs = seqs.translate()

    # Reverse complement
    if args.revcomp:
        seqs = seqs.revcomp()

    # Return altered sequences
    return seqs

################################################################################################

def print_summary(seqs, args):

    if args.s_nam:
        print("# Sequence names:")
        for name in sorted(seqs.seqnamelist):
            print(name)
        print()

    if args.s_num:
        print("Number of sequences: {:10d}".format(len(seqs)))

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
        avg, std = seqs.sequence_diversity(ignoregaps=False)
        print("Nucleotide diversity pi (average pairwise sequence distance)")
        print("    Mean:               {:.5f}".format(avg))
        print("    Standard dev:       {:.5f}\n".format(std))

    if args.s_sit:
        numvar = 0
        numgap = 0
        numambig = 0
        for col in seqs.columns():
            columnset = set(col)
            if len(columnset) != 1:        # nvalues == 1 <=> conserved column
                numvar += 1
            if "-" in columnset:
                numgap += 1
            if (seqs.ambigsymbols & columnset):
                numambig += 1
        print("Site summary (note: variable, gappy, and ambiguous sites may overlap)")
        print("    No. variable sites:  {:>6d}".format(numvar))
        print("    No. gappy sites:     {:>6d}".format(numgap))
        print("    No. ambiguous sites: {:>6d}\n".format(numambig))


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

    if args.s_seqcom:
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
    elif args.outformat == "nexusgap":
        filehandle.write(seqs.nexusgap(width = args.width) + '\n')
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
            raise seqlib.SeqError("File with name {} already exists. Aborting to avoid overwriting.".format(outname))
        fh = open(outname, "w")
        print_seqs(alignment, args, filehandle=fh)
        fh.close()

################################################################################################

if __name__ == "__main__":
    main()
