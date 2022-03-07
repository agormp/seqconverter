#!/usr/bin/env python3

# Anders Gorm Pedersen, agpe@dtu.dk, 2012-2021
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
        seqs = read_seqs(args)

        if args.dupseqfilter:
            seqs = filterdupseqs(seqs)

        if args.filterpos:
            seqs = positionfilter(seqs, args)

        seqs = change_seqs(seqs, args)

        if args.multifile:
            # Output to multiple files: one partition per file, in nexus format
            write_partitions(seqs, args)
        else:
            # Output of one single sequence file to stdout.
            if args.summary:
                print_summary(seqs, args)
            elif args.summarynames:
                for name in sorted(seqs.seqnamelist):
                    print(name)
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

    formatg.add_argument("--nocomments", action="store_true", dest="nocomments",
                        help="Do not include comments in output (only print seqnames)")

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
                          help="For output from HMMer's hmmalign: remove columns corresponding to insert states")

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

    summaryg = parser.add_argument_group("Summaries")

    summaryg.add_argument("--summary", action="store_true", dest="summary",
                      help="Print summary of data set (names, number, lengths, composition, etc.). No sequences are output.")

    summaryg.add_argument("--names", action="store_true", dest="summarynames",
                      help="Print names of sequences in data set.")

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

    # Sanity check #1: option --degap cannot be used together with args --remgapcols, --remallgapcols, or --remconscols
    if args.degap and args.remgapcols:
        raise seqlib.SeqError("Options --degap and --remgapcols cannot be used simultaneously")
    if args.degap and args.remallgapcols:
        raise seqlib.SeqError("Options --degap and --remallgapcols cannot be used simultaneously")
    if args.degap and args.remconscols:
        raise seqlib.SeqError("Options --degap and --remconscols cannot be used simultaneously")

    # Sanity check #2: option --degap cannot be used if output is an alignment
    if args.degap and alignout:
        raise seqlib.SeqError("Removal of all gap characters (--degap) cannot be performed on aligned sequences")

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
            seqfile = seqlib.Genbankfilehandle(filename, namefromfeatures=args.gbname)
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

    return seqs

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

    # Print number of seqs
    print("Number of sequences: {:10d}".format(len(seqs)))

    # Unaligned sequences
    if not args.aligned:
        # Print minimum, maximum, and average lengths
        seqlenlist = []
        lensum = 0
        for seq in seqs:
            seqlen = len(seq)
            seqlenlist.append(seqlen)
            lensum += seqlen
        print("Minimum length:   {:13d}".format(min(seqlenlist)))
        print("Maximum length:   {:13d}".format(max(seqlenlist)))
        print("Average length:   {:13.2f}\n\n".format(lensum/len(seqs)))

    # Alignments:
    else:
        # Print alignment length
        print("Alignment length:    {:10d}\n".format( seqs.alignlen() ))

        # Print distance summary: nucleotide diversity (= average pairwise distance) + its standard error
        avg, sem = seqs.sequence_diversity()
        print("Nucleotide diversity pi (average pairwise sequence distance)")
        print("    Mean:               {:.5f}".format(avg))
        print("    Standard error:     {:.5f}\n".format(sem))

        # Print number of variable sites, number of gapcontaining sites, and number of ambiguity containing sites
        numvar = 0
        numgap = 0
        numambig = 0
        for i, col in enumerate(seqs.columns()):
            columnset = set(col)
            if len(columnset) != 1:        # nvalues == 1 <=> conserved column
                numvar += 1
            if "-" in columnset:
                numgap += 1
            if (seqs.ambigsymbols & columnset):
                numambig += 1
        print("Site summary (note: variable, gappy, and ambiguous sites may overlap)")
        print("    No. variable sites:  {:5d}".format(numvar))
        print("    No. gappy sites:     {:5d}".format(numgap))
        print("    No. ambiguous sites: {:5d}\n".format(numambig))

    # Print composition #########################################################################################
    compositiondict = seqs.composition()

    # If any ambiguity symbols have zero counts: remove from output
    remset = set()
    for symbol in seqs.ambigsymbols:
        if compositiondict[symbol][0] == 0:
            remset.add(symbol)
    symbols = compositiondict.keys() - remset
    symbols = sorted(list(symbols))

    print("Composition:")
    print("    {:^7s}{:^8s}{:^8s}".format("Symbol", "Freq", "Count"))
    tot_count = 0
    for symbol in symbols:
        tot_count += compositiondict[symbol][0]
    for symbol in symbols:
        print("    {:^7s}{:^8.3f}{:^8d} / {}".format(symbol, compositiondict[symbol][1], compositiondict[symbol][0], tot_count))

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
        filehandle.write(seqs.tab( nocomments=args.nocomments ) + '\n')
    elif args.outformat == "fasta":
        filehandle.write(seqs.fasta( nocomments=args.nocomments ) + '\n')
    elif args.outformat == "how":
        filehandle.write(seqs.how( nocomments=args.nocomments ) + '\n')
    elif args.outformat == "nexus":
        if any([args.paste, args.overlap, args.charset]) and (not args.multifile):
            parts = True
        else:
            parts = False
        filehandle.write(seqs.nexus(print_partitioned = parts) + '\n')
    elif args.outformat == "nexusgap":
        filehandle.write(seqs.nexusgap() + '\n')
    elif args.outformat == "phylip":
        filehandle.write(seqs.phylip() + '\n')
    elif args.outformat == "clustal":
        filehandle.write(seqs.clustal() + '\n')

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
