#!/usr/bin/env python3
# Anders Gorm Pedersen, agpe@dtu.dk,
# Section for Bioinformatics, DTU Health Technology, Technical University of Denmark

# Converts between different sequence file formats. Performs various manipulations on sequences

import sys
import os.path
import sequencelib as seqlib
#import statistics as st NOTE: only works on 3.4 and upwards. Wait for CBS compatibility...
from optparse import OptionParser

################################################################################################

def build_parser():
    parser = OptionParser(usage="usage: %prog [options] [SEQFILE [SEQFILE...]]",
                          version="2.1")

    parser.add_option("-I", type="choice", dest="informat",
                      choices=["auto", "nexus", "phylip", "fasta", "clustal", "raw", "tab", "genbank", "how"], metavar="FORMAT",
                      help="Input format: auto, fasta, tab, raw, genbank, how, clustal, phylip, nexus [default: auto]")

    parser.add_option("-O", type="choice", dest="outformat",
                      choices=["nexus", "phylip", "fasta", "clustal", "raw", "tab", "nexusgap", "how"], metavar="FORMAT",
                      help="Output format: fasta, tab, raw, how, clustal, phylip, nexus, nexusgap, nexusmulti [default: fasta]")

    parser.add_option("--nocomments", action="store_true", dest="nocomments",
                        help="Do not print comments in output (only print seqnames)")

    parser.add_option("--subseq", action="store", type="string", dest="subseq", metavar="START,STOP",
                          help="Extract subsequence, positions START to STOP, from alignment")

    parser.add_option("--subseqrename", action="store_true", dest="subseqrename",
                          help="When extracting sub-sequences: add '_START_STOP' to seqnames")

    parser.add_option("--windows", action="store", type="int", dest="wsize", metavar="WSIZE",
                          help="Extract all overlapping sequence windows of size WSIZE")

    parser.add_option("--subsample", action="store", type="int", dest="samplesize", metavar="N",
                        help="Randomly extract N sequences from sequence set")

    parser.add_option("--subset", action="store", type="string", dest="namefile", metavar="NAMEFILE",
                        help="Extract sequences listed in NAMEFILE")

    parser.add_option("--remseqs", action="store", type="string", dest="remfile", metavar="NAMEFILE",
                        help="Remove sequences listed in NAMEFILE")

    parser.add_option("--filterpos", action="store", type="string", dest="filterpos", metavar="VARIANT[,VARIANT,...]",
                          help="Extract sequences containing specific residues on specific positions. Syntax is: <POS><RESIDUE>, possibly in a comma-separated lidt. Example: 484K,501Y")

    parser.add_option("--gff", action="store", type="string", dest="gff", metavar="FILE",
                          help="Get annotation info for sequences from GFF-formatted FILE")

    parser.add_option("--gffsymbol", action="store", type="string", dest="gffsymbol", metavar="ANNOT",
                          help="Use the character ANNOT in the sequence annotation strings derived from GFF file")

    parser.add_option("--degap", action="store_true", dest="degap",
                          help="Remove all gap characters from sequences")

    parser.add_option("--remcols", action="store", type="string", dest="remcols", metavar="INDEX LIST",
                          help="Remove listed columns from alignment. Columns can be indicated as comma-separated list of indices, and as ranges. Example: --remcols=10,15,22-40,57")

    parser.add_option("--remambigcols", action="store_true", dest="remambigcols",
                          help="Remove columns where one or more residues are ambiguity symbols (e.g., N for nucleotides)")

    parser.add_option("--remgapcols", action="store_true", dest="remgapcols",
                          help="Remove columns where one or more residues are gaps")

    parser.add_option("--remallgapcols", action="store_true", dest="remallgapcols",
                        help="Remove columns that are all-gaps")

    parser.add_option("--remfracgapcols", action="store", type="float", dest="frac", metavar="FRAC",
                        help="Remove columns that contain > FRAC fraction gaps")

    parser.add_option("--remconscols", action="store_true", dest="remconscols",
                          help="Remove conserved columns from alignment")

    parser.add_option("--paste", action="store_true", dest="paste",
                          help="Concatenate identically named sequences from separate input files. " +
                          "Sequences are pasted end to end in the same order as the input files. " +
                          "All input files must contain same number of sequences, and sequences " +
                          "in different files must have same name." +
                          "(To see partitions choose nexus output, or output to multiple partition files).")

    parser.add_option("--overlap", action="store_true", dest="overlap",
                          help="Similar to --paste, but for input alignments that overlap partly. " +
                          "Overlap is discovered automatically and partition boundaries are then set " +
                          "such that each partition is covered by a unique set of genes. " +
                          "(To see partitions choose nexus output, or output to multiple partition files).")

    parser.add_option("--minoverlap", action="store", type="int", dest="minoverlap", metavar="N",
                        help="Minimum overlap required for merging input alignments (default: set automatically based on seq lengths)")

    parser.add_option("--multifile", action="store_true", dest="multifile",
                        help="Outputs to multiple files (one per partition) instead of stdout. " +
                          "Partitions are generated automatically based on other options.")

    parser.add_option("--charset", action="store_true", dest="charset",
                        help="Appends Nexus form charset block listing partitions in data (forces output in Nexus format). " +
                          "Charsets and partitions are generated automatically based on other options.")

    parser.add_option("--mbpartblock", action="store_true", dest="mbpartblock",
                        help="Appends MrBayes block with commands for running partitioned analysis (forces output in Nexus format). " +
                          "Charsets and partitions are generated automatically based on other options.")

    parser.add_option("--bestblock", action="store_true", dest="bestblock",
                          help="Appends MrBayes block with commands for running BEST analysis to output file (forces output in Nexus format). " +
                            "Charsets and partitions are generated automatically so they correspond to input files " +
                            "(one partition for each file).")

    parser.add_option("--filterdupseq", action="store_true", dest="dupseqfilter",
                          help="Remove duplicate sequences (keeping one of each); print names of removed sequences on stderr.")

    parser.add_option("--filterdupname", action="store_true", dest="dupnamefilter",
                          help="Remove sequences with duplicate names (keeping one of each). If this option is not set (default): stop execution on duplicate names.")

    parser.add_option("--gbname", action="store", type="string", dest="gbname", metavar="FIELD1[,FIELD2,FIELD3,...]",
                        help="For Genbank input: construct sequence names from the list of named fields, in the specified order")

    parser.add_option("--appendnumber", action="store_true", dest="appendnumber",
                          help="Append numbering at end of existing sequence names (SeqA_001, SeqXYZ_002, ...")

    parser.add_option("--rename", action="store", type="string", dest="rename", metavar="OLD,NEW",
                          help="Rename sequence from OLD to NEW")

    parser.add_option("--renamenumber", action="store", type="string", dest="renamenumber", metavar="BASENAME",
                          help="Rename sequences to this form: BASENAME_001, ...")

    parser.add_option("--renameregexp", action="store", type="string", dest="renameregexp", metavar='"REGEXP"',
                          help="Rename sequences by deleting parts of names matching regular expression in REGEXP")

    parser.add_option("--regdupfix", action="store_true", dest="fixdupnames",
                        help="Fix duplicate names, created by regexp, by appending numbers to duplicates (seqA, seqA_2, ...)")

    parser.add_option("--savenames", action="store", type="string", dest="savenamefile", metavar="FILE",
                          help="Save renaming information in FILE for later use")

    parser.add_option("--restorenames", action="store", type="string", dest="restorenames", metavar="FILE",
                          help="Restore original names using info previously saved in FILE")

    parser.add_option("--revcomp", action="store_true", dest="revcomp",
                      help="Return reverse complement of sequence(s).")

    parser.add_option("--translate", action="store_true", dest="translate",
                      help="Translate DNA into amino acid sequences (requires sequences to be DNA, in frame, and length multiple of 3)")

    parser.add_option("--summary", action="store_true", dest="summary",
                      help="Print summary of data set (names, number, lengths, composition, etc.). No sequences are output.")

    parser.add_option("--names", action="store_true", dest="summarynames",
                      help="Print names of sequences in data set.")

    parser.add_option("--debug", action="store_true", dest="debug",
                      help="Print longer error messages")

    parser.set_defaults(informat="auto", outformat="fasta", nocomments=False, degap=False, remcols=None, remambigcols=False, remgapcols=False, remallgapcols=False,
                        frac=None, remconscols=False, paste=False, overlap=False, minoverlap=-1,
                        charset=False, multifile=False, mbpartblock=False, bestblock=False, gff=None, gffsymbol=None,
                        dupnamefilter=False, fixdupnames=False, dupseqfilter=False, subseq=None, wsize=None, samplesize=0, savenamefile=None, subseqrename=False,
                        gbname=None, appendnumber=False, rename=None, renamenumber=None, renameregexp=False, namefile=None, remfile=None, filterpos=None,
                        restorenames=None, revcomp=False, translate=False, summary=False, summarynames=False, debug=False)

    return parser

################################################################################################

def check_commandline(options, args):

    if options.informat == "auto":
        options.informat = "autodetect"         # Long name required by seqlib, short better for user interface...

    # Set flags indicating whether input and/or output is aligned, and whether file should be read as alignment
    # Note: options that require manipulation of columns, and --paste and --bestblock also implies aligned sequences
    # (reading as alignment results in check of equal sequence lengths)
    alignin = alignout = options.aligned = False
    if (options.informat == "nexus" or options.informat == "phylip" or options.informat == "clustal"):
        alignin = True
    if (options.outformat == "nexus" or options.outformat == "nexusgap" or options.outformat == "phylip" or options.outformat == "clustal"):
        alignout = True
    if (alignin or alignout or options.remcols or options.remambigcols or options.remgapcols or options.remallgapcols or
                            options.remconscols or options.subseq or options.paste or options.bestblock or
                            options.overlap):
        options.aligned = True
    if options.frac is not None:        # Note: value may be zero, and bool(0) = False!
        options.aligned = True

    # If option --gbname is set: force input format to "genbank"
    if options.gbname:
        options.informat = "genbank"

    # If options.frac is set: convert from string to float
    if options.frac is not None:
        options.frac = float(options.frac)

    # Perform sanity checks on combination of options
    # (Note: if autodetection of filetype is requested, then many checks are impossible - defer error checking to seqlib)

    # Sanity check #1: option --degap cannot be used together with options --remgapcols, --remallgapcols, or --remconscols
    if options.degap and options.remgapcols:
        raise seqlib.SeqError("Options --degap and --remgapcols cannot be used simultaneously")
    if options.degap and options.remallgapcols:
        raise seqlib.SeqError("Options --degap and --remallgapcols cannot be used simultaneously")
    if options.degap and options.remconscols:
        raise seqlib.SeqError("Options --degap and --remconscols cannot be used simultaneously")

    # Sanity check #2: option --degap cannot be used if output is an alignment
    if options.degap and alignout:
        raise seqlib.SeqError("Removal of all gap characters (--degap) cannot be performed on aligned sequences")

    # Sanity check #3: option --subseq cannot be used together with options --remambigcols, --remgapcols, --remallgapcols or --remconscols
    if options.subseq and options.remambigcols:
        raise seqlib.SeqError("Options --subseq and --remambigcols cannot be used simultaneously")
    if options.subseq and options.remgapcols:
        raise seqlib.SeqError("Options --subseq and --remgapcols cannot be used simultaneously")
    if options.subseq and options.remallgapcols:
        raise seqlib.SeqError("Options --subseq and --remallgapcols cannot be used simultaneously")
    if options.subseq and options.remconscols:
        raise seqlib.SeqError("Options --subseq and --remconscols cannot be used simultaneously")

    # Sanity check #4: option --subseqrename can only be used in combination with --subseq
    if options.subseqrename and not options.subseq:
        raise seqlib.SeqError("Option --subseqrename (add '_START_STOP' to seqnames) requires option --subseq")

    # Sanity check #5: option --savenames requires option --renamenumber or --renameregexp
    if options.savenamefile and not (options.renamenumber or options.renameregexp):
        raise seqlib.SeqError("Option --savenames requires option --renamenumber or --renameregexp")

    # Sanity check #6: options --renamenumber and --restorenames are incompatible
    if options.renamenumber and  options.restorenames:
        raise seqlib.SeqError("Option --renamenumber and --restorenames cannot be used simultaneously")

    # Sanity check #7: options --renameregexp and --restorenames are incompatible
    if options.renameregexp and  options.restorenames:
        raise seqlib.SeqError("Option --renameregexp and --restorenames cannot be used simultaneously")

    # Sanity check #8: options --renameregexp and --paste are incompatible
    if options.renameregexp and  options.paste:
        raise seqlib.SeqError("Option --renameregexp and --paste cannot be used simultaneously: first rename, then paste resulting files")

    # Sanity check #9: option --bestblock is not meaningful unless several input files are provided
    if options.bestblock and len(args) == 1:
        raise seqlib.SeqError("Option --bestblock requires multiple input files (one for each partition)")

    # Sanity check #9c: option --overlap is not meaningful unless several input files are provided
    if options.overlap and len(args) == 1:
        raise seqlib.SeqError("Option --overlap requires multiple input files")

    # Sanity check #10: option --gffsymbol requires option --gff
    if options.gffsymbol and not options.gff:
        raise seqlib.SeqError("Option --gffsymbol requires option --gff")

    return (options, args)

################################################################################################

def read_seqs(options, args):

    # Check filenames (if no filename is given: assume stdin)
    if len(args) < 1:
        filenames = ["-"]
    else:
        filenames = args
        for filename in filenames:
            if not os.path.isfile(filename):
                raise seqlib.SeqError("File %s not found." % filename)

    # Read sequences from all files
    seqlist = []          # List of either Seq_set or Seq_alignment objects
    for filename in filenames:
        if options.gbname:
            seqfile = seqlib.Genbankfilehandle(filename, namefromfeatures=options.gbname)
        else:
            seqfile = seqlib.Seqfile(filename, options.informat)

        # If filetype was autodetected: check if this should be read as an alignment:
        if options.informat == "autodetect" and isinstance(seqfile, seqlib.Alignfile_reader):
            options.aligned = True

        if options.aligned:
            seqlist.append(seqfile.read_alignment(options.dupnamefilter))
        else:
            seqlist.append(seqfile.read_seqs(options.dupnamefilter))

    # Post processing to check if this is alignment, despite other flags not indicating it
    # This can happen if alignment is in e.g. fasta format, and none of the alignment status-triggering flags are used
    # Assume sequences are aligned if all sequences have same length:
    if not options.aligned:
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
            options.aligned = True

    # If automatic overlap detection has been requested: find overlaps and concatenate or merge sub-alignments
    if options.overlap:

        # Create consensus sequence from each seq in seqset
        consensus_list = []
        alignment_dict = {}
        for alignment in seqlist:
            consensus_list.append( alignment.consensus() )
            alignment_dict[alignment.name] = alignment

        # Find overlap between consensus sequences, get list of contigs
        ra = seqlib.Read_assembler( consensus_list )
        contiglist = ra.assemble(minoverlap=options.minoverlap)

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
            if options.paste or options.bestblock:
                seqs.appendalignment(seqset)
            else:
                seqs.addseqset(seqset, options.dupnamefilter)

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

def positionfilter(seqs, options):

    # Extract positions and residues from input of the form: 484K,501Y,987W.
    # Save as list of (pos, residue) tuples
    varlist = []
    items = options.filterpos.split(",")
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

def change_seqs(seqs, options):

    # Extract random subsample of sequences
    if options.samplesize:
        seqs = seqs.subsample(options.samplesize)

    # Extract sequences named in options.namefile
    if options.namefile:
        namelist = []
        namef = open(options.namefile)
        for line in namef:
            namelist.append(line.strip())
        namef.close()
        seqs = seqs.subset(namelist)

    # Remove sequences named in options.remfile
    if options.remfile:
        namelist = []
        namef = open(options.remfile)
        for line in namef:
            namelist.append(line.strip())
        namef.close()
        seqs.remseqs(namelist)

    # Removal of gaps from individual sequences
    if options.degap:
        # If sequences have been read as alignment: first convert to Seq_set:
        if options.aligned:
            newseqs = seqlib.Seq_set()
            for seq in seqs:
                newseqs.addseq(seq)
            seqs = newseqs
        seqs.remgaps()

    # Removal of listed columns
    if options.remcols:
        remlist = []
        items = options.remcols.split(",")
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
    if options.remambigcols:
        seqs.remambigcol()

    # Removal of gappy columns
    if options.remgapcols:
        seqs.remgapcol()

    # Removal of all-gap columns
    if options.remallgapcols:
        seqs.remallgapcol()

    # Removal of columns with > FRAC fraction gaps
    if options.frac is not None:
        seqs.remfracgapcol(options.frac)

    # Removal of conserved columns
    if options.remconscols:
        seqs.remconscol()

    # Extraction of subsequence
    # Note: it is assumed that indexing starts at 1, and that stop is included ("slicesyntax=False")
    if options.subseq:
        (start, stop) = options.subseq.split(",")
        start = int(start)
        stop = int(stop)
        seqs = seqs.subseq(start, stop, slicesyntax=False, rename=options.subseqrename)

    # Extraction of sequence windows
    if options.wsize:
        newseqs = seqlib.Seq_set()
        for seq in seqs:
            for seqwin in seq.windows(options.wsize, rename=True):
                newseqs.addseq(seqwin)
        seqs = newseqs

    # Renaming
    if options.rename:
        (oldname, newname) = options.rename.split(",")
        seqs.changeseqname(oldname, newname)
    elif options.renamenumber:
        seqs.rename_numbered(options.renamenumber, options.savenamefile)
    elif options.renameregexp:
        seqs.rename_regexp(options.renameregexp, options.savenamefile, options.dupnamefilter, options.fixdupnames)
    # Note: option appendnumber is not under "else" clause: In principle it can be combined with previous renaming!
    if options.appendnumber:
        seqs.appendnumber(options.savenamefile)

    # Restore names
    if options.restorenames:
        seqs.transname(options.restorenames)


    # Translation
    if options.translate:
        seqs = seqs.translate()

    # Reverse complement
    if options.revcomp:
        seqs = seqs.revcomp()

    # Return altered sequences
    return seqs

################################################################################################

def print_summary(seqs, options):

    # Print number of seqs
    print("Number of sequences: {:10d}".format(len(seqs)))

    # Unaligned sequences
    if not options.aligned:
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
        dm = seqs.distmatrix(dist="pdist_ignoregaps")
        dist_summary = dm.distsummary()
        print("Nucleotide diversity pi (average pairwise sequence distance)")
        print("    Mean:               {:.5f}".format(dist_summary[0]))
        print("    Standard error:     {:.5f}\n".format(dist_summary[1]))

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
        print("Site summary (note: variable, gappy, and ambiguous sites may overlap")
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

def print_seqs(seqs, options, filehandle=sys.stdout):

    # NOTE: default is to write to stdout (print to terminal) but other filehandle can be set

    # Some options implies Nexus format output
    if options.charset or options.bestblock or options.mbpartblock:
        options.outformat = "nexus"

    # Print sequences in requested format
    if options.outformat == "raw":
        filehandle.write(seqs.raw() + '\n')
    elif options.outformat == "tab":
        filehandle.write(seqs.tab( nocomments=options.nocomments ) + '\n')
    elif options.outformat == "fasta":
        filehandle.write(seqs.fasta( nocomments=options.nocomments ) + '\n')
    elif options.outformat == "how":
        filehandle.write(seqs.how( nocomments=options.nocomments ) + '\n')
    elif options.outformat == "nexus":
        if (options.paste or options.overlap or options.bestblock or options.charset) and (not options.multifile):
            parts = True
        else:
            parts = False
        filehandle.write(seqs.nexus(print_partitioned = parts) + '\n')
    elif options.outformat == "nexusgap":
        filehandle.write(seqs.nexusgap() + '\n')
    elif options.outformat == "phylip":
        filehandle.write(seqs.phylip() + '\n')
    elif options.outformat == "clustal":
        filehandle.write(seqs.clustal() + '\n')

################################################################################################

def write_partitions(seqs, options):
    part_alignments = seqs.partitions_as_seqalignments()
    for alignment in part_alignments:
        outname = "partition_{alname}.{suffix}".format(alname=alignment.name, suffix=options.outformat)
        if os.path.isfile(outname):
            raise seqlib.SeqError("File with name {} already exists. Aborting to avoid overwriting.".format(outname))
        fh = open(outname, "w")
        print_seqs(alignment, options, filehandle=fh)
        fh.close()

################################################################################################

def use_gff_file(seqs, options):
    """Reads file with GFF information, extracts annotation info, and adds to corresponding sequence objects"""

    # Check that GFF file exists
    if not os.path.isfile(options.gff):
        raise seqlib.SeqError("File %s not found." % options.gff)

    # Construct empty annotation lists for all sequences in sequence set
    for seq in seqs:
        seq.annotlist = ["."] * len(seq)
        seq.comments += "'.'=Nothing; "
        seq.annotsymbols = {"."}

    # Set annotation symbol
    if options.gffsymbol:
        annotsymbol = options.gffsymbol
    else:
        annotsymbol = "A"

    # Parse GFF file, add annotation to sequence objects along the way
    for line in open(options.gff):
        words = line.split()
        seqname = words[0]
        annot = words[2]
        start = int(words[3])
        stop = int(words[4])
        seq = seqs.getseq(seqname)
        seq.annotlist[start - 1:stop] = [annotsymbol] * (stop - start + 1)
        if annotsymbol not in seq.annotsymbols:
            seq.comments += "'{}'={}; ".format(annotsymbol, annot)
            seq.annotsymbols.add(annotsymbol)

    # Postprocess sequence set: construct seq.annotation strings from seq.annotlist lists
    for seq in seqs:
        seq.annotation = "".join(seq.annotlist)

################################################################################################

def main():
    parser = build_parser()
    (options, args) = parser.parse_args()
    try:
        (options, args) = check_commandline(options, args)
        seqs = read_seqs(options, args)

        if options.gff:
            use_gff_file(seqs, options)

        if options.dupseqfilter:
            seqs = filterdupseqs(seqs)

        if options.filterpos:
            seqs = positionfilter(seqs, options)

        seqs = change_seqs(seqs, options)

        if options.multifile:
            # Output to multiple files: one partition per file, in nexus format
            write_partitions(seqs, options)
        else:
            # Output of one single sequence file to stdout.
            if options.summary:
                print_summary(seqs, options)
            elif options.summarynames:
                for name in sorted(seqs.seqnamelist):
                    print(name)
            else:
                print_seqs(seqs, options)

            if options.bestblock:
                print("")
                print(seqs.bestblock())

            if options.charset:
                print("")
                print(seqs.charsetblock())

            if options.mbpartblock:
                print("")
                print(seqs.mbpartblock())


    except seqlib.SeqError as exc:
        if options.debug:
            import traceback
            traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write("Error: {}\n".format(exc.errormessage))

################################################################################################

if __name__ == "__main__":
    main()
