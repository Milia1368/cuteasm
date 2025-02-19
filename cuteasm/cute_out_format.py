import pickle
import time,re,os,logging
from collections import defaultdict
from cute_candidate import *


def sorted_nicely(vcf_entries):
    """ Sort the given vcf entries (in the form ((contig, start, end), vcf_string, sv_type)) in the way that humans expect.
        e.g. chr10 comes after chr2
        Algorithm adapted from https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/"""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    tuple_key = lambda entry: ( alphanum_key(str(entry[0][0])), entry[0][1], entry[0][2] )
    return sorted(vcf_entries, key = tuple_key)


def write_final_vcf(candidates,
                    contig_names, 
                    contig_lengths,
                    reference,
                    argv,options):
    if not options.work_dir.endswith('/'):
        options.work_dir += '/'

    file_out = open(options.work_dir + options.output , 'w')


    # General header
    file_out.write("##fileformat=VCFv4.2\n")
	#file.write("##source=cuteSV-%s\n"%(VERSION))
    import time
    file_out.write("##fileDate=%s\n"%(time.strftime('%Y-%m-%d %H:%M:%S %w-%Z',time.localtime())))
    for contig_name, contig_length in zip(contig_names, contig_lengths):
        print("##contig=<ID={0},length={1}>".format(contig_name, contig_length), file=file_out)

	# Specific header
	# ALT
    file_out.write("##ALT=<ID=INS,Description=\"Insertion of novel sequence relative to the reference\">\n")
    file_out.write("##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n")
    file_out.write("##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the reference\">\n")
    if not options.tandem_duplications_as_insertions  or \
       not options.interspersed_duplications_as_insertions :
        print("##ALT=<ID=DUP,Description=\"Duplication\">", file=file_out)
    if not options.tandem_duplications_as_insertions :
        print("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">", file=file_out)
    if not options.interspersed_duplications_as_insertions :
        print("##ALT=<ID=DUP:INT,Description=\"Interspersed Duplication\">", file=file_out)
    file_out.write("##ALT=<ID=INV,Description=\"Inversion of reference sequence\">\n")
    file_out.write("##ALT=<ID=BND,Description=\"Breakend of translocation\">\n")

    # INFO
    
    
    file_out.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    
    file_out.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    file_out.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n")
    file_out.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    file_out.write("##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description=\"Genomic origin of interspersed duplication seems to be deleted\">\n")
    file_out.write("##INFO=<ID=Reverse,Number=0,Type=Flag,Description=\"Cut-paste forward or reverse\">\n")
    file_out.write("##INFO=<ID=Dest_start,Number=1,Type=Integer,Description=\"Dest start on reference\">\n")
    file_out.write("##INFO=<ID=Dest_chr,Number=1,Type=String,Description=\"Dest chr name\">\n")
    print("##INFO=<ID=READS,Number=.,Type=String,Description=\"Names of all supporting reads\">", file=file_out)
    file_out.write("##INFO=<ID=STRAND,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n")
    # FORMAT Reverse
    # file.write("\n")
    file_out.write("##INFO=<ID=Support,Number=0,Type=Integer,Description=\"Support reads\">\n")
    # file_out.write("##INFO=<ID=Reverse,Number=0,Type=Flag,Description=\"Reverse or not\">\n")
    print("##FILTER=<ID=not_fully_covered,Description=\"Tandem duplication is not fully covered by a contig\">", file=file_out)
    print("##FILTER=<ID=incomplete_inversion,Description=\"Only one inversion breakpoint is supported\">", file=file_out)
    file_out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    if not options.tandem_duplications_as_insertions:
        file_out.write("##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of tandem duplication (e.g. 2 for one additional copy)\">\n")
    
    file_out.write("##CommandLine=\"cuteAsm %s\"\n"%(" ".join(argv)))
    file_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+options.sample+"\n")
    # Prepare VCF entries depending on command-line parameters
    vcf_entries = []
    for candidate in candidates:
        if candidate.type=="DEL" :
            contig, start, end = candidate.get_key()
            vcf_entries.append(((contig, max(1, start), end), candidate.get_vcf_entry(reference=reference), "DEL"))
        elif  candidate.type=="INV":
            contig, start, end = candidate.get_key()
            vcf_entries.append(((contig, start+1, end), candidate.get_vcf_entry(reference=reference), "INV"))
        elif  candidate.type=="INS":
            contig, start, end = candidate.get_source()
            vcf_entries.append(((contig, max(1, start), end), candidate.get_vcf_entry(reference=reference), "INS"))
        elif  candidate.type=='DUP_TAN':
            if options.tandem_duplications_as_insertions:
                vcf_entries.append(((candidate.source_contig, candidate.source_start+1, candidate.source_end), candidate.get_vcf_entry_as_ins(reference=reference), "INS"))
            else:
                vcf_entries.append(((candidate.source_contig, candidate.source_start+1, candidate.source_end), candidate.get_vcf_entry_as_dup(), "DUP"))
        elif candidate.type=='DUP_INT':
            # vcf_entries.append(((candidate.get_source()[0], candidate.get_source()[1]+1, candidate.get_source()[1] + 2), candidate.get_vcf_entry(reference=reference), "BND"))
            contig, start, end = candidate.get_source()
            vcf_entries.append(((contig, start+1, end), candidate.get_vcf_entry_as_dup(options.query_names), "DUP_INT"))

        else:
            vcf_entries.append(((candidate.get_source()[0], candidate.get_source()[1]+1, candidate.get_source()[1] + 2), candidate.get_vcf_entry(reference=reference), "BND"))


    # Sort and write entries to VCF
    svtype_counter = defaultdict(int)
    for source, entry, svtype in sorted_nicely(vcf_entries):
        variant_id = "cute_asm.{svtype}.{number}".format(svtype = svtype, number = svtype_counter[svtype] + 1)
        entry_with_id = entry.replace("PLACEHOLDERFORID", variant_id, 1)
        svtype_counter[svtype] += 1
        print(entry_with_id, file=file_out)

    file_out.close()
