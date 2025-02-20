import sys
import os
import logging
import argparse

VERSION = '0.1.1'

class cuteAsmdp(object):
	'''
	Detailed descriptions of cuteSV version and its parameters.
	'''

	USAGE="""\
		
	Current version: v%s
	Author: Tao Jiang
	Contact: tjiang@hit.edu.cn

	If you use cuteSV in your work, please cite:
		Jiang T et al. Long-read-based human genomic structural variation detection with cuteSV. 
		Genome Biol 21,189(2020). https://doi.org/10.1186/s13059-020-02107-y


	Suggestions:

	For PacBio CLR data:
		--max_cluster_bias_INS		100
		--diff_ratio_merging_INS	0.3
		--max_cluster_bias_DEL	200
		--diff_ratio_merging_DEL	0.5

	For PacBio CCS(HIFI) data:
		--max_cluster_bias_INS		1000
		--diff_ratio_merging_INS	0.9
		--max_cluster_bias_DEL	1000
		--diff_ratio_merging_DEL	0.5

	For ONT data:
		--max_cluster_bias_INS		100
		--diff_ratio_merging_INS	0.3
		--max_cluster_bias_DEL	100
		--diff_ratio_merging_DEL	0.3


	"""%(VERSION)

	# MinSizeDel = 'For current version of cuteSV, it can detect deletions larger than this size.'

class Options:
    # def __init__(self, args):
    #     # 基本参数
    #     self.input = args.input
    #     self.reference = args.reference
    #     self.refbam=args.refbam
    #     self.output = args.output
    #     self.work_dir = args.work_dir
    #     self.debug = args.debug
    #     self.threads = args.threads
    #     self.sample = args.sample
    #     self.retain_work_dir = args.retain_work_dir
    #     self.report_readid = args.report_readid
    #     self.bam1 = args.input[0]
    #     self.bam2 = args.input[1] if len(args.input) > 1 else None
    #     # SV签名参数
    #     self.max_split_parts = args.max_split_parts
    #     self.min_mapq = args.min_mapq
    #     self.min_read_len = args.min_read_len
    #     self.include_bed = args.include_bed
    #     self.min_sv_size = args.min_sv_size
    #     self.max_sv_size = args.max_sv_size
    #     self.max_edit_distance=args.max_edit_distance
    #     self.partition_max_distance=args.partition_max_distance

    #     # 基因型参数
    #     self.genotype = args.genotype

    #     # 输出参数
    #     self.symbolic_alleles = args.symbolic_alleles
    #     self.tandem_duplications_as_insertions = args.tandem_duplications_as_insertions
    #     self.interspersed_duplications_as_insertions=args.interspersed_duplications_as_insertions
    #     self.query_names = args.query_names
    # def __init__(self):
        
    #     # 基本参数
    #     # self.input = ['/io/wuxiaomia/Titan/data/aln/Assembly/NA24385_shiJie/ch38_no_samll_hp1.bam ','/io/wuxiaomia/Titan/data/aln/Assembly/NA24385_shiJie/ch38_no_samll_hp2.bam']
    #     self.reference = '/io/wuxiaomia/Titan/data/GRCh38_full_analysis_set_plus_decoy_hla.fa'#/io/wuxiaomia/Titan/data/GRCh38_chromosomes.fa'
    #     self.output = 'add_bnd.vcf'
    #     self.work_dir ='/io/wuxiaomia/Titan/data/aln/Assembly/simulations/Assembly_shijie/cuteasm/new1'
    #     self.debug = False
    #     self.threads = 1
    #     self.refbam=None#'/io/wuxiaomia/Titan/data/GRCh38_chr1_chr2.bam'##'/io/wuxiaomia/Titan/data/GRCh38_full_analysis_set_plus_decoy_hla.bam'##'/io/wuxiaomia/Titan/data/aln/Assembly/simulations/visor/small.reference.bam'
    #     self.sample = 'sample'
    #     # self.retain_work_dir = args.retain_work_dir
    #     # self.report_readid = args.report_readid
    #     self.bam1 = '/io/wuxiaomia/Titan/data/aln/Assembly/simulations/Assembly_shijie/h1.sorted.bam'
    #     self.bam2 = '/io/wuxiaomia/Titan/data/aln/Assembly/simulations/Assembly_shijie/h2.sorted.bam'#None#
    #     # SV签名参数
    #     # self.max_split_parts = args.max_split_parts
    #     self.min_mapq = 20
    #     # self.min_read_len = args.min_read_len
    #     # self.include_bed = args.include_bed
    #     # self.min_sv_size = 40
    #     self.max_sv_size = 100000
    #     self.max_edit_distance=10
    #     self.partition_max_distance=1000

    #     # 基因型参数
    #     self.genotype = True

    #     # 输出参数
    #     # self.symbolic_alleles = args.symbolic_alleles
    #     self.tandem_duplications_as_insertions = False
    #     self.interspersed_duplications_as_insertions=False
    #     self.query_names = True
    def __init__(self, args=None):
        if args:
            # 基本参数
            self.input = args.input
            self.reference = args.reference
            self.refbam = args.refbam
            self.output = args.output
            self.work_dir = args.work_dir
            self.debug = args.debug
            self.threads = args.threads
            self.sample = args.sample
            self.retain_work_dir = args.retain_work_dir
            self.report_readid = args.report_readid
            self.bam1 = args.input[0]
            self.bam2 = args.input[1] if len(args.input) > 1 else None
            # SV签名参数
            self.max_split_parts = args.max_split_parts
            self.min_mapq = args.min_mapq
            self.min_read_len = args.min_read_len
            self.include_bed = args.include_bed
            self.min_sv_size = args.min_sv_size
            self.max_sv_size = args.max_sv_size
            self.max_edit_distance = args.max_edit_distance
            self.partition_max_distance = args.partition_max_distance

            # 基因型参数
            self.genotype = args.genotype

            # 输出参数
            self.symbolic_alleles = args.symbolic_alleles
            self.tandem_duplications_as_insertions = args.tandem_duplications_as_insertions
            self.interspersed_duplications_as_insertions = args.interspersed_duplications_as_insertions
            self.query_names = args.query_names
        else:#自定义参数
            # 基本参数
            self.input = ['/io/wuxiaomia/Titan/data/aln/Assembly/NA24385_shiJie/ch38_no_samll_hp1.bam ',
                          '/io/wuxiaomia/Titan/data/aln/Assembly/NA24385_shiJie/ch38_no_samll_hp2.bam']
            self.reference = '/io/wuxiaomia/Titan/data/GRCh38_chromosomes.fa'
            self.output = 'test_.vcf'
            self.work_dir = '/io/wuxiaomia/Titan/data/aln/Assembly/NA24385_shiJie/cuteasm/new'
            self.debug = False
            self.threads = 1
            self.refbam = None
            self.sample = 'sample'
            self.bam1 = '/io/wuxiaomia/Titan/data/aln/Assembly/NA24385_shiJie/ch38_no_samll_hp1.bam'
            self.bam2 = '/io/wuxiaomia/Titan/data/aln/Assembly/NA24385_shiJie/ch38_no_samll_hp2.bam'
            # SV签名参数
            self.min_mapq = 20
            self.max_sv_size = 100000
            self.max_edit_distance = 10
            self.partition_max_distance = 1000

            # 基因型参数
            self.genotype = True

            # 输出参数
            self.tandem_duplications_as_insertions = False
            self.interspersed_duplications_as_insertions = False
            self.query_names = True


def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="cuteAsm", 
                                     description=cuteAsmdp.USAGE, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--version', '-v', 
                        action='version', 
                        version='%(prog)s {version}'.format(version=VERSION))

    # **************Parameters of input******************
    parser.add_argument('-b', '--input', type=str, nargs='+', required=True,
                        help="Input BAM file(s) (space-separated list for multiple files)")
    parser.add_argument('-r',"--reference",  required=True,
                        type=str, 
                        help="The reference genome in fasta format.")
    parser.add_argument('-o','--output' ,
                        type=str, 
                        default='varints.vcf',
                        help="Output VCF format file.")
    parser.add_argument('-w','--work_dir',  required=True,
                        type=os.path.abspath, 
                        help="Work-directory for distributed jobs")
    
    parser.add_argument('--refbam',  required=False,default=None,help="Reference bam file.")
    
    # Verbosity flag
    parser.add_argument('--debug', 
                        action='store_true', 
                        help='Enable more verbose logging (default: %(default)s)')

    # ************** Other Parameters******************
    parser.add_argument('-t', '--threads', 
                        help="Number of threads to use. [%(default)s]", 
                        default=16, 
                        type=int)
    parser.add_argument('--retain_work_dir', 
                        help="Enable to retain temporary folder and files.", 
                        action="store_true")
    parser.add_argument('--report_readid', 
                        help="Enable to report supporting read ids for each SV.", 
                        action="store_true")

    # ************** Parameters in signatures collection ******************
    GroupSignaturesCollect = parser.add_argument_group('Collection of SV signatures')
    GroupSignaturesCollect.add_argument('-p', '--max_split_parts', 
                                        help="Maximum number of split segments a read may be aligned before it is ignored. [%(default)s]", 
                                        default=10, 
                                        type=int)
    GroupSignaturesCollect.add_argument('-q', '--min_mapq', 
                                        help="Minimum mapping quality value of alignment to be considered. [%(default)s]", 
                                        default=20, 
                                        type=int)
    GroupSignaturesCollect.add_argument('-l', '--min_read_len', 
                                        help="Ignores reads that report alignments no longer than this length (default: %(default)s)", 
                                        default=500, 
                                        type=int)
    GroupSignaturesCollect.add_argument('-m', '--max_edit_distance', 
                                        help="edit (default: %(default)s)", 
                                        default=10, 
                                        type=int)
    GroupSignaturesCollect.add_argument('-d', '--partition_max_distance', 
                                        help="edit (default: %(default)s)", 
                                        default=1000, 
                                        type=int)
    GroupSignaturesCollect.add_argument('-include_bed', 
                                        help="Optional given bed file. Only detect SVs in regions in the BED file. [NULL]", 
                                        default=None, 
                                        type=str)
    GroupSignaturesCollect.add_argument('--min_sv_size',
                                        type=int,
                                        default=30,
                                        help="Minimum SV size to detect (default: %(default)s).")
    GroupSignaturesCollect.add_argument('--max_sv_size',
                                        type=int,
                                        default=100000,
                                        help="Maximum SV size to detect (default: %(default)s).")

    # ************** Parameters in genotyping ******************
    GroupGenotype = parser.add_argument_group('Computing genotypes')
    GroupGenotype.add_argument('--genotype',
                               help="Enable to generate genotypes.",
                               action="store_true")

    # ************** Parameters in Output ******************
    group_output = parser.add_argument_group('OUTPUT')
    group_output.add_argument('--sample', 
                              type=str, 
                              default="Sample", 
                              help='Sample ID to include in output vcf file (default: %(default)s)')
    group_output.add_argument('--symbolic_alleles', 
                              action='store_true', 
                              help="Use symbolic alleles, such as <DEL> or <INV> in the VCF output.")
    group_output.add_argument('--tandem_duplications_as_insertions', 
                              action='store_true', 
                              help="Represent tandem duplications as insertions in output VCF.")
    group_output.add_argument('--interspersed_duplications_as_insertions', 
                              action='store_true', 
                              help="Represent interspersed duplications as insertions in output VCF.")
    group_output.add_argument('--query_names', 
                              action='store_true', 
                              help="Output names of supporting query sequences in INFO tag of VCF.")
    
    args = parser.parse_args(argv)
    return args