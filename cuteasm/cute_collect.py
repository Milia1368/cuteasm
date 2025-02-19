import logging
import pysam

from cute_intra import parse_read
from cute_inter import *
# from cute_filter import *
# from cluster import cluster_dupint


def retrieve_other_alignments(main_alignment, bam):
    """Reconstruct other alignments of the same read for a given alignment from the SA tag"""
    #reconstructing other alignments from SA tag does not work if sequence of main_alignment is hard-clipped
    if main_alignment.get_cigar_stats()[0][5] > 0:#判断cigar统计符号数中，H是否存在；如果有硬剪裁则直接跳过 ！！！为什么硬剪裁直接跳过
        return []
    try:
        sa_tag = main_alignment.get_tag("SA").split(";")#SA:Z:chr2,50,+，9M,30,0;
        #SA:Z:chr12,22924,-,19101S56105M4131I43910S,60,9493;chr1,108374,-,79361S33053M38678D10833S,60,38967;
    except KeyError:
        return []
    other_alignments = []
    # For each other alignment encoded in the SA tag
    for element in sa_tag:#对SA每一条记录进行提取
        # Read information from the tag
        fields = element.split(",")
        if len(fields) != 6:
            continue
        rname = fields[0]
        pos = int(fields[1])
        strand = fields[2]
        # CIGAR string encoded in SA tag is shortened
        cigar = fields[3]
        mapq = int(fields[4])
        nm = int(fields[5])#编辑距离

        # Generate an aligned segment from the information
        #提取SA中信息转换为一条记录
        a = pysam.AlignedSegment()
        a.query_name = main_alignment.query_name#和最初的alignment的编号一样，但是其他位置设置自己的信息
        a.query_sequence = ''
        if strand == "+":#
            a.flag = 2048
        else:
            a.flag = 2064
        a.reference_id = bam.get_tid(rname)
        a.reference_start = pos - 1
        try:
            a.mapping_quality = mapq
        except OverflowError:
            a.mapping_quality = 0
        try:
            a.cigarstring = cigar
        except OverflowError:
            logging.error("OverflowError while retrieving supplementary CIGAR string. Read name: {0}, Position: {1}, CIGAR: {2}".format(rname, pos, cigar))
            continue
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
        a.query_qualities = ''
        a.set_tags([("NM", nm, "i")])

        other_alignments.append(a)
    return other_alignments
import time
from overlap_collect import collect,plot_deta,plot_d
from cute_inter import aln_filter
from realign import paste_del_match,cluster_dupint,intra_ins_del_match
def startic(primary, supplementaries,deta_list,d_list,header):
    query=primary.query_sequence
    query_length=primary.infer_read_length()
    alignments = [primary] + supplementaries
    alignment_list = []
    for alignment in alignments:#
        try:
            refename=alignment.reference_name
        except:
            refename= header['SQ'][alignment.reference_id]['SN']
        alignment_list.append([refename, alignment.query_alignment_start,alignment.query_alignment_end,alignment.reference_start,alignment.reference_end,alignment.is_reverse])
    #filter在内部的alignment等
    alignment_list=aln_filter(alignment_list)
    # alignment_list.sort(key=lambda x:(x[1],x[2]))
    read_name=primary.query_name
    collect(alignment_list,deta_list,read_name,d_list)
def analyze_alignment_file_coordsorted(bam,reference,current_chromosome,options,inter_sv_candidates,intra_candidates):#options 是命令行参数
    header = bam.header
    # deta_list=[]
    # d_list=[]
    # print("Processing chromosome {0}...".format(current_chromosome))
    alignment_it = bam.fetch(contig = current_chromosome)#找到bam中对应的记录
    paste=[]
    dup_tan=[]
    dup_int=[]
    cut_paste=[]
    bnd=[]
    merged_inter = {svtype: [] for svtype in types_to_output}
    merged_intra = {svtype: [] for svtype in ['DEL', 'INS']}
    # if 'chr1' !=current_chromosome:
    #     continue
    start_time=time.time()
    num=0
    
    while True:
        try:
            num+=1
            current_alignment = next(alignment_it)#迭代对象 对于每一条记录 ！！！是指每一条记录吗
            # if current_alignment.query_name!='h1tg000090l':#'h2tg000056l':#'h2tg000008l':#'h2tg000022l' 'h2tg000076l'!!!
            #     continue
            if current_alignment.is_unmapped or current_alignment.is_secondary or current_alignment.mapping_quality < options.min_mapq:
                continue #过滤掉没比对上0X4\secondary0X100（256\比对质量低的 直接不分析这些reads0X800(2048)
                #secondary 主比对之外的比对 次级比对，可能会在判定复杂区域或者评估比对软件质量有用
            if current_alignment.is_supplementary:#补充比对 用于表示跨度较大的比对
                #！#因为不是主比对，次要比对，所以在没有主比对信息之前是没有用的，所以只关注它本身存在的变异indel
                parse_read(current_alignment,merged_intra,current_alignment.reference_name,reference)
                # pass
            else:#其他比对 其他比对指的是什么
                #！#这里包含主比对（包含SA标签），和完全好的比对
                supplementary_alignments = retrieve_other_alignments(current_alignment, bam)#检索其他比对信息：对该align的SA进行进一步提取
                good_suppl_alns = [aln for aln in supplementary_alignments if not aln.is_unmapped and aln.mapping_quality >= options.min_mapq]#对质量好的进行保留
                #保留质量好的比对 然后进行分析
                parse_read(current_alignment,merged_intra,current_alignment.reference_name,reference)
                a,b,c=analysis_split_read(current_alignment, good_suppl_alns,merged_inter,current_alignment.query_name,header,options)
                paste.extend(a)
                dup_tan.extend(b)
                bnd.extend(c)
                
                # startic(current_alignment, good_suppl_alns,deta_list,d_list,header)
        except StopIteration:
            break
    # dup_int_clusters, trans_clusters=cluster_dupint(dup_int)
    #对paster匹配
    del_list=[]
    del_list.extend(merged_inter['DEL']+merged_intra['DEL'])
    del_list=paste_del_match(paste,del_list,reference,dup_tan,dup_int,cut_paste)

    #对每个染色体内部的intra和inter进行匹配
    ins_list=[]
    ins_list.extend(merged_inter['INS']+merged_intra['INS'])
    intra_candidates['INS'],intra_candidates['DEL'],cutpastein=intra_ins_del_match(ins_list,del_list,reference)
    cut_paste.extend(cutpastein)
    # dup_int=cluster_dupint(dup_int)
    dup_tan.extend(merged_inter[types_to_output[duptan_num]])
    #bnd
    
    
    inter_sv_candidates[types_to_output[duptan_num]].extend(cluster_dup(dup_tan))
    inter_sv_candidates[types_to_output[dupint_num]].extend(dup_int+cut_paste)
    inter_sv_candidates[types_to_output[ins_num]]=[]
    inter_sv_candidates[types_to_output[del_num]]=[]

    for type in [types_to_output[bnd_num],types_to_output[inv_num]]:
        inter_sv_candidates[type].extend(merged_inter[type])
    end_time=time.time()
    logging.info("Processing chromosome {0}...".format(current_chromosome))
    logging.info("num {0}...".format(num))
    logging.info("time:{}".format(end_time-start_time))
    # plot_deta(deta_list)export
    # plot_d(d_list)
    # with open("deta_overlap.txt",'w') as f:
    #     for line in deta_list:

    #         print(line,file=f)
    return bnd