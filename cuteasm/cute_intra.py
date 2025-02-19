from __future__ import print_function

import sys
import pysam
from cute_candidate import *
from realign import ins_filter
def generate_combine_sigs(sigs, Chr_name, read_name, svtype, candidate, merge_dis):
    if len(sigs) == 0:
        pass
    elif len(sigs) == 1:
        if sigs[0][1]>30:

            if svtype == 'INS':
                candidate.append([1,Chr_name, sigs[0][0], sigs[0][1],read_name,sigs[0][2]])
                
            else:
                candidate.append([0,Chr_name, sigs[0][0], sigs[0][1],read_name])
    else:
        temp_sig = sigs[0]
        if svtype == "INS":
            temp_sig += [sigs[0][0]]
            for i in sigs[1:]:
                if i[0] - temp_sig[3] <= merge_dis:
                    temp_sig[1] += i[1]
                    temp_sig[2] += i[2]
                    temp_sig[3] = i[0]
                else:
                    if temp_sig[1]>30:
                        candidate.append([1,Chr_name, temp_sig[0], temp_sig[1],read_name, temp_sig[2]])
                    temp_sig = i
                    temp_sig.append(i[0])
            if temp_sig[1]>30:
                candidate.append([1,Chr_name, temp_sig[0], temp_sig[1],read_name, temp_sig[2]])
        else:
            temp_sig += [sum(sigs[0])]
            # merge_dis_bias = max([i[1]] for i in sigs)
            for i in sigs[1:]:
                if i[0] - temp_sig[2] <= merge_dis:
                    temp_sig[1] += i[1]
                    temp_sig[2] = sum(i)
                    print('del')
                else: 
                    if temp_sig[1]>30:
                        candidate.append([0,Chr_name, temp_sig[0], temp_sig[1], read_name])
                    temp_sig = i
                    temp_sig.append(i[0])
            if temp_sig[1]>30:
                candidate.append([0,Chr_name, temp_sig[0], temp_sig[1], read_name])
OPLIST=[
    pysam.CBACK,
    pysam.CDEL,
    pysam.CDIFF,
    pysam.CEQUAL,
    pysam.CHARD_CLIP,
    pysam.CINS,
    pysam.CMATCH,
    pysam.CPAD,
    pysam.CREF_SKIP,
    pysam.CSOFT_CLIP
]
RefChangeOp=set([0,2,7,8])

#QUERY CHANGE, REF CHANGE
CHANGETABLE={
    pysam.CMATCH:     (True,True),
    pysam.CINS:       (True,False),
    pysam.CDEL:       (False,True),
    pysam.CREF_SKIP:  (False,True),
    pysam.CPAD:       (False,False),
    pysam.CEQUAL:     (True,True),
    pysam.CDIFF:      (True,True)
}
CHANGEOP=[CHANGETABLE[i] if i in CHANGETABLE.keys() else (False,False) for i in range(max(OPLIST)+1)]
REFCHANGEOP=[CHANGETABLE[i][1] if i in CHANGETABLE.keys() else False for i in range(max(OPLIST)+1)]
INDELOP=[(i==pysam.CDEL or i==pysam.CINS) for i in range(max(OPLIST)+1)]
def is_simalar(a,b):
    if a[0]==b[0] and a[1]==b[1] and a[2]==b[2]:
        return 1
    return 0
def is_simalar2(a,b):
    if a[0]==b[0] and a[1]==b[1] :
        return 1#500,10,0,100
def parse_read(read, sv_candidates,Chr_name,ref, min_read_len=500, min_siglength=10, merge_del_threshold=10, merge_ins_threshold=10):
    del_candidate=sv_candidates['DEL']
    ins_candidate=[]
    if read.query_length < min_read_len:
        return []
    Combine_sig_in_same_read_ins = list()
    Combine_sig_in_same_read_del = list()

    #new start
   
    pos_start = read.reference_start # 0-based
    pos_end = read.reference_end
    sig_start=pos_start
    softclip_left = 0
    softclip_right = 0
    hardclip_left = 0
    hardclip_right = 0
    shift_ins_read = 0
    if read.cigar[0][0] == 4:
        softclip_left = read.cigar[0][1]
    elif read.cigar[0][0] == 5:
        hardclip_left = read.cigar[0][1]
    
    shift_ins_read=-hardclip_left
    for op, oplen in read.cigartuples:
        # calculate offset of an ins sig in read
        if op != 2:#might be fixed later
            shift_ins_read += oplen
        if oplen >= min_siglength and INDELOP[op]:
            if op==2:
                Combine_sig_in_same_read_del.append([sig_start, oplen])
                sig_start += oplen
            else:
                Combine_sig_in_same_read_ins.append([sig_start, oplen,
                    str(read.query_sequence[shift_ins_read-oplen:shift_ins_read])])
        else:
            # if op in RefChangeOp:
            if REFCHANGEOP[op]:
                sig_start += oplen

    
    if read.cigar[-1][0] == 4:
        softclip_right = read.cigar[-1][1]
    elif read.cigar[-1][0] == 5:
        hardclip_right = read.cigar[-1][1]

    if hardclip_left != 0:
        softclip_left = hardclip_left
    if hardclip_right != 0:
        softclip_right = hardclip_right

    # ************Combine signals in same read********************
    generate_combine_sigs(Combine_sig_in_same_read_ins, Chr_name, read.query_name, "INS", ins_candidate, merge_ins_threshold)
    
    generate_combine_sigs(Combine_sig_in_same_read_del, Chr_name, read.query_name, "DEL", del_candidate, merge_del_threshold)
    # dup_list,ins_list=ins_filter(ins_candidate,ref)
    # sv_candidates['DUP_TAN'].extend(dup_list)
    sv_candidates['INS'].extend(ins_candidate)

def analyze_cigar_indel(tuples, min_length):
    """Parses CIGAR tuples (op, len) and returns Indels with a length > minLength"""
    pos_ref = 0
    pos_read = 0
    indels = []
    #M I D N S H P = X
    #0 1 2 3 4 5 6 7 8
    #pos_ref 和pos_read 分别标记变异出现相对于起始的位置
    for operation, length in tuples:
        if operation == 0:                     # alignment match
            pos_ref += length
            pos_read += length
        elif operation == 1:                   # insertion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "INS"))
            pos_read += length
        elif operation == 2:                   # deletion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "DEL"))
            pos_ref += length
        elif operation == 4:                   # soft clip
            pos_read += length
        elif operation == 7 or operation == 8:        # match or mismatch
            pos_ref += length
            pos_read += length
    return indels

def analyze_alignment_indel(alignment, bam, query_name, options,sv_candidates):
    #获得alignment的参考基因组的名称及起始位点
    ref_chr = bam.getrname(alignment.reference_id)
    ref_start = alignment.reference_start
    del_candidate=sv_candidates['DEL']
    ins_candidate=sv_candidates['INS']
    #从cigar中获取indel信息
    indels = analyze_cigar_indel(alignment.cigartuples, options.min_sv_size)
    for pos_ref, pos_read, length, typ in indels:
        if typ == "DEL":#！！！为什么都要加length
            del_candidate.append([0,ref_chr, ref_start + pos_ref, length, query_name])
            # sv_candidates[typ].append([ref_chr, ref_start + pos_ref, length, query_name])
        elif typ == "INS":
            insertion_seq = alignment.query_sequence[pos_read:pos_read+length]
            ins_candidate.append([1,ref_chr, ref_start + pos_ref, length, query_name, insertion_seq])
            # sv_candidates[typ].append([ref_chr, ref_start + pos_ref, length, query_name, insertion_seq])#是不是有点问题 ins ref位点为ref--ref+length
    # sv_candidates['INS'].extend(cluster_ins(ins_candidate))
    # sv_candidates['DEL'].extend(cluster_del(del_candidate))
    # return sv_candidates


