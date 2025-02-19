"""仿照forcecalling将ins区域在区间范围内重比对"""

from edlib import align
from collections import Counter
# #输入 insertion
#输出 将一部分insertion重比对作为dup,如果比对落在其信号附近
def dup_find_nearby(query,reference,dup_list):
    
    length=query[3]
    readname=query[4]
    queryseq=query[5]
    start=query[2]
    chr=query[1]
    ref_length=chr_length = reference.get_reference_length(chr)
    flanking=min(int(length*1.5),25000)
    region_start=max(0,start-flanking)
    targetseq=reference.fetch(chr, region_start, min(start+flanking,ref_length)).upper()
    align_ans=align(queryseq, targetseq, mode="HW", task="locations", k=0)

    edi_dis=align_ans["editDistance"]
    if edi_dis!=-1:#说明可以比对上
        edii_lo=align_ans["locations"]
        ed_cigar=align_ans["cigar"]
        new_pos_start=region_start+edii_lo[0][0]
        new_pos_end=region_start+edii_lo[0][1]
        dup_list.append([chr, new_pos_start, new_pos_end, 1, True, readname])
        return True
    # 否则认为新插入序列
    return False
def ins_filter(ins_list,reference):
    dup_list=[]
    ins_final=[]
    for ins_can in ins_list:
       if not  dup_find_nearby(ins_can,reference,dup_list):
           ins_final.append(ins_can)
    return dup_list,ins_final

#ins 和 del 匹配 以查找潜在的cut-paste tra
#ins[]=[1,chr, pos_start, length, readname,queryseq]
#del[]=[0,chr, pos_start,  length, readname]
#匹配的信号满足，1

#对ins建立数据结构 排序后的inss 以50bplenth为桶，位点排序后的（不排序也行）的inss
def ins_bucket(sorted_ins_list):
    # ins_list.sort(key=lambda x: (x[1], x[2]))
    num=0
    ins_bucket={}
    for ins_can in sorted_ins_list:
        new_ins=[num,ins_can[1],ins_can[2],ins_can[3]]
        length=ins_can[3]
        bucket=length//50
        if bucket not in ins_bucket:
            ins_bucket[bucket]=[]
        ins_bucket[bucket].append(new_ins)
        num+=1
    return ins_bucket
#查找在某一length范围内的某一区域的ins
def find_ins_in_bucket(ins_bucket, chr, pos_start,pos_end, length):
    bucket = length // 50
    result = []
    for offset in [-1, 0, 1]:
        target_bucket = bucket + offset
        if target_bucket in ins_bucket:
            ins_list = ins_bucket[target_bucket]
            for ins_can in ins_list:
                if ins_can[1] == chr and pos_end+200 >= ins_can[2] >= pos_start-200:
                    result.append(ins_can)
    return result
#del buget
def del_bucket(sorted_del_list):
    # ins_list.sort(key=lambda x: (x[1], x[2]))
    num=0
    del_bucket={}
    for del_can in sorted_del_list:
        new_del=[num,del_can[1],del_can[2],del_can[3]]
        length=del_can[3]
        bucket=length//50
        if bucket not in del_bucket:
            del_bucket[bucket]=[]
        del_bucket[bucket].append(new_del)
        num+=1
    return del_bucket
#查找在某一length范围内的某一区域的del
def find_del_in_bucket(del_bucket, chr, pos_start,pos_end, length):
    bucket = length // 50
    result = []
    try:
        for offset in [-1, 0, 1]:
            target_bucket = bucket + offset
            if target_bucket in del_bucket:
                del_list = del_bucket[target_bucket]
                for del_can in del_list:
                    if del_can[1] == chr and pos_end+200 >= del_can[2] >= pos_start-200:#只匹配相同染色体上的cut_pasete
                        result.append(del_can)
        return result
    except:
        # 如果出现异常，返回一个默认值，防止del_result未赋值
        return []
#查找在长度范围内的del
def find_del_in_bucket_forintra(del_bucket, chr, length):
    bucket = length // 50
    result = []
    try:
        for offset in [-1, 0, 1]:
            target_bucket = bucket + offset
            if target_bucket in del_bucket:
                del_list = del_bucket[target_bucket]
                for del_can in del_list:
                    if del_can[1] == chr :#只匹配相同染色体上的cut_pasete
                        result.append(del_can)
        return result
    except:
        # 如果出现异常，返回一个默认值，防止del_result未赋值
        return []
#匹配
from  Bio.Seq import Seq
def simi_k(queryseq,del_can,reference):
    flag=False
    chrom=del_can[1]
    region_start=del_can[2]
    region_end=del_can[3]+region_start
    targetseq=reference.fetch(chrom, region_start, region_end).upper()
    align_ans=align(queryseq, targetseq, mode="NW", task="distance",k=20)
    edi_dis=align_ans["editDistance"]
    if edi_dis==-1:
        query_rc = str(Seq(queryseq).reverse_complement())
        align_ans1=align(query_rc, targetseq, mode="SHW", task="distance",k=20)
        if align_ans1["editDistance"]==-1:
            return (-1,flag)
        edi_dis=align_ans1["editDistance"]
        flag=True
    if edi_dis==-1:
        edi_dis=100
   
    return (edi_dis,flag)

# [1,alignment_list[index][al_chr_name],flag,alignment_list[index][al_r_start], alignment_list[index][al_r_end],ele[al_chr_name],ele[al_is_reverse],ele[al_r_start],ele[al_r_end],detlength,readname]
#[duptan_num,ele2[al_chr_name],ele2[al_r_start],ele2[al_r_end],False,ele2[al_is_reverse],ele2[al_q_start],readname,3]
#([bnd_num,ele1[al_chr_name],ele1[al_r_start],ele2[al_chr_name],ele2[al_r_end],readname,'Nor']
# [del_num,ele1[al_chr_name],ele1[al_r_end],-delta_length,readname]
# def paste_del_match(paste,del_list,reference,dup_tan,dup_int,rm_cut):
#     # dup_tan=[]
#     # dup_int=[]
#     cut_paste=[]
#     if len(paste)==0:
#         return del_list
#     del_list.sort(key=lambda x: (x[1], x[2]))
#     del_buget=del_bucket(del_list)
#     del_delindex=[]
#     index=0
#     for can in paste:
        
#         flag=False
#         readname=can[-1]
#         detalen=can[-2]
#         if can[0]==1:
#             chrom=can[1]
#             ins_start=can[3]
#             ins_end=can[4]
#             last_start=can[7]
#             last_end=can[8]
#             ins_flag=can[2]
#             last_flag=can[6]
#         else:
#             chrom=can[5]
#             ins_start=can[7]
#             ins_end=can[8]
#             last_start=can[3]
#             last_end=can[4]
#             ins_flag=can[6]
#             last_flag=can[2]
#         del_result=find_del_in_bucket(del_buget,chrom,ins_start,ins_end,ins_end-ins_start)
#         if len(del_result)<5:
#             queryseq=reference.fetch(chrom, ins_start, ins_end).upper()
#             k_list = []
#             bool_list = []
#             for del_can in del_result:
#                 k, is_true = simi_k(queryseq, del_can, reference)
#                 k_list.append(k)
#                 bool_list.append(is_true)
#             if len(k_list)>0:
#                 min_k = min(k_list)
#                 if min_k < 20 :
#                     for index, k in enumerate(k_list):
#                         if k == min_k:
#                             del_can = del_result[index]
#                             del_index = del_can[0]
#                             is_true = bool_list[index]
#                             del_delindex.append(del_index)
#                                             # [-1,del_can[1],del_can[2],del_can[3]+del_can[2],chrom,ins_start,ins_end,-1,readname,True]
#                                             # [-1,del_can[1],del_can[2],del_can[3]+del_can[2],chrom,ins_start,ins_end,ins_flag,readname,True]
#                             cut_paste.append([del_index,k,del_can[1], index,del_can[2], del_can[3]+del_can[2], chrom, last_end+detalen,last_end+detalen+ins_end-ins_start, is_true, readname, True])
#                             flag = True
#                             break
                
#         if not flag:
#             print("Not find pair del,indicates a paste")
#             #处理是dup tan还是dup end  还有tra?
#             if ins_end-last_end>-10000:
#                 #默认为dup tan
#                 dup_tan.append([-1,chrom,ins_start,ins_end,True,ins_flag,-1,readname,'delpair'])
#             else:
#                 #默认为dup int
#                 qs=last_end+detalen
#                 length=ins_end-ins_start
#                 dup_int.append([-1,chrom,ins_start,ins_end,ins_flag,qs,length,readname,1,'pari',False])
#         index+=1
#     # 使用Counter统计频次
#     counter = Counter(del_delindex)

#     # 筛选出出现次数大于1的元素
#     rm={}
#     bnd=[]
#     result = [num for num in counter if counter[num] > 1]
#     for cut in cut_paste:
#         del_index=cut[0]
#         k=cut[1]
#         if del_index in result:
#             if del_index not in rm:
#                 rm[del_index]=[]
#             rm[del_index].append(cut)
#         else:
#             rm_cut.append([-1]+cut[3:])
#     flag=False
#     for key,value in rm:
#         min_k=min([cut[1] for cut in value])
#         for i  in value:
#             if i[1]==min_k and not flag:
#                 rm_cut.append([-1]+cut[3:])
#                 flag=True

#     new_del_list=[del_list[index] for index in range(len(del_list)) if index not in del_delindex]
    
#     return new_del_list


from collections import Counter


def paste_del_match(paste, del_list, reference, dup_tan, dup_int, rm_cut):
    if not paste:
        return del_list

    del_list.sort(key=lambda x: (x[1], x[2]))
    del_budget = del_bucket(del_list)

    del_delindex = []
    cut_paste = []
    index=0
    for candidate in paste:
        flag = False
        readname = candidate[-1]
        detalen = candidate[-2]
        chrom, start, end, last_start, last_end, ins_flag, last_flag = parse_candidate(candidate)

        del_result = find_del_in_bucket(del_budget, chrom, start, end, end - start)
        if len(del_result) < 5:
            query_seq = reference.fetch(chrom, start, end).upper()
            k_and_bool = [(simi_k(query_seq, del_can, reference)) for del_can in del_result]
            if k_and_bool:
                k_values, bool_values = zip(*k_and_bool)
                min_k = min(k_values)
                if min_k < 20:
                    min_k_index = k_values.index(min_k)
                    del_can = del_result[min_k_index]
                    del_index = del_can[0]
                    is_true = bool_values[min_k_index]
                    del_delindex.append(del_index)
                    cut_paste.append([del_index, min_k, index,del_can[1],  del_can[2], del_can[3] + del_can[2], chrom,
                                      last_end + detalen, last_end + detalen + end - start, is_true, readname, True])
                    flag = True

        if not flag:
            handle_dup(last_end, dup_tan, dup_int, chrom, start, end, detalen,ins_flag, readname)
        index+=1

    repeated_indices = {num for num in Counter(del_delindex).keys() if Counter(del_delindex)[num] > 1}
    for cut in cut_paste:
        del_index = cut[0]
        if del_index not in repeated_indices:
            rm_cut.append([-1] + cut[3:])

    new_del_list = [dels for ind, dels in enumerate(del_list) if ind not in del_delindex]
    return new_del_list
    # return del_list


def parse_candidate(candidate):
    if candidate[0] == 1:
        return candidate[1], candidate[3], candidate[4], candidate[7], candidate[8], candidate[2], candidate[6]
    return candidate[5], candidate[7], candidate[8], candidate[3], candidate[4], candidate[6], candidate[2]

# from cute_inter import find_duptan

def handle_dup(last_end, dup_tan, dup_int, chrom, start, end, detalen,ins_flag, readname):
    # flag=find_duptan(prev_aln[:], current_aln[:],  query_length, readname,duptans, dupinters)
    if end-last_end < -10000:
    #     # dup_tan.append([-1, chrom, start, end, True, ins_flag, -1, readname, 'delpair'])
    # else:
        qs = last_end+detalen
        length = end - start
        dup_int.append([-1, chrom, start, end, ins_flag, qs, length, readname, 1, 'pari', False])
    pass 
#对dupintcluster
def cluster_dupint(dup_int):
    if len(dup_int)==0:
        return dup_int
    dup_int.sort(key=lambda x:(x[1],x[2]))
    result = []
    for item in dup_int:
        merged = False
        for res in result:
            if item[1] == res[1] and item[4]==res[4]:
                if abs(item[5]-res[5])<=200:
                    if max(abs(item[2]-res[2]),abs(item[3]-res[3]))  <= 50 :
                        res[2] = min(item[2], res[2])
                        res[3] = max(item[3], res[3])
                        res[5]=int((item[5]+res[5])/2)
                        res[6]=res[3]-res[2]
                        res[-3] += 1
                        merged = True
                        break
        if not merged:
            result.append(item)

    # print('clusetrdupint',result)
    return result
#对cut paste cluster
#[-1, del_can[1], del_can[2], del_can[3]+del_can[2], chrom, last_end+detalen,last_end+detalen+ins_end-ins_start, is_true, readname, True]
def cluster_cut_paste(cut_paste):
    if len(cut_paste)==0:
        return cut_paste
    cut_paste.sort(key=lambda x:(x[1],x[2],x[5]))
    result = []
    for item in cut_paste:
        merged = False
        for res in result:
            if item[4] == res[4] and item[1]==res[1]:#chr
                if item[2]==res[2] and item[3]-res[3]<20:
                    if item[5]-res[5]<20:
                        res[3] = max(item[3], res[3])
                        res[5]=int((item[5]+res[5])/2) 
                        res[6]=res[3]-res[2]+res[5]
                        merged = True
                        break
        if not merged:
            result.append(item)
    return result

#intra pair
# [ins_num,ele2[al_chr_name],int((ele2[al_r_start]+ele1[al_r_end])/2),delta_length,readname,str(query1[ele1[2]+int(dis_ref/2):ele2[1]-int((dis_ref)/2)])
def intra_ins_del_match(ins_list,del_list,reference):
 #对每个intra ins进行匹配
 #intra排序
 #del 建桶
 #对于大于5000的ins进行匹配
    new_del_list=[]
    new_ins_list=[]
    cut_paste=[]
    if len(ins_list)==0 or len(del_list)==0:
        return ins_list,del_list,cut_paste
    ins_list.sort(key=lambda x: (x[3],[1], x[2]))
    del_list.sort(key=lambda x: (x[1], x[2]))
    del_buget=del_bucket(del_list)
    del_delindex=[]
    for ins_can in ins_list:
        flag=False
        if ins_can[3]>2000:
            chrom=ins_can[1]
            ins_start=ins_can[2]
            ins_end=ins_can[3]+ins_start
            readname=ins_can[4]
            del_result=find_del_in_bucket_forintra(del_buget,chrom,ins_can[3])
            if 0<len(del_result)<5 :
                queryseq=ins_can[5]
                k_list=[]
                bool_list = []
                for del_can in del_result:
                    k, is_true = simi_k(queryseq, del_can, reference)
                    k_list.append(k)
                    bool_list.append(is_true)
                if len(k_list)>0:
                    min_k = min(k_list)
                    if min_k < 15 :
                        for index, k in enumerate(k_list):
                            if k == min_k:
                                del_can = del_result[index]
                                del_index = del_can[0]
                                is_true = bool_list[index]
                                if del_index not in del_delindex:#只保留第一次匹配成功的
                                    del_delindex.append(del_index)
                                    cut_paste.append([-1, del_can[1], del_can[2], del_can[3]+del_can[2], chrom, ins_start, ins_end, is_true, readname, True])
                                    flag = True
                                break

        if not flag:          
            new_ins_list.append(ins_can)
    new_del_list=[del_list[index] for index in range(len(del_list)) if index not in del_delindex]
    # return new_ins_list,new_del_list,cut_paste
    return ins_list,del_list,cut_paste
