from __future__ import print_function

import sys
from statistics import mean
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster


from cute_candidate import *

#全局 repeat区域
#split_read read    q_start q_end   mapq    r_start r_end   is_reversed chr_name    aln_num ref_num
#             0         1     2       3         4     5         6           7           8     9
# from test1 import plot
from collections import Counter
from enum import Enum



from statistics import mean
import numpy as np



def is_similar(chr1, start1, end1, chr2, start2, end2):
    if chr1 == chr2 and abs(start1 - start2) < 20 and abs(end1 - end2) < 20:
        return True
    else:
        return False



from cluster import cluster_del,cluster_ins,cluster_dup

Max_ref_gap=500#500
Max_read_gap=500#500
Max_ref_overlap=3000#3000
Max_read_overlap=300#300#ins del
Max_dup_readoverlap=3000
Min_sv_size=30
Max_sv_Size=200000
al_r_start=3
al_r_end=4
al_q_start=1
al_q_end=2
al_ref_num=9
al_aln_num=8
al_is_reverse=5
al_chr_name=0
types_to_output=['DEL','INS','INV','DUP_TAN','DUP_INT','BND']
#                   #0  #1    #2    #3    #4
del_num=0
inv_num=2
ins_num=1
duptan_num=3
dupint_num=4
bnd_num=5
    #0      1   2       3           4       5
#   #ins    chr bp      svlen   read    sequence
    #del    chr bp      svlen   read
    #inv    chr start   end     read
    #dup_tan    chr start   end     read
    #dup_int    chr start   end     read
    #bnd    chr start   chr     start   read
#[refename, alignment.query_alignment_start,alignment.query_alignment_end,alignment.reference_start,alignment.reference_end,alignment.is_reverse]
def aln_filter1(new_alns):#去掉有重叠的 最大化比对
    # 假设 new_alns 是输入的列表，newalns 是输出的列表
    newalns = []
    alni = 1
    new_alns.sort(key=lambda x:(x[1],x[2]))
    # with open('aln_before.txt','a') as f:
    #     wite_alns=[(aln[0],aln[6],aln[2]-aln[1],aln[1],aln[2],aln[3],aln[4],aln[5]) for aln in new_alns]
    #     print(len(wite_alns),file=f)
    #     print(wite_alns,file=f)
    # plot(new_alns,'h2tg56lbefore')
    newalns.append(new_alns[0])
    while alni < len(new_alns):
        if alni+1<len(new_alns) and new_alns[alni+1][1]-newalns[-1][2]<500:
            alni+=1
        else:
            if new_alns[alni][2]-newalns[-1][2]<100 :
                alni+=1
                continue
            newalns.append(new_alns[alni])
            alni+=1
    # plot(newalns,'h2tg56lafter')
    # with open('aln_after.txt','a') as f:
    #     wite_alns=[(aln[0],aln[6],aln[2]-aln[1],aln[1],aln[2],aln[3],aln[4],aln[5]) for aln in newalns]
    #     print(len(wite_alns),file=f)
    #     print(wite_alns,file=f)
    return newalns
def aln_filter(new_alns):#针对INV发现做出的改进 以ref作为底
    #假设1：一条contig上的比对在此区间是尽可能连续的 没有复杂的易位等
    # 假设 new_alns 是输入的列表，newalns 是输出的列表
    new_alns.sort(key=lambda x:(x[1],x[2]))
    newalns = []
    alni = 1
    newalns.append(new_alns[0])
    while alni < len(new_alns):
        #如果在read上有重叠
        if new_alns[alni][2] < newalns[-1][2] :
            if new_alns[alni][0]==newalns[-1][0] :
                if new_alns[alni][3]>newalns[-1][3] :
                    if alni+1 < len(new_alns)-1:
                        if new_alns[alni][4]<new_alns[alni+1][4]:#且在ref上有重叠
                            newalns.append(new_alns[alni]+[new_alns[alni][2]-new_alns[alni][1],new_alns[alni][4]-new_alns[alni][3],new_alns[alni][2]-new_alns[alni][1]-new_alns[alni][4]+new_alns[alni][3]])
                            alni+=1
                            continue
                if alni<len(new_alns)-1 :
                    if new_alns[alni][5]!=newalns[-1][5] or new_alns[alni+1][5]!=newalns[-1][5]:
                        newalns.append(new_alns[alni]+[new_alns[alni][2]-new_alns[alni][1],new_alns[alni][4]-new_alns[alni][3],new_alns[alni][2]-new_alns[alni][1]-new_alns[alni][4]+new_alns[alni][3]])
                        alni+=1
                        continue
                    
          
                # print(new_alns[alni])
                alni+=1
                continue
                
        newalns.append(new_alns[alni]+[new_alns[alni][2]-new_alns[alni][1],new_alns[alni][4]-new_alns[alni][3],new_alns[alni][2]-new_alns[alni][1]-new_alns[alni][4]+new_alns[alni][3]])
        alni+=1

    return newalns


def nearby_collect(ele1,ele2,query,query_length, readname,ins_candidates,del_candidates,duptans,dupinters):
    #same chr same strand

    if ele1[al_is_reverse]==ele2[al_is_reverse]:#strand
                query1=query
                if ele2[al_is_reverse]:#反转 修改坐标a 
                    if ele1[al_r_start] > ele2[al_r_start]: 
                        q_start=query_length-ele1[al_q_end]
                        q_end=query_length-ele1[al_q_start]
                        ele1[al_q_start],ele1[al_q_end]=q_start,q_end

                        q_start=query_length-ele2[al_q_end]
                        q_end=query_length-ele2[al_q_start]
                        ele2[al_q_start],ele2[al_q_end]=q_start,q_end
                        e=ele1[:]
                        ele1=ele2[:]
                        ele2=e[:]
                        query1=query[::-1]
                dis_ref = ele2[al_r_start] - ele1[al_r_end]
                dis_read = ele2[al_q_start] - ele1[al_q_end]
                delta_length = dis_read-dis_ref
                
                overlap_rlength= ele1[al_r_end]-ele2[al_r_start]
                if -Max_read_overlap<=dis_read:#无 overlap on read
                    if overlap_rlength <=Max_ref_overlap:#无overlap on ref 
                        if Max_sv_Size>=delta_length >=Min_sv_size:#ins
                            if dis_ref<=Max_ref_gap:#ref 无gap
                                ins_candidates.append([ins_num,ele2[al_chr_name],int((ele2[al_r_start]+ele1[al_r_end])/2),delta_length,readname,str(query1[ele1[2]+int(dis_ref/2):ele2[1]-int((dis_ref)/2)])])
                        elif Max_sv_Size>=-delta_length >=Min_sv_size:#del
                            if dis_read<Max_read_gap:#read 无gap
                                del_candidates.append([del_num,ele1[al_chr_name],ele1[al_r_end],-delta_length,readname])
                        # elif abs(delta_length)>Max_sv_Size:
                        #trans same#或者是直接为参考基因组上平常区域
                            # dupinters.append( [ele1[al_is_reverse],ele2[al_is_reverse],ele1[al_chr_name], ele1[al_r_end],ele1[al_chr_name], ele2[al_r_start],readname,'del/ins'])
                        #存不存在一种 ref差距不大 dup大的情况 是否可以看作dup
def find_duptan(ele1,ele2,query_length, readname,duptans,dupinters):
    dup_tan=[]
    if ele1[al_is_reverse]==ele2[al_is_reverse]:#strand
                
        if ele2[al_is_reverse]:#反转 修改坐标a 
            if ele1[al_r_start] > ele2[al_r_start]: 
                q_start=query_length-ele1[al_q_end]
                q_end=query_length-ele1[al_q_start]
                ele1[al_q_start],ele1[al_q_end]=q_start,q_end

                q_start=query_length-ele2[al_q_end]
                q_end=query_length-ele2[al_q_start]
                ele2[al_q_start],ele2[al_q_end]=q_start,q_end
                e=ele1[:]
                ele1=ele2[:]
                ele2=e[:]

        dis_ref = ele2[al_r_start] - ele1[al_r_end]
        dis_read = ele2[al_q_start] - ele1[al_q_end]
        delta_length = dis_read-dis_ref
        
        # if -Max_read_overlap<=dis_read:#无 overlap on read
        if Max_read_gap>=dis_read>=-Max_dup_readoverlap:#无overlap on read 无gap
            if dis_ref<=0:
                if delta_length>=Min_sv_size:
                    
                    if ele1[al_r_start]<ele2[al_r_end]:
                        #case1 dup_complete
                        if ele2[al_r_end]<=ele1[al_r_end]:
                            dup_tan.append([duptan_num,ele1[al_chr_name],ele2[al_r_start],ele2[al_r_end],True,ele2[al_is_reverse],ele2[al_q_start],readname,1])
                        else:#case2 dup_after
                            if Max_sv_Size>=ele1[al_r_end]-ele2[al_r_start]>=Min_sv_size:
                                dup_tan.append([duptan_num,ele1[al_chr_name],ele2[al_r_start],ele1[al_r_end],False,ele2[al_is_reverse],ele2[al_q_start],readname,2])
                    else:
                        #足够近标为dup;否则则是trans
                        if dis_ref>=-Max_sv_Size:#case3 dup_before
                            dup_tan.append([duptan_num,ele2[al_chr_name],ele2[al_r_start],ele2[al_r_start]+delta_length,False,ele2[al_is_reverse],ele2[al_q_start],readname,3])
                        # else:#认为可能发生了易位
                        #     dupinters.append([ele1[al_is_reverse],ele2[al_is_reverse],ele1[al_chr_name],ele1[al_r_start],ele2[al_chr_name],ele2[al_r_start],readname,'dup_int_before'])
    if len(dup_tan)!=0:
        duptans.extend(dup_tan)
        return True
    return False


def analysis_split_read(primary, supplementaries,sv_candidates,  readname,header,options):
    bnd_candidates=[]
    inv_candidates=[]
    invs=[]#存储inv的中间结果 用以配对
    duptans=[]
    dupinters=[]
    ins_candidates=[]
    del_candidates=[]
    query=primary.query_sequence
    query_length=primary.infer_read_length()
    alignments = [primary] + supplementaries
    alignment_list = []

  
    
    for alignment in alignments:#
        try:
            refename=alignment.reference_name
        except:
            refename= header['SQ'][alignment.reference_id]['SN']
        alignment_list.append([refename, alignment.query_alignment_start,alignment.query_alignment_end,alignment.reference_start,alignment.reference_end,alignment.is_reverse,alignment.mapping_quality])
    # filter在内部的alignment等
    # aln_filter2(alignment_list)
    alignment_list=aln_filter(alignment_list)
    if len(alignment_list)>100:#80
        return [],[],[]
    flag=False
    max_refpos=0
    max_alnpos=0
    thred=200#默认为下一个片段至少要大多少为可信高
    max_thred=2000000#否则就认为是易位了
        #进入复杂变异处理模式
    # alignment_list.sort(key=lambda x:(x[al_q_start],x[al_q_end]))

    #初始定义
    #记录遍历max位置
    refpos_dic={}
    paste=[]
    bnd=[]
    chr_list=[aln[0] for aln in alignment_list]
    chr_count = Counter(chr_list)
    # 获取出现次数最多的元素及其出现次数
    most_common_chr = chr_count.most_common(1)[0][0]
    for index in range(1, len(alignment_list)):
        current_aln = alignment_list[index]
        flag=False
        

        if current_aln[0] not in refpos_dic:
            refpos_dic[current_aln[0]] = []
            refpos_dic[current_aln[0]].append((current_aln[3], current_aln[4], index, current_aln[0]))
        else:
            # if last_ref!=current_aln[0]:
            #     flag=True
            lastindex=refpos_dic[current_aln[0]][-1][2]
            prev_aln = alignment_list[lastindex]
            max_refpos = refpos_dic[current_aln[0]][-1][1]
            same_chr = current_aln[0]==alignment_list[index-1][0]
            main_chr=True if most_common_chr==current_aln[0] else False
            # max_alnpos=alignment_list[index][al_q_end]
            # 聚焦于aln的变化
            aln_growth = current_aln[al_q_end] > prev_aln[al_q_end]+thred
            

            if main_chr:
                ref_growth = current_aln[al_r_end] > max_refpos+thred
                
                ele = prev_aln
                if aln_growth:
                    
                    detlength = current_aln[al_q_start] - ele[al_q_end]
                    if ref_growth:
                    
                        if current_aln[al_r_start] - max_refpos - (current_aln[al_q_start] - max_alnpos) > max_thred:
                            # aln + ref++ 易位
                            
                            paste.append([0, ele[al_chr_name], ele[al_is_reverse], ele[al_r_start], ele[al_r_end],
                                        current_aln[al_chr_name], current_aln[al_is_reverse], current_aln[al_r_start], current_aln[al_r_end],
                                        detlength, readname])
                            bnd.append([1,readname, ele[al_chr_name], ele[al_r_end]+detlength,index, current_aln[al_chr_name], current_aln[al_r_start], 
                                    current_aln[al_r_end],'++', current_aln[al_r_end]- current_aln[al_r_start]])
                        else:
                            nearby_collect(prev_aln[:], current_aln[:], query, query_length, readname, ins_candidates,
                                        del_candidates, duptans, dupinters)
                            refpos_dic[current_aln[0]].append((current_aln[3], current_aln[4], index, current_aln[0]))
                    else:  # ref 变小
                        # flag=find_duptan(prev_aln[:], current_aln[:],  query_length, readname,duptans, dupinters)
                        # if not flag:
                        paste.append([1, current_aln[al_chr_name], current_aln[al_is_reverse], current_aln[al_r_start], current_aln[al_r_end],
                                    ele[al_chr_name], ele[al_is_reverse], ele[al_r_start], ele[al_r_end], detlength,
                                    readname])
                        bnd.append([0, readname, current_aln[al_chr_name], current_aln[al_r_start], current_aln[al_r_end],
                                        lastindex, ele[al_chr_name], ele[al_r_end]+detlength,'+-', current_aln[al_r_end]- current_aln[al_r_start]])
                    max_alnpos=current_aln[al_q_end]
                else:  # aln 变小
                    deata=current_aln[al_q_start]-ele[al_q_start]
                    if ref_growth:
                        bnd.append([1,readname, ele[al_chr_name], ele[al_r_start]+deata,index, current_aln[al_chr_name], current_aln[al_r_start], 
                                    current_aln[al_r_end],'-+', current_aln[al_r_end]- current_aln[al_r_start]])
                        
                    else:  # ref 变小
                        
                        bnd.append([0, readname, current_aln[al_chr_name], current_aln[al_r_start], current_aln[al_r_end],
                                        lastindex, ele[al_chr_name], ele[al_r_start]+deata,'--', current_aln[al_r_end]- current_aln[al_r_start]])
            if not same_chr:  # 不一样染色体
                if current_aln[al_q_end]>max_alnpos+thred:
                    max_alnpos=current_aln[al_q_end]
       
    for index in range(len(alignment_list)-1):
        if flag:
            flag=False
            continue

        ele1 = alignment_list[index][:]
        ele2 = alignment_list[index + 1][:]
        if ele2[al_chr_name]==ele1[al_chr_name]:#chr
            if ele1[al_is_reverse]==ele2[al_is_reverse]:#strand
                find_duptan(ele1,ele2,query_length,readname,duptans,dupinters)
            else:#inv 只判断出-+-  +-+类型 其他则被认为是trans inv 及后续聚类+---+ -+++-类型
                #new inv 规律性 大块内会有一小部分再次比对到中间或者附近 将它作为正确的指示信号
                #非常远的倒转作为FP -----过滤
                if ele1[al_is_reverse]:
                    if index+2<len(alignment_list):
                        ele3=alignment_list[index+2][:]
                        if  ele3[al_is_reverse] and ele3[al_chr_name]==ele1[al_chr_name]:
                            
                            if ele2[al_q_start] + 0.5 * (ele1[al_r_start] - ele3[al_r_end]) >= ele3[al_q_end] and ele1[al_q_start] + 0.5 * (ele1[al_r_start] - ele3[al_r_end]) >= ele2[al_q_end]:#正常在中间的iver
                                # No overlaps in split reads

                                if ele2[al_r_start]-ele3[al_r_end]>= -Min_sv_size and ele1[al_r_start] - ele2[al_r_end]>=-Min_sv_size:#这个和下边的不太一样
                                    inv_candidates.append([inv_num,ele1[al_chr_name],ele3[al_r_end],ele2[al_r_end],'-+-',readname])
                                else:#说明是移位反转
                                    bnd_candidates.append([bnd_num,ele1[al_chr_name],ele3[al_r_end],ele2[al_chr_name],ele2[al_r_end],readname,'INV bnd -+-'])
                                    inv_candidates.append([inv_num,ele2[al_chr_name],ele2[al_r_start],ele2[al_r_end],'-+-',readname])
                            else:#易位
                                # dupinters.append([ele1[al_is_reverse],ele2[al_is_reverse],ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_start],readname,'inv_bnd -+-',index])
                                invs.append((ele2[3],ele2[4],False,index+1,ele1[0]))
                                flag=True
                        else:#-+++---
                            # invs.append((ele1[3],ele1[4],True,index,ele1[0]))#置空
                            refs=refpos_dic[ele1[al_chr_name]]
                            refs_index=[i[2] for i in refs]
                            if index in refs_index :
                                if index+1 not in refs_index:
                                    invs.append((ele2[3],ele2[4],False,index+1,ele1[0]))
                                    flag=True
                            else:
                                if index+1 in refs_index:
                                    invs.append((ele1[3],ele1[4],True,index,ele1[0]))
                        continue         

                        
                    if abs(ele1[al_q_end]-ele2[al_q_start])<=100 and ele2[al_r_start]-ele1[al_r_start]<=Max_sv_Size:
                        if ele1[al_r_start]<ele2[al_r_start]:
                            inv_candidates.append([inv_num,ele1[al_chr_name],ele1[al_r_start],ele2[al_r_start],'-+',readname])
                        # else:
                        #     print('inv',end='\t')
                        #     print([inv_num,ele1[al_chr_name],ele1[al_r_start],ele2[al_r_start],'-+',readname])
                    # else:
                    #     dupinters.append([ele1[al_is_reverse],ele2[al_is_reverse],ele1[al_chr_name],ele1[al_r_start],ele2[al_chr_name],ele2[al_r_start],readname,'inv_-+',index])

                elif ele2[al_is_reverse]:#+-
                    if index+2<len(alignment_list):
                        ele3=alignment_list[index+2][:]
                        if not ele3[al_is_reverse] and ele3[al_chr_name]==ele1[al_chr_name]:
                            #如果内涵
                            if ele1[1]<ele2[1] and ele2[2]<ele1[2] and  ele2[3]>ele1[3] and ele2[4]<ele3[4]:
                                inv_candidates.append([inv_num,ele2[al_chr_name],ele2[al_r_start],ele2[al_r_end],'+-+ in',readname])
                                flag=True
                                continue
                                
                            if ele2[al_q_start] + 0.5 * (ele3[al_r_start] - ele1[al_r_end]) >= ele1[al_q_end] and ele3[al_q_start] + 0.5 * (ele3[al_r_start] - ele1[al_r_end]) >= ele2[al_q_end]:#两边没有太大的gap
                                # No overlaps in split reads

                                if ele2[al_r_start]-ele1[al_r_end]>= -Min_sv_size and ele3[al_r_start] - ele2[al_r_end]>=-Min_sv_size:#保证比对也在中间
                                    inv_candidates.append([inv_num,ele1[al_chr_name],ele1[al_r_end],ele2[al_r_end],"+-+",readname])
                                    
                                   
                                else:#说明是移位反转
                                    bnd_candidates.append([bnd_num,ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_start],readname,'INV bnd +-+'])
                                    inv_candidates.append([inv_num,ele2[al_chr_name],ele2[al_r_start],ele2[al_r_end],"+-+",readname])
                            else:#易位
                                # dupinters.append([ele1[al_is_reverse],ele2[al_is_reverse],ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_start],readname,'inv_bnd +-+',index])
                                invs.append((ele2[3],ele2[4],True,index+1,ele1[0]))
                                flag=True
                        else:#+--++++
                            refs=refpos_dic[ele1[al_chr_name]]
                            refs_index=[i[2] for i in refs]
                            if index in refs_index :
                                if index+1 not in refs_index:
                                    invs.append((ele2[3],ele2[4],True,index+1,ele1[0]))
                                    flag=True
                            else:
                                if index+1 in refs_index:
                                    invs.append((ele1[3],ele1[4],False,index,ele1[0]))
                        continue
                    if abs(ele1[al_q_end]-ele2[al_q_start])<=100 and ele2[al_r_end]-ele1[al_r_end]<=Max_sv_Size:
                        if ele1[al_r_end] < ele2[al_r_end]:

                            inv_candidates.append([inv_num,ele1[al_chr_name],ele1[al_r_end],ele2[al_r_end],'+-',readname])
                    #     else:
                    #         print('inv',end='\t')
                    #         print([inv_num,ele1[al_chr_name],ele1[al_r_end],ele2[al_r_end],'+-',readname])
                    # else:
                    #     dupinters.append([ele1[al_is_reverse],ele2[al_is_reverse],ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_end],readname,'inv_+-',index])

                    # print("+--")
               
               
        else:#dif chr trans
            length=min(ele1[al_r_end]-ele1[al_r_start],ele2[al_r_end]-ele2[al_r_start])
            if abs(ele1[al_q_end]-ele2[al_q_start])<=1000:
                
                bnd_candidates.append([5,ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_start],readname,length,False,False])
            elif  abs(ele2[al_q_end]-ele1[al_q_start])<=1000:
                bnd_candidates.append([5,ele1[al_chr_name],ele1[al_r_start],ele2[al_chr_name],ele2[al_r_end],readname,length,False,False])
            elif ele2[al_q_end]<=ele1[al_q_end]:#内部的bnd
                if index+2<len(alignment_list):
                    if alignment_list[index+2][al_chr_name]==ele1[al_chr_name]:
                        bnd.append([index+1,ele2[al_chr_name],ele2[al_r_start],ele2[al_r_end],ele1[al_chr_name],ele1[al_r_start]+ele2[al_q_start]-ele1[al_q_start],readname,2,ele2[al_r_end]-ele2[al_r_start]])
            #             flag=True
            # else:
                
            #     # bnd.append([index,ele1[al_chr_name],ele1[al_r_end]+ele2[al_q_start]-ele1[al_q_end],ele2[al_chr_name],ele2[al_r_start],readname,2,ele2[al_r_end]-ele2[al_r_start]])
            #     print("else bnd",index)#跨过的bnd#['chr6_GL000253v2_alt', 178623, 1382982, 178623, 1382982, False, 60, 1204359, 1204359, 0]['chr6', 1214527, 2041800, 29953766, 30731549, False, 60, 827273, 777783, 49490]
    #invs处理
    if len(invs)>0:
        for inv in invs:
        
            start,end,flag,index,chr=inv
            refs_refe = refpos_dic[chr]
            refs_index=[i[2] for i in refs_refe]
            ref_index=0
            for s,e,i,c in refs_refe:
                if start<s :
                    break
                ref_index+=1
            if end-e<Max_ref_overlap:
                id = None
                for id in range(0, index)[::-1]:
                    if id in refs_index:
                        break
                if id is not None:
                    if alignment_list[id][al_r_start]<start and alignment_list[id][al_r_end]>end:#在上一个比对中间
                        if alignment_list[id][al_is_reverse]!=flag:
                            inv_candidates.append([inv_num,chr,start,end,'invs',readname])
                            continue
                    if ref_index-1>=0:
                        if alignment_list[refs_refe[ref_index-1][2]][al_is_reverse] == alignment_list[i][al_is_reverse] :#内部相反的 +-+ -+- 处理不了易位反转
                            if flag!=alignment_list[i][al_is_reverse]:
                                # inv_candidates.append([inv_num,chr,start,end,'invs',readname])
                                if i-1!=id:#有易位
                                    ele=alignment_list[index]
                                    ele_before=alignment_list[id]
                                    delath=ele[al_q_start]-ele_before[al_q_start]
                                    bnd.append([index,chr,start,end,chr,ele_before[al_r_start]+delath,readname,3,end-start])

                                
    # if len(dupinters)>0:
    #     print("dupinters",len(dupinters))
    # if len(bnd)>0:
    #     print("bnd",len(bnd))
    
    sv_candidates[types_to_output[del_num]].extend(cluster_del(del_candidates))
    sv_candidates[types_to_output[ins_num]].extend(cluster_ins(ins_candidates))
    
    sv_candidates[types_to_output[bnd_num]].extend(bnd_candidates)
    sv_candidates[types_to_output[inv_num]].extend(inv_candidates)
    return paste,duptans,bnd
