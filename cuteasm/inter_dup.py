from __future__ import print_function

import sys
from statistics import mean
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster


from cute_candidate import *



def analysis_inv(ele_1, ele_2, read_name, candidate, SV_size,bam):
    if not ele_1[5] :#什么意思啊
        # +--
        if ele_1[3] - ele_2[3] >= SV_size:
            if ele_2[0] + 0.5 * (ele_1[3] - ele_2[3]) >= ele_1[1]:
                
                candidate[types_to_output[inv_num]].append([ele_1[4],ele_2[3],ele_1[3],"++",read_name])
                # head-to-head
                # 5'->5'
        if ele_2[3] - ele_1[3] >= SV_size:
            if ele_2[0] + 0.5 * (ele_2[3] - ele_1[3]) >= ele_1[1]:
             
                candidate[types_to_output[inv_num]].append([ele_1[4],ele_1[3],ele_2[3],"++",read_name])
                # head-to-head
                # 5'->5'
    else:
        # -++
        if ele_2[2] - ele_1[2] >= SV_size:
            if ele_2[0] + 0.5 * (ele_2[2] - ele_1[2]) >= ele_1[1]:
               
                candidate[types_to_output[inv_num]].append([ele_1[4],ele_1[2],ele_2[2],"--",read_name])
                # tail-to-tail
                # 3'->3'
        if ele_1[2] - ele_2[2] >= SV_size:
            if ele_2[0] + 0.5 * (ele_1[2] - ele_2[2]) >= ele_1[1]:
              
                candidate[types_to_output[inv_num]].append([ele_1[4],ele_2[2],ele_1[2],"--",read_name])
                # tail-to-tail
                # 3'->3'


def analysis_bnd(ele_1, ele_2, read_name, candidate,bam):
    '''
    *********Description*********
    *	TYPE A:		N[chr:pos[	*
    *	TYPE B:		N]chr:pos]	*
    *	TYPE C:		[chr:pos[N	*
    *	TYPE D:		]chr:pos]N	*
    *****************************
    '''
    a=ele_1[-1].reference_id
    b=ele_2[-1].reference_id
    if a==b:
        a=ele_1[2]
        b=ele_2[2]
    if ele_2[0] - ele_1[1] <= 100:
        if not ele_1[5]:
            if not ele_2[5]:
                # +&+
                if a < b:
                    candidate.append([ele_1[4],ele_1[3],ele_2[4],ele_2[2],'A',read_name])
                    # N[chr:pos[
                else:
                    candidate.append([ele_2[4],ele_2[2],ele_1[4],ele_1[3],'D',read_name])
                    # ]chr:pos]N
            else:
                # +&-
                if a < b:
                    
                    candidate.append([ele_1[4],ele_1[3],ele_2[4],ele_2[3],'B',read_name])
                    # N]chr:pos]
                else:
                    
                    candidate.append([ele_2[4],ele_2[3],ele_1[4],ele_1[3],'B',read_name])
                    # N]chr:pos]
        else:
            if not  ele_2[5]:
                # -&+
                if a <b:
                  
                    candidate.append([ele_1[4],ele_1[2],ele_2[4],ele_2[2],'C',read_name])
                    # [chr:pos[N
                else:
                    
                    candidate.append([ele_2[4],ele_2[2],ele_1[4],ele_1[2],'C',read_name])
                    # [chr:pos[N
            else:
                # -&-
                if a < b:
                    
                    candidate.append([ele_1[4],ele_1[2],ele_2[4],ele_2[3],'D',read_name])
                    # ]chr:pos]N
                else:
                    
                    candidate.append([ele_2[4],ele_2[3],ele_1[4],ele_1[2],'A',read_name])
                    # N[chr:pos[
def analysis_split_read(split_read,RLength,read_name,SV_size, MaxSize, query,bam):
    #转化为输入list
    #输出 candidate SVs
    '''
    read_start	read_end	ref_start	ref_end	chr	strand aln
    #0			#1			#2			#3		#4	#5     #6
    '''
    SP_list=sorted(split_read,key=lambda x:x[0])
    trigger_INS_TRA = 0
    candidate=[]
   
   
    # Store Strands of INV
    if len(SP_list) == 2:
        ele_1 = SP_list[0]
        # q_start_e1=ele_1[0]
        # q_end_e1=ele_1[1]
        # ref_start_e1=ele_1[2]
        # ref_start_e1=ele_1[3]
        # chr_e1=ele_1[4]
        # strand_e1=ele_1[5]
        # aln_e1=ele_1[6]
        ele_2 = SP_list[1]
        # q_start_e2=ele_2[0]
        # q_end_e2=ele_2[1]
        
        if ele_1[4] == ele_2[4]:#chr
            if ele_1[5] != ele_2[5]:#strand
                analysis_inv(ele_1, 
                                ele_2, 
                                read_name, 
                                candidate,
                                SV_size,bam)

            else:
                # dup & ins & del 
                a = 0
                if ele_1[5] :
                    ele_1 = [RLength-SP_list[a+1][1], RLength-SP_list[a+1][0]]+SP_list[a+1][2:]
                    ele_2 = [RLength-SP_list[a][1], RLength-SP_list[a][0]]+SP_list[a][2:]
                    query = query[::-1]
#     read_start	read_end	ref_start	ref_end	chr	strand aln
#     #0			#1			#2			#3		#4	#5     #6
                if ele_1[3] - ele_2[2] >= SV_size:#ref 间隔 前后有交叠
                    # if ele_2[1] - ele_1[1] >= ele_1[3] - ele_2[2]:
                    if ele_2[0] - ele_1[1] >= ele_1[3] - ele_2[2]:#query上的overlap大于ref上
                        candidate.append(CandidateInsertion(ele_2[4],(ele_1[3]+ele_2[2])//2,ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1],read_name,str(query[ele_1[1]+int((ele_1[3]-ele_2[2])/2):ele_2[0]-int((ele_1[3]-ele_2[2])/2)])))
                    else:
                        
                        candidate.append(CandidateDuplication(ele_2[4],ele_2[2],ele_1[3],read_name,bam))

                delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                    if ele_2[2] - ele_1[3] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                        candidate.append(CandidateInsertion(ele_2[4],(ele_2[2]+ele_1[3])//2,delta_length,read_name,str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)])))
                delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                    if ele_2[0] - ele_1[1] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                        candidate.append(CandidateDeletion(ele_2[4],ele_1[3],delta_length,read_name))
    
        else:
            trigger_INS_TRA = 1
            analysis_bnd(ele_1, ele_2, read_name, candidate,bam)

    else:
        # over three splits
        
        for a in range(len(SP_list[1:-1])):
            ele_1 = SP_list[a]
            ele_2 = SP_list[a+1]
            ele_3 = SP_list[a+2]
#     read_start	read_end	ref_start	ref_end	chr	strand aln
#     #0			#1			#2			#3		#4	#5     #6
# '''
            if ele_1[4] == ele_2[4]:#chr
                if ele_2[4] == ele_3[4]:#chr
                    if ele_1[5] == ele_3[5] and ele_1[5] != ele_2[5]:#有strand不一样
                        if ele_2[5] :
                            # +-+
                            if ele_2[0] + 0.5 * (ele_3[2] - ele_1[3]) >= ele_1[1] and ele_3[0] + 0.5 * (ele_3[2] - ele_1[3]) >= ele_2[1]:
                                # No overlaps in split reads

                                if ele_2[2] >= ele_1[3] and ele_3[2] >= ele_2[3]:
                                    candidate.append(CandidateInversion(ele_1[4],ele_1[3],ele_2[3],"++",read_name,bam))
                                    # head-to-head
                                    # 5'->5'
                                    candidate.append(CandidateInversion(ele_1[4],ele_2[2],ele_3[2],"--",read_name,bam))
                                    # tail-to-tail
                                    # 3'->3'
                        else:
                            # -+-
                            if ele_1[1] <= ele_2[0] + 0.5 * (ele_1[2] - ele_3[3]) and ele_3[0] + 0.5 * (ele_1[2] - ele_3[3]) >= ele_2[1]:
                                # No overlaps in split reads

                                if ele_2[2] - ele_3[3] >= -50 and ele_1[2] - ele_2[3] >= -50:
                                   
                                    candidate.append(CandidateInversion(ele_1[4],ele_3[3],ele_2[3],"++",read_name,bam))
                                    # head-to-head
                                    # 5'->5'
                                    candidate.append(CandidateInversion(ele_1[4],ele_2[2],ele_1[2],"--",read_name,bam))
                                    # tail-to-tail
                                    # 3'->3'	

                    if len(SP_list) - 3 == a:
                        if ele_1[5] != ele_3[5]:
                            if ele_2[5] == ele_1[5]:
                                # ++-/--+
                                analysis_inv(ele_2, 
                                                ele_3, 
                                                read_name, 
                                                candidate, 
                                                SV_size,bam)
                            else:
                                # +--/-++
                                analysis_inv(ele_1, 
                                                ele_2, 
                                                read_name, 
                                                candidate, 
                                                SV_size,bam)

                    if ele_1[5] == ele_3[5] and ele_1[5] == ele_2[5]:#
                        # dup & ins & del 
                        if ele_1[5] :
                            ele_1 = [RLength-SP_list[a+2][1], RLength-SP_list[a+2][0]]+SP_list[a+2][2:]
                            ele_2 = [RLength-SP_list[a+1][1], RLength-SP_list[a+1][0]]+SP_list[a+1][2:]
                            ele_3 = [RLength-SP_list[a][1], RLength-SP_list[a][0]]+SP_list[a][2:]
                            query = query[::-1]

                        if ele_2[3] - ele_3[2] >= SV_size and ele_2[2] < ele_3[3]:
                         
                            candidate.append(CandidateDuplication(ele_2[4],ele_3[2],ele_2[3],read_name,bam))


                        if a == 0:
                            if ele_1[3] - ele_2[2] >= SV_size:
                                candidate.append(CandidateDuplication(ele_2[4],ele_2[2],ele_1[3],read_name,bam))

                        delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                        if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                            if ele_2[2] - ele_1[3] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                if ele_3[2] >= ele_2[3]:
                                    candidate.append(CandidateInsertion(ele_2[4],(ele_2[2]+ele_1[3])//2,delta_length,read_name,str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)])))
                        delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                        if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                            if ele_2[0] - ele_1[1] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                if ele_3[2] >= ele_2[3]:
                                    candidate.append(CandidateDeletion(ele_2[4],ele_1[3],delta_length,read_name))
                        
                        if len(SP_list) - 3 == a:
                            ele_1 = ele_2
                            ele_2 = ele_3

                            delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                            if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                                if ele_2[2] - ele_1[3] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                    candidate.append(CandidateInsertion(ele_2[4],(ele_2[2]+ele_1[3])//2,delta_length,read_name,str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)])))

                            delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                            if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
                                if ele_2[0] - ele_1[1] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                   
                                    candidate.append(CandidateDeletion(ele_2[4],ele_1[3],delta_length,read_name))

                    if len(SP_list) - 3 == a and ele_1[5] != ele_2[5] and ele_2[5] == ele_3[5]:
                        ele_1 = ele_2
                        ele_2 = ele_3
                        ele_3 = None
                    if ele_3 == None or (ele_1[5] == ele_2[5] and ele_2[5] != ele_3[5]):
                        if ele_1[5]:
                            ele_1 = [RLength-SP_list[a+1][1], RLength-SP_list[a+1][0]]+SP_list[a+1][2:]
                            ele_2 = [RLength-SP_list[a][1], RLength-SP_list[a][0]]+SP_list[a][2:]
                            query = query[::-1]
                        delta_length = ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]
                        if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                            if ele_2[2] - ele_1[3] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                                candidate.append(CandidateInsertion(ele_2[4],(ele_2[2]+ele_1[3])//2,delta_length,read_name,str(query[ele_1[1]+int((ele_2[2]-ele_1[3])/2):ele_2[0]-int((ele_2[2]-ele_1[3])/2)])))

                        delta_length = ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
                        if ele_1[3] - ele_2[2] < max(SV_size, delta_length/5) and delta_length >= SV_size:
                            if ele_2[0] - ele_1[1] <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):
                               
                                candidate.append(CandidateDeletion(ele_2[4],ele_1[3],delta_length,read_name))

            else:
                trigger_INS_TRA = 1
                analysis_bnd(ele_1, ele_2, read_name, candidate,bam)

                if len(SP_list) - 3 == a:
                    if ele_2[4] != ele_3[4]:
                        analysis_bnd(ele_2, ele_3, read_name, candidate,bam)

    if len(SP_list) >= 3 and trigger_INS_TRA == 1:#假设中间几段是碎的
        if SP_list[0][4] == SP_list[-1][4]:
            if SP_list[0][5] != SP_list[-1][5]:
                pass
            else:
                if SP_list[0][5] == '+':
                    ele_1 = SP_list[0]
                    ele_2 = SP_list[-1]
                else:
                    ele_1 = [RLength-SP_list[-1][1], RLength-SP_list[-1][0]]+SP_list[-1][2:]
                    ele_2 = [RLength-SP_list[0][1],RLength-SP_list[0][0]]+SP_list[0][2:]
                    query = query[::-1]
                dis_ref = ele_2[2] - ele_1[3]
                dis_read = ele_2[0] - ele_1[1]
                if dis_ref < 100 and dis_read - dis_ref >= SV_size and (dis_read - dis_ref <= MaxSize or MaxSize == -1):
                    candidate.append(CandidateInsertion(ele_2[4],min(ele_2[2], ele_1[3]),dis_read - dis_ref,read_name,str(query[ele_1[1]+int(dis_ref/2):ele_2[0]-int(dis_ref/2)])))

                if dis_ref <= -SV_size:
                    candidate.append(CandidateDuplication(ele_2[4],ele_2[2],ele_1[3],read_name,bam))
    return candidate
def print_pos_sort(sorted_reads):
    for i in sorted_reads:
        print(i[-1],end="\t")
    print("Hhh")


def max_antichain(sequence):
    # 获取序列的长度
    n = len(sequence)
    
    # 创建一个动态规划数组，初始化为1
    dp = [1] * n
    
    # 计算每个元素的最大偏序长度
    for i in range(n):
        for j in range(i):
            if sequence[i] > sequence[j]:  # 如果可以比较
                dp[i] = max(dp[i], dp[j] + 1)
    
    # 找到最大偏序长度
    max_length = max(dp)
    
    # 反向构建最大偏序子序列
    antichain = []
    for i in range(n - 1, -1, -1):
        if dp[i] == max_length:
            antichain.append(sequence[i])
            max_length -= 1
    
    # 由于是反向构建，最后需要反转结果
    return antichain[::-1]



def analysis_trans(a1,a2,read_name,ref):
    candidate=[]
    ele_1=[a1[1],a1[2],a1[4],a1[5],a1[7],a1[6],a1[0]]
    ele_2=[a2[1],a2[2],a2[4],a2[5],a2[7],a2[6],a2[0]]
    #     read_start	read_end	ref_start	ref_end	chr	strand aln
    #     #0			#1			#2			#3		#4	#5     #6
   
    analysis_bnd(ele_1,ele_2,read_name,candidate,ref)
    
    return candidate

def detect_reverse(revese_list,is_reverse):
    list_len=len(revese_list)
    most_common_element = max(set(revese_list), key=revese_list.count)
    count = revese_list.count(most_common_element)
    if count/list_len>0.5:
        r_l=[0 if i == most_common_element else 1 for i in revese_list]
    else :
        r_l=[0 if i == is_reverse else 1 for i in revese_list]
    return r_l



#全局 repeat区域
#split_read read    q_start q_end   mapq    r_start r_end   is_reversed chr_name    aln_num ref_num
#             0         1     2       3         4     5         6           7           8     9
from test1 import plot
from collections import Counter
from edlib import align
from enum import Enum
 #ref_start  ref_end ref_num has_aln aln_num_list
    #   0           1       2       3       (4)
class POINT_TYPE(Enum):
    START = 1
    END = 2

class Point:
    def __init__(self, v: int, t: POINT_TYPE):
        self.val = v
        self.type = t
def resolve_dup(reads,dup):
    #需要对里面相邻的两个判断是否是有del或ins 以及dup
    #之后需要对相同的聚合 dup
    #case 1:<=3
    if len(reads)<=3:
        #case 1 存在部分交叠
        #case 2 全部相等
        #case 3 包含
        pass
    #寻找最短的交叠区域吗
    dupunit=[]
    reads.sort(key=lambda x:x[0] )
    res=0
    # for i in range(len(reads)-1):
    #     if 

        

    data=[(x[4],x[5]) for x in reads]
    tuple_counts = Counter(data)
    for tuple_key, count in tuple_counts.items():
        if count>=2:
            dup.append(tuple_key)
        



from statistics import mean
import numpy as np



def is_similar(chr1, start1, end1, chr2, start2, end2):
    if chr1 == chr2 and abs(start1 - start2) < 20 and abs(end1 - end2) < 20:
        return True
    else:
        return False


def reciprocal_overlap_distance(inversion1, inversion2):
    start1, end1, direction1 = inversion1
    start2, end2, direction2 = inversion2
    #Inversion breakpoints with same direction cannot be joined
    if direction1 == direction2:
        return 1
    #Non-overlapping inversion breakpoints cannot be joined
    if start2 >= end1:
        return 1
    if start1 >= end2:
        return 1
    
    if start2 >= start1:
        overlap = min(end1, end2) - start2
    else:
        overlap = min(end1, end2) - start1
    
    relative_overlap1 = overlap / float(end1 - start1)
    relative_overlap2 = overlap / float(end2 - start2)
    minimum_relative_overlap = min(relative_overlap1, relative_overlap2)
    return 1 - minimum_relative_overlap


def process_overlapping_inversions(active_inversions, query_name, bam):
    if len(active_inversions) < 2:
        clusters = [active_inversions]
    else:
        data = np.array( [[inversion[1], inversion[2], 0 if inversion[3].split("_")[0] == "left" else 1] for inversion in active_inversions])
        Z = linkage(data, method = "complete", metric = reciprocal_overlap_distance)
        cluster_indices = list(fcluster(Z, 0.3, criterion='distance'))
        clusters = [[] for i in range(max(cluster_indices))]
        for inversion_index, cluster_index in enumerate(cluster_indices):
            clusters[cluster_index-1].append(active_inversions[inversion_index])

    inversion_candidates = []
    for cluster in clusters:
        chrom = cluster[0][0]
        start = max([i[1] for i in cluster])
        end = min([i[2] for i in cluster])
        complete = True if len(cluster) > 1 else False
        inversion_candidates.append([chrom, start, end, query_name, complete])
    return inversion_candidates


from cluster import cluster_del,cluster_ins,cluster_dupint,cluster_dup
def analyze_read_segments(primary, supplementaries, bam, options,sv_candidates):
    # sv_candidates={}
    # for svtype in types_to_output:
    #     sv_candidates[svtype]=[]
    read_name = primary.query_name
    alignments = [primary] + supplementaries#将所有的信息聚合在一起 这是本来在一条alignment里的信息
    alignment_list = []
    for alignment in alignments:
        #correct query coordinates for reversely mapped reads 纠正反转的查询坐标
        if alignment.is_reverse:#flag:0X10
            q_start = alignment.infer_read_length() - alignment.query_alignment_end#!!!为什么要这样计算
            q_end = alignment.infer_read_length() - alignment.query_alignment_start
        else:
            q_start = alignment.query_alignment_start
            q_end = alignment.query_alignment_end

        new_alignment_dict = {  'q_start': q_start, 
                                'q_end': q_end, 
                                'ref_id': alignment.reference_id, 
                                'ref_start': alignment.reference_start, 
                                'ref_end': alignment.reference_end,
                                'is_reverse': alignment.is_reverse  }
        alignment_list.append(new_alignment_dict)

    sorted_alignment_list = sorted(alignment_list, key=lambda aln: (aln['q_start'], aln['q_end']))
    #inferred_read_length = alignments[0].infer_read_length()
    tandem_duplications = []
    translocations = []
    inversions = []

    for index in range(len(sorted_alignment_list) - 1):
        alignment_current = sorted_alignment_list[index]
        alignment_next = sorted_alignment_list[index + 1]

        distance_on_read = alignment_next['q_start'] - alignment_current['q_end']#在queryreads上前后断点相隔距离
        #对每一个align的信号，根据前后信号是否在同一个染色体及方向来划分类别
        #Same chromosome
        if alignment_current['ref_id'] == alignment_next['ref_id']:
            ref_chr = bam.get_reference_name(alignment_current['ref_id'])
            #Same orientation
            if alignment_current['is_reverse'] == alignment_next['is_reverse']:
                #Compute distance on reference depending on orientation
                if alignment_current['is_reverse']:#如果方向相反，则
                    distance_on_reference = alignment_current['ref_start'] - alignment_next['ref_end']
                else:
                    distance_on_reference = alignment_next['ref_start'] - alignment_current['ref_end']
                #No overlap on read
                if distance_on_read >= -options.query_overlap_tolerance:#没有overlap
                    #No overlap on reference
                    if distance_on_reference >= -options.reference_overlap_tolerance:
                        deviation = distance_on_read - distance_on_reference
                        #INS candidate
                        if deviation >= options.min_sv_size:
                            #No gap on reference
                            if distance_on_reference <= options.reference_gap_tolerance:
                                if not alignment_current['is_reverse']:
                                    insertion_seq = primary.query_sequence[alignment_current['q_end']:alignment_current['q_end']+deviation]#他把这一部分作为插入序列会不会不全
                                    sv_candidates[types_to_output[ins_num]].append([ref_chr, alignment_current['ref_end'],deviation, read_name, insertion_seq])
                                else:
                                    insertion_seq = primary.query_sequence[primary.infer_read_length() - alignment_next['q_start']:primary.infer_read_length() - alignment_next['q_start'] + deviation]
                                    sv_candidates[types_to_output[ins_num]].append([ref_chr, alignment_current['ref_start'],deviation, read_name, insertion_seq])
                        #DEL candidate
                        elif -options.max_sv_size <= deviation <= -options.min_sv_size:
                            #No gap on read
                            if distance_on_read <= options.query_gap_tolerance:
                                if not alignment_current['is_reverse']:
                                    sv_candidates[types_to_output[del_num]].append([ref_chr, alignment_current['ref_end'], -deviation,read_name])
                                else:
                                    sv_candidates[types_to_output[del_num]].append([ref_chr, alignment_next['ref_end'],  - deviation, read_name])
                        #Either very large DEL or TRANS#潜在的trans
                        elif deviation < -options.max_sv_size:
                            #No gap on read
                            if distance_on_read <= options.query_gap_tolerance:
                                if not alignment_current['is_reverse']:
                                    sv_candidates[types_to_output[bnd_num]].append([ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_start'], 'fwd', read_name])
                                    translocations.append(('fwd', 'fwd', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_start']))
                                else:
                                    sv_candidates[types_to_output[bnd_num]].append([ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_end'] - 1, 'rev', read_name])
                                    translocations.append(('rev', 'rev', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_end'] - 1))
                    #overlap on reference
                    else:#主要处理重复区域
                        #No gap on read
                        if distance_on_read <= options.query_gap_tolerance:
                            deviation = distance_on_read - distance_on_reference
                            #Tandem Duplication
                            if deviation >= options.min_sv_size:
                                if not alignment_current['is_reverse']:
                                    #Tandem Duplication (fully covered)
                                    if alignment_next['ref_end'] > alignment_current['ref_start']:
                                        tandem_duplications.append((ref_chr, alignment_next['ref_start'], alignment_next['ref_start'] + deviation, True, True))
                                    #Tandem duplication (not fully covered)
                                    elif distance_on_reference >= -options.max_sv_size:
                                        tandem_duplications.append((ref_chr, alignment_next['ref_start'], alignment_next['ref_start'] + deviation, False, True))
                                    #Either very large TANDEM or TRANS
                                    else:
                                        sv_candidates[types_to_output[bnd_num]].append([ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_start'], 'fwd', read_name])
                                        translocations.append(('fwd', 'fwd', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_start']))
                                else:
                                    #Tandem Duplication
                                    if alignment_next['ref_start'] < alignment_current['ref_end']:
                                        tandem_duplications.append((ref_chr, alignment_current['ref_start'], alignment_current['ref_start'] + deviation, True, False))
                                    #Large tandem duplication
                                    elif distance_on_reference >= -options.max_sv_size:
                                        tandem_duplications.append((ref_chr, alignment_current['ref_start'], alignment_current['ref_start'] + deviation, False, False))
                                    #Either very large TANDEM or TRANS
                                    else:
                                        sv_candidates[types_to_output[bnd_num]].append([ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_end'] - 1, 'rev', read_name])
                                        translocations.append(('rev', 'rev', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_end'] - 1))
            #Different orientations
            else:#处理inversion：前面是处理内部的is_reverse不直接标为iversion,而是标为ins或del
                #Normal to reverse
                #这些处理什么意思
                if not alignment_current['is_reverse'] and alignment_next['is_reverse']:
                    distance_on_reference = alignment_next['ref_end'] - alignment_current['ref_end']
                    deviation = distance_on_read - distance_on_reference
                    if -options.query_overlap_tolerance <= distance_on_read <= options.query_gap_tolerance:
                        if alignment_next['ref_start'] - alignment_current['ref_end'] >= -options.reference_overlap_tolerance: # Case 1
                            #INV candidate
                            if options.min_sv_size <= -deviation <= options.max_sv_size:
                                inversions.append((ref_chr, alignment_current['ref_end'], alignment_current['ref_end'] - deviation, "left_fwd"))
                                #transitions.append(('inversion', 'left_fwd', ref_chr, alignment_current['ref_end'], alignment_next['ref_end']))
                            #Either very large INV or TRANS 它将比较大的都给作为trans
                            else:
                                sv_candidates[types_to_output[bnd_num]].append([ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_end'] - 1, 'rev', read_name])
                                translocations.append(('fwd', 'rev', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_end'] - 1))
                        elif alignment_current['ref_start'] - alignment_next['ref_end'] >= -options.reference_overlap_tolerance: # Case 3
                            #INV candidate
                            if options.min_sv_size <= deviation <= options.max_sv_size:
                                inversions.append((ref_chr, alignment_next['ref_end'], alignment_next['ref_end'] + deviation, "left_rev"))
                                #transitions.append(('inversion', 'left_rev', ref_chr, alignment_next['ref_end'], alignment_current['ref_end']))
                            #Either very large INV or TRANS
                            else:
                                sv_candidates[types_to_output[bnd_num]].append([ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_end'] - 1, 'rev', read_name])
                                translocations.append(('fwd', 'rev', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_end'] - 1))
                    else:#他直接pass掉，可能是没有这样的信号？
                        pass
                        #print("Overlapping read segments in read", read_name)
                #Reverse to normal
                if alignment_current['is_reverse'] and not alignment_next['is_reverse']:
                    distance_on_reference = alignment_next['ref_start'] - alignment_current['ref_start'] 
                    deviation = distance_on_read - distance_on_reference
                    if -options.query_overlap_tolerance <= distance_on_read <= options.query_gap_tolerance:
                        if alignment_next['ref_start'] - alignment_current['ref_end'] >= -options.reference_overlap_tolerance: # Case 2
                            #INV candidate
                            if options.min_sv_size <= -deviation <= options.max_sv_size:
                                inversions.append((ref_chr, alignment_current['ref_start'], alignment_current['ref_start'] - deviation, "right_fwd"))
                                #transitions.append(('inversion', 'right_fwd', ref_chr, alignment_current['ref_start'], alignment_next['ref_start']))
                            #Either very large INV or TRANS
                            else:
                                sv_candidates[types_to_output[bnd_num]].append([ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_start'], 'fwd', read_name])
                                translocations.append(('rev', 'fwd', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_start']))
                        elif alignment_current['ref_start'] - alignment_next['ref_end'] >= -options.reference_overlap_tolerance: # Case 4
                            #INV candidate
                            if options.min_sv_size <= deviation <= options.max_sv_size:
                                inversions.append((ref_chr, alignment_next['ref_start'], alignment_next['ref_start'] + deviation, "right_rev"))
                                #transitions.append(('inversion', 'right_rev', ref_chr, alignment_next['ref_start'], alignment_current['ref_start']))
                            #Either very large INV or TRANS
                            else:
                                sv_candidates[types_to_output[bnd_num]].append([ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_start'], 'fwd', read_name])
                                translocations.append(('rev', 'fwd', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_start']))
                    else:
                        pass
                        #print("Overlapping read segments in read", read_name)
        #Different chromosomes
        else:#Trans
            ref_chr_current = bam.getrname(alignment_current['ref_id'])
            ref_chr_next = bam.getrname(alignment_next['ref_id'])
            #Same orientation
            if alignment_current['is_reverse'] == alignment_next['is_reverse']:
                #No overlap on read
                if distance_on_read >= -options.query_overlap_tolerance:
                    #No gap on read
                    if distance_on_read <= options.query_gap_tolerance:
                        if not alignment_current['is_reverse']:
                            sv_candidates[types_to_output[bnd_num]].append([ref_chr_current, alignment_current['ref_end'] - 1, 'fwd', ref_chr_next, alignment_next['ref_start'], 'fwd', read_name])
                            translocations.append(('fwd', 'fwd', ref_chr_current, alignment_current['ref_end'] - 1, ref_chr_next, alignment_next['ref_start']))
                        else:
                            sv_candidates[types_to_output[bnd_num]].append([ref_chr_current, alignment_current['ref_start'], 'rev', ref_chr_next, alignment_next['ref_end'] - 1, 'rev', read_name])
                            translocations.append(('rev', 'rev', ref_chr_current, alignment_current['ref_start'], ref_chr_next, alignment_next['ref_end'] - 1))
                #Overlap on read
                else:
                    pass
                    #print("Overlapping read segments in read", read_name)
            #Different orientation
            else:
                #No overlap on read
                if distance_on_read >= -options.query_overlap_tolerance:
                    #No gap on read
                    if distance_on_read <= options.query_gap_tolerance:
                        if not alignment_current['is_reverse']:
                            sv_candidates[types_to_output[bnd_num]].append([ref_chr_current, alignment_current['ref_end'] - 1, 'fwd', ref_chr_next, alignment_next['ref_end'] - 1, 'rev', read_name])
                            translocations.append(('fwd', 'rev', ref_chr_current, alignment_current['ref_end'] - 1, ref_chr_next, alignment_next['ref_end'] - 1))
                        else:
                            sv_candidates[types_to_output[bnd_num]].append([ref_chr_current, alignment_current['ref_start'], 'rev', ref_chr_next, alignment_next['ref_start'], 'fwd', read_name])
                            translocations.append(('rev', 'fwd', ref_chr_current, alignment_current['ref_start'], ref_chr_next, alignment_next['ref_start']))
                #Overlap on read
                else:
                    pass
                    #print("Overlapping read segments in read", read_name)
    # if len(inversions)>0:
    #     print(len(inversions))

    #Handle tandem duplications
    #处理串联重复，将类似is_simiar 串联重复信号进行合并 +1
    current_chromosome = None
    current_starts = []
    current_ends = []
    current_copy_number = 0
    current_fully_covered = []
    for tandem_duplication in tandem_duplications:
        if current_chromosome == None:
            current_chromosome = tandem_duplication[0]
            current_starts.append(tandem_duplication[1])
            current_ends.append(tandem_duplication[2])
            current_copy_number = 1
            current_fully_covered.append(tandem_duplication[3])
            current_direction = tandem_duplication[4]
        else:
            if is_similar(current_chromosome, mean(current_starts), mean(current_ends), tandem_duplication[0], tandem_duplication[1], tandem_duplication[2]) and current_direction == tandem_duplication[4]:
                current_starts.append(tandem_duplication[1])
                current_ends.append(tandem_duplication[2])
                current_copy_number += 1
                current_fully_covered.append(tandem_duplication[3])
            else:
                fully_covered = True if sum(current_fully_covered) else False
                sv_candidates[types_to_output[duptan_num]].append([current_chromosome, int(mean(current_starts)), int(mean(current_ends)), current_copy_number, fully_covered, read_name])
                current_chromosome = tandem_duplication[0]
                current_starts =[tandem_duplication[1]]
                current_ends =[tandem_duplication[2]]
                current_copy_number = 1
                current_fully_covered = [tandem_duplication[3]]
    if current_chromosome != None:
        fully_covered = True if sum(current_fully_covered) else False
        sv_candidates[types_to_output[duptan_num]].append([current_chromosome, int(mean(current_starts)), int(mean(current_ends)), current_copy_number, fully_covered, read_name])

    #Handle interspersed duplications
    for this_index in range(len(translocations)):
        this_dir1 = translocations[this_index][0]
        this_dir2 = translocations[this_index][1]
        this_chr1 = translocations[this_index][2]
        this_pos1 = translocations[this_index][3]
        this_chr2 = translocations[this_index][4]
        this_pos2 = translocations[this_index][5]

        for before_dir1, before_dir2, before_chr1, before_pos1, before_chr2, before_pos2 in translocations[:this_index]:
            #Same direction at destination and origin
            if before_dir1 == this_dir2 and before_dir2 == this_dir1:
                #Same position at destination
                if is_similar(before_chr1, before_pos1, 0, this_chr2, this_pos2, 0):
                    #Same chromosome for origin
                    if before_chr2 == this_chr1:
                        #INS_DUP candidate
                        if before_dir2 == before_dir1:
                            if before_dir1 == 'fwd':
                                length = this_pos1 + 1 - before_pos2
                                if options.min_sv_size <= length <= options.max_sv_size:
                                    sv_candidates[types_to_output[dupint_num]].append([before_chr2, before_pos2, this_pos1 + 1, before_chr1, int(mean([before_pos1 + 1, this_pos2])), int(mean([before_pos1 + 1, this_pos2])) + length, read_name])
                            elif before_dir1 == 'rev':
                                length = before_pos2 + 1 - this_pos1
                                if options.min_sv_size <= length <= options.max_sv_size:
                                    sv_candidates[types_to_output[dupint_num]].append([before_chr2, this_pos1, before_pos2 + 1, before_chr1, int(mean([before_pos1, this_pos2 + 1])), int(mean([before_pos1, this_pos2 + 1])) + length, read_name])
                        #INV_INS_DUP candidate
                        else:
                            pass
    #合并inversion信号 
    #Handle inversions (simple inversions produce two novel adjacencies that need to be merged for a complete candidate)
    sorted_inversions = sorted(inversions, key=lambda inversion: (inversion[0], inversion[1], inversion[2])) 
    active_inversions = []
    for inversion in sorted_inversions:
        chrom, start, end, direction = inversion
        if len(active_inversions) == 0:
            active_inversions.append(inversion)
        else:
            #If current inversion overlaps one of the active inversions
            if chrom == active_inversions[-1][0] and start < max([i[2] for i in active_inversions]):
                active_inversions.append(inversion)
            else:
                #Cluster inversions
                sv_candidates[types_to_output[inv_num]].extend(process_overlapping_inversions(active_inversions, read_name, bam))
                active_inversions = []
    if len(active_inversions) > 0:
        sv_candidates[types_to_output[inv_num]].extend(process_overlapping_inversions(active_inversions, read_name, bam))
    # sv_candidate['DEL'].extend(cluster_del(sv_candidates['DEL']))  
    # sv_candidate['INS'].extend(cluster_ins(sv_candidates['INS']))  

    # return sv_candidate
Max_ref_gap=500#500
Max_read_gap=500#500
Max_ref_overlap=3000#3000
Max_read_overlap=300#300#ins del
Max_dup_readoverlap=500
Min_sv_size=40
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
def aln_filter(new_alns):
    # 假设 new_alns 是输入的列表，newalns 是输出的列表
    newalns = []
    alni = 1
    new_alns.sort(key=lambda x:(x[1],x[2]))
    # plot(new_alns,'h2tg56lbefore')
    newalns.append(new_alns[0])
    while alni < len(new_alns):
        
        if new_alns[alni][2]-newalns[-1][2]<100:
            alni+=1
        else:
            newalns.append(new_alns[alni])
            alni+=1
    # plot(newalns,'h2tg56lafter')
    return newalns

def analysis_split_read2(primary, supplementaries,sv_candidates,  readname,header,options):
    bnd_candidates=[]
    inv_candidates=[]
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
        alignment_list.append([refename, alignment.query_alignment_start,alignment.query_alignment_end,alignment.reference_start,alignment.reference_end,alignment.is_reverse])
    # filter在内部的alignment等
    alignment_list=aln_filter(alignment_list)
    # alignment_list.sort(key=lambda x:(x[al_q_start],x[al_q_end]))
    for index in range(len(alignment_list) - 1):
        ele1 = alignment_list[index][:]
        ele2 = alignment_list[index + 1][:]
        if ele2[al_chr_name]==ele1[al_chr_name]:#chr
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
                        elif abs(delta_length)>Max_sv_Size:
                        #trans same#或者是直接为参考基因组上平常区域
                            dupinters.append( [ele1[al_is_reverse],ele2[al_is_reverse],ele1[al_chr_name], ele1[al_r_end],ele1[al_chr_name], ele2[al_r_start],readname,'del/ins'])
                            print(delta_length)
                        #存不存在一种 ref差距不大 dup大的情况 是否可以看作dup
                    else:#dup
                        if Max_read_gap>=dis_read>=-Max_dup_readoverlap:#无overlap on read 无gap
                            if dis_ref<=0:
                                if delta_length>=Min_sv_size:
                                    
                                    if ele1[al_r_start]<=ele2[al_r_start]:
                                        #case1 dup_complete
                                        if ele2[al_r_end]<=ele1[al_r_end]:
                                            duptans.append([duptan_num,ele1[al_chr_name],ele2[al_r_start],ele2[al_r_end],True,ele2[al_is_reverse],ele2[al_q_start],readname,1])
                                        else:#case2 dup_after
                                            duptans.append([duptan_num,ele1[al_chr_name],ele2[al_r_start],ele1[al_r_end],True,ele2[al_is_reverse],ele2[al_q_start],readname,2])
                                    else:
                                        #足够近标为dup;否则则是trans
                                        if ele1[al_r_start]-ele2[al_r_end]<Max_ref_gap:#case3 dup_before
                                            duptans.append([duptan_num,ele2[al_chr_name],ele2[al_r_start],ele2[al_r_end],True,ele2[al_is_reverse],ele2[al_q_start],readname,3])
                                        else:#认为可能发生了易位
                                           dupinters.append([ele1[al_is_reverse],ele2[al_is_reverse],ele1[al_chr_name],ele1[al_r_start],ele2[al_chr_name],ele2[al_r_start],readname,'dup_int_before'])

            else:#inv
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
                        else:
                            # print("this -++")
                            pass
                        continue         

                        
                    if abs(ele1[al_q_end]-ele2[al_q_start])<=100 and ele2[al_r_start]-ele1[al_r_start]<=Max_sv_Size:
                        if ele1[al_r_start]<ele2[al_r_start]:
                            inv_candidates.append([inv_num,ele1[al_chr_name],ele1[al_r_start],ele2[al_r_start],'-+',readname])
                        else:
                            print('inv',end='\t')
                            print([inv_num,ele1[al_chr_name],ele1[al_r_start],ele2[al_r_start],'-+',readname])
                    else:
                        dupinters.append([ele1[al_is_reverse],ele2[al_is_reverse],ele1[al_chr_name],ele1[al_r_start],ele2[al_chr_name],ele2[al_r_start],readname,'inv_-+'])

                    # print("-++")
                    # pass
                elif ele2[al_is_reverse]:#+-
                    if index+2<len(alignment_list):
                        ele3=alignment_list[index+2][:]
                        if not ele3[al_is_reverse] and ele3[al_chr_name]==ele1[al_chr_name]:
                            if ele2[al_q_start] + 0.5 * (ele3[al_r_start] - ele1[al_r_end]) >= ele1[al_q_end] and ele3[al_q_start] + 0.5 * (ele3[al_r_start] - ele1[al_r_end]) >= ele2[al_q_end]:#两边没有太大的gap
                                # No overlaps in split reads

                                if ele2[al_r_start]-ele1[al_r_end]>= -Min_sv_size and ele3[al_r_start] - ele2[al_r_end]>=-Min_sv_size:#保证比对也在中间
                                    inv_candidates.append([inv_num,ele1[al_chr_name],ele1[al_r_end],ele2[al_r_end],"+-+",readname])
                                    
                                   
                                else:#说明是移位反转
                                    bnd_candidates.append([bnd_num,ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_start],readname,'INV bnd +-+'])
                                    inv_candidates.append([inv_num,ele2[al_chr_name],ele2[al_r_start],ele2[al_r_end],"+-+",readname])
                        else:
                            # print("this +--")    
                            pass
                        continue
                    if abs(ele1[al_q_end]-ele2[al_q_start])<=100 and ele2[al_r_end]-ele1[al_r_end]<=Max_sv_Size:
                        if ele1[al_r_end] < ele2[al_r_end]:

                            inv_candidates.append([inv_num,ele1[al_chr_name],ele1[al_r_end],ele2[al_r_end],'+-',readname])
                        else:
                            print('inv',end='\t')
                            print([inv_num,ele1[al_chr_name],ele1[al_r_end],ele2[al_r_end],'+-',readname])
                    else:
                        dupinters.append([ele1[al_is_reverse],ele2[al_is_reverse],ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_end],readname,'inv_+-'])

                    # print("+--")
                    # pass  
                 # if ele1[al_is_reverse]:
                #     dis_ref = ele2[al_r_start] - ele1[al_r_start]
                #     dis_read = ele2[al_q_start] - ele1[al_q_end]
                #     delta_length = dis_read-dis_ref
                #     length=min(ele1[2]-ele1[1],(ele2[2]-ele2[1]))
                #     ratio=dis_read/length
                #     # print("ratio",ratio)
                #     # if  Max_read_gap>=dis_read>=-3000:
                #     if -0.2<=ratio<=0.8:
                #     # if  Max_read_gap>=dis_read>=-3000:
                #         if ele2[al_r_start] - ele1[al_r_end]>=-Max_ref_overlap: # Case 1 -+
                #             #INV candidate
                #             if Min_sv_size <= -delta_length <= Max_sv_Size:
                #                 inversions.append((ele1[al_chr_name], ele1[al_r_start], ele1[al_r_start] - delta_length,index, "-+"))
                                
                #             #Either very large INV or TRANS
                #             else:
                #                 inversions.append((ele1[al_chr_name], ele1[al_r_start], ele2[al_r_start], index,"--+"))
                
                #         else: # Case 2:trans inv -
                #             #INV candidate
                #             if Min_sv_size <= -delta_length <= Max_sv_Size:
                #                 inversions.append((ele1[al_chr_name], ele2[al_r_start], ele2[al_r_start] + delta_length, index,"-+ trans"))
                #             #Either very large INV or TRANS
                #             else:
                #                 inversions.append((ele1[al_chr_name], ele1[al_r_start], ele2[al_r_start],index, "--+ trans"))
                #     # else:#trans inv
                #     #     print("Overlapping read segments in read",dis_read)
                
                # elif ele2[al_is_reverse]:
                #     dis_ref = ele2[al_r_end] - ele1[al_r_end]
                #     dis_read = ele2[al_q_start] - ele1[al_q_end]
                #     length=min(ele1[2]-ele1[1],(ele2[2]-ele2[1]))
                #     ratio=dis_read/length
                #     # print("ratio",ratio)
                #     delta_length = dis_read-dis_ref
                #     # if  Max_read_gap>=dis_read>=-3000:
                #     if -0.2<=ratio<=0.8:
                #         if ele2[al_r_start] - ele1[al_r_end] >= -Max_ref_overlap: # Case 1 -+
                #             #INV candidate
                #             if Min_sv_size <= -delta_length <= Max_sv_Size:
                #                 inversions.append((ele1[al_chr_name], ele1[al_r_end], ele1[al_r_end] - delta_length,index, "+-"))
                                
                #             #Either very large INV or TRANS
                #             else:
                #                 inversions.append((ele1[al_chr_name], ele1[al_r_end], ele2[al_r_end],index, "+--"))
                
                #         else: # Case 2:trans inv -
                #             #INV candidate
                #             if Min_sv_size <= -delta_length <= Max_sv_Size:
                #                 inversions.append((ele1[al_chr_name], ele2[al_r_end], ele2[al_r_end] + delta_length,index, "+- trans"))
                #             #Either very large INV or TRANS
                #             else:
                #                 inversions.append((ele1[al_chr_name], ele1[al_r_end], ele2[al_r_end],index, "+-- trans"))
                #     # else:
                #     #     print("Overlapping read segments in read",dis_read)
                    
        else:#dif chr trans
            if abs(ele1[al_q_end]-ele2[al_q_start])<=1000:
                bnd_candidates.append([bnd_num,ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_start],readname,'Nor'])
            elif  abs(ele2[al_q_end]-ele1[al_q_start])<=1000:
                bnd_candidates.append([bnd_num,ele1[al_chr_name],ele1[al_r_start],ele2[al_chr_name],ele2[al_r_end],readname,'Nor'])
    #汇入总的
    # for can in candidates:
    #     num=can[0]
    #     sv_candidates[types_to_output[num]].append(can)
    # sv_candidates[types_to_output[dupint_num]].extend(dupinters)
    # sv_candidates[types_to_output[duptan_num]].extend(duptans)
    # sv_candidates[types_to_output[bnd_num]].extend(translocations)
    #打印输出
    # if len(dupinters)!=0 or len(duptans)!=0 :
    #     print('dupinter'+str(len(dupinters)),end='\t')
    #     plot(alignment_list,readname)
    #     print('duptans'+str(len(duptans)),end='\t')
    #     print('/n')  
    #对dup inter 区域进行潜在的聚类 因为会有可能断点区域；不是的话就为bnd
    #生成的duprtan列表需要怎么处理  太多了
    #两部聚类 减少相同的？
    #提速
    #实验###
    #我把del/ins中关于复杂区域中 在inter中的 识别给过滤掉了
    #cluster
    
    sv_candidates[types_to_output[del_num]].extend(cluster_del(del_candidates))
    sv_candidates[types_to_output[ins_num]].extend(cluster_ins(ins_candidates))
    can1,can2=cluster_dupint(dupinters)
    sv_candidates[types_to_output[dupint_num]].extend(can1)
    sv_candidates[types_to_output[duptan_num]].extend(cluster_dup(duptans))
    # sv_candidates[types_to_output[bnd_num]].extend(can2)
    sv_candidates[types_to_output[bnd_num]].extend(bnd_candidates)
    sv_candidates[types_to_output[inv_num]].extend(inv_candidates)