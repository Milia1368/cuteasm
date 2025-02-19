
import numpy as np
from collections import Counter
from statistics import mean

max_cluster_bias=1000
diff_ratio_merging_INS=0.9
max_cluster_bias_DEL=1000
diff_ratio_merging_DEL=0.5


Max_ref_gap=500
Max_read_gap=500
Max_ref_overlap=3000
Max_read_overlap=300#ins del
Max_dup_readoverlap=500
Min_sv_size=40
Max_sv_Size=10000
al_r_start=3



#2 start 3 length
def cluster_del(sig_list, max_shift=20, min_overlap_ratio=0.5, min_size_similarity=0.5):#将有重叠的选择最长的
    # 使用 NumPy 数组初始化聚类数组
    cluster = -1 * np.ones(len(sig_list), dtype=int)  # 初始化聚类数组为 -1
    sig_list=sorted(sig_list,key=lambda x:(x[2],x[3]))
    def calculate_overlap_ratio(start1, end1, start2, end2):
        """计算重叠比例"""
        overlap_length = max(0, min(end1, end2) - max(start1, start2))
        return overlap_length / min(end1 - start1, end2 - start2)

    while np.any(cluster == -1):  # 当还有未分配的聚类时
        for i, sig1 in enumerate(sig_list):
            if cluster[i] == -1:
                cluster[i] = i  # 将当前信号标记为其自身的聚类
                start1, end1 = int(sig1[2]), int(sig1[2]) + int(sig1[3])

                for j in range(i + 1, len(sig_list)):
                    sig2 = sig_list[j]
                    if cluster[j] == -1 and sig1[1]==sig2[1]:
                        start2, end2 = int(sig2[2]), int(sig2[2]) + int(sig2[3])
                        shift = abs(start1 - start2)
                        # if shift>max_cluster_bias:
                        overlap_ratio = calculate_overlap_ratio(start1, end1, start2, end2)
                        size_similarity = min(sig1[3], sig2[3]) / max(sig1[3], sig2[3])

                        if shift <= max_shift and overlap_ratio >= min_overlap_ratio and size_similarity >= min_size_similarity:
                            cluster[j] = cluster[i]  # 将 sig2 分配到 sig1 的聚类

    # 统计有效聚类
    valid_clusters = [idx for idx, count in Counter(cluster).items() if count >= 1]
    


    final_sig_list = []
    # 假设这里的valid_clusters、cluster、sig_list等相关变量已经有合适的值（示例中未完整给出其赋值情况）
    for cluster_idx in valid_clusters:
        sig_idxs = np.where(cluster == cluster_idx)[0]
        if len(sig_idxs) > 0:
            # 提取对应索引位置的信号列表
            selected_sigs = [sig_list[idx] for idx in sig_idxs]
            
            # 初始化用于计算平均值的变量（根据信号元素的结构和需要平均的属性来确定）
            avg_sig = [0] * len(selected_sigs[0])
            for i in range(len(selected_sigs[0])):
                if i in [2, 3]:  # 这里假设对索引为2和3的元素求平均，你可按需调整
                    avg_sig[i] =int(np.mean([sig[i] for sig in selected_sigs]))
                else:
                    # 对于其他元素，可选择取第一个信号的对应元素值（也可根据实际需求改变处理方式）
                    avg_sig[i] = selected_sigs[0][i]
            final_sig_list.append(avg_sig)
    # if len(final_sig_list)<len(sig_list):
    #     for sig in sig_list:
    #         if sig not in final_sig_list:
    #             print('del',end='\t')
    #             print(sig)

    return final_sig_list
def cluster_ins(sig_list, max_shift=20, min_size_similarity=0.5):#将有重叠的选择最长的
    # 使用 NumPy 数组初始化聚类数组
    cluster = -1 * np.ones(len(sig_list), dtype=int)  # 初始化聚类数组为 -1

    while np.any(cluster == -1):  # 当还有未分配的聚类时
        for i, sig1 in enumerate(sig_list):
            if cluster[i]==-1:
                sig1 = sig_list[i]
                cluster[i] = i
                for j in range(i + 1, len(sig_list)):
                    if cluster[j]==-1:
                        sig2 = sig_list[j]
                        if sig1[1]!=sig2[1]:
                            continue
                        size_similarity = min(sig1[3],sig2[3])/max(sig1[3],sig2[3])
                        shift = abs(sig1[2]-sig2[2])
                        
                        if (shift<=max_shift) and\
                        (size_similarity >= min_size_similarity):
                            cluster[j]=cluster[i]    
    # 统计有效聚类
    valid_clusters = [idx for idx, count in Counter(cluster).items() if count >= 1]
    
    # 选择每个有效聚类中的最佳信号
    final_sig_list = []
    for cluster_idx in valid_clusters:
        sig_idxs = np.where(cluster == cluster_idx)[0]
        best_sig = sig_list[sig_idxs[0]]
        
        # 使用 NumPy 的向量化操作来找到最佳信号
        best_sig = max((sig_list[idx] for idx in sig_idxs), key=lambda x: x[3])  # 选择长度最大的信号
                
        final_sig_list.append(best_sig)
    # if len(final_sig_list)<len(sig_list):
    #     for sig in sig_list:
    #         if sig not in final_sig_list:
    #             print('ins',end='\t')
    #             print(sig)

    return final_sig_list



def is_similar(chr1, start1, end1, chr2, start2, end2):
    if chr1 == chr2 and abs(start1 - start2) < 200 and abs(end1 - end2) < 200:
        return True
    else:
        return False
def cluster_dupint(translocations):
    # 先按照位置3进行排序
    sorted_translocations = sorted(translocations, key=lambda x: x[3])

    # 第一步，按照位置3距离小于1000进行初步聚类，并取每个聚类的平均值作为代表
    distance_threshold_1 = 1000
    position_clusters_1 = []
    for signal in sorted_translocations:
        assigned = False
        for cluster in position_clusters_1:
            first_signal = cluster[0]
            if abs(signal[3] - first_signal[3]) < distance_threshold_1:
                cluster.append(signal)
                assigned = True
                break
        if not assigned:
            position_clusters_1.append([signal])

    # 计算每个初步聚类的平均值作为代表（根据实际信号元素个数和需要平均的属性调整）
    averaged_clusters_1 = []
    for cluster in position_clusters_1:
        num_signals = len(cluster)
        avg_signal = [0] * len(cluster[0])  # 根据信号元素个数初始化用于平均计算的元素列表
        for i in range(len(cluster[0])):
            if i in [3, 5]:  # 对位置3和位置5对应的属性进行平均计算
                avg_signal[i] = int(sum(signal[i] for signal in cluster) / num_signals)
            else:
                avg_signal[i] = cluster[0][i]  # 其他属性直接取第一个信号的对应属性（可按需调整）
        averaged_clusters_1.append(avg_signal)

    # 第二步，基于位置5以10000为距离再次聚类
    distance_threshold_2 = 10000
    position_clusters_2 = []
    for signal in averaged_clusters_1:
        assigned = False
        for cluster in position_clusters_2:
            first_signal = cluster[0]
            if abs(signal[5] - first_signal[5]) < distance_threshold_2:
                cluster.append(signal)
                assigned = True
                break
        if not assigned:
            position_clusters_2.append([signal])

    # 第三步，根据聚类簇大小判断类型并返回结果
    dup_int_clusters = []
    trans_clusters = []
    for cluster in position_clusters_2:
        if len(cluster) > 2:
            dup_int_clusters.extend(cluster)
        else:
            trans_clusters.extend(cluster)
    trans_clusters=[[5]+tan[2:] for tan in trans_clusters]
    dup_int_clusters=[[4]+tan[2:] for tan in dup_int_clusters]
    return dup_int_clusters, trans_clusters


    
def cluster_dup(tandem_duplications):
    duptan_can=[]
    tandem_duplications=sorted(tandem_duplications,key=lambda inversion: (inversion[3], inversion[5]))
    current_chromosome = None
    current_starts = []
    current_ends = []
    current_copy_number = 0
    current_fully_covered = []
    for tandem_duplication in tandem_duplications:
        if current_chromosome == None:
            current_chromosome = tandem_duplication[1]
            current_starts.append(tandem_duplication[2])
            current_ends.append(tandem_duplication[3])
            current_copy_number = 1
            current_fully_covered.append(tandem_duplication[4])
            current_direction = tandem_duplication[5]
            read_name=tandem_duplication[7]
        else:
            if is_similar(current_chromosome, mean(current_starts), mean(current_ends), tandem_duplication[1], tandem_duplication[2], tandem_duplication[3]) and current_direction == tandem_duplication[5]:
                current_starts.append(tandem_duplication[2])
                current_ends.append(tandem_duplication[3])
                current_copy_number += 1
                current_fully_covered.append(tandem_duplication[4])
            else:
                fully_covered = True if sum(current_fully_covered) else False
                duptan_can.append([current_chromosome, int(mean(current_starts)), int(mean(current_ends)), current_copy_number, fully_covered, read_name])
                current_chromosome = tandem_duplication[1]
                current_starts =[tandem_duplication[2]]
                current_ends =[tandem_duplication[3]]
                current_copy_number = 1
                current_fully_covered = [tandem_duplication[4]]
                current_direction = tandem_duplication[5]
    if current_chromosome != None:
        fully_covered = True if sum(current_fully_covered) else False
        duptan_can.append([current_chromosome, int(mean(current_starts)), int(mean(current_ends)), current_copy_number, fully_covered, read_name])
    # print('duptan',end='\t')
    # print(len(tandem_duplications),end='\t')
    # print(len(duptan_can))
    return duptan_can
#inv去重
#[inv_num,chr,start,end,'invs',readname]
def remove_invs(inversions):
    inversions=sorted(inversions,key=lambda inversion: (inversion[1],inversion[2], inversion[3]))#sort
    invs=[]
    if len(inversions)==0:
        return invs
    invs.append(inversions[0])
    for inv in inversions[1:]:
        if inv[1]==invs[-1][1]:#chr
            if inv[2]==invs[-1][2] and inv[3]==invs[-1][3]:#start end
                # print('inv',end='\t')
                # print(inv)
                continue
        invs.append(inv)
    return invs

#bnd merge
def merge_bnd(bndlist):
    #bnd里有一部分是已经匹配成果的bnd，一部分是
    invs=[]
    nor=[]
    same_bnd={}
    same_bnd[0]=[]
    same_bnd[1]=[]
    # bnd=[]
    bnd_can=[]
    for bndcan in bndlist:
        current_type=bndcan[-2]
        if current_type==3:#invs
            invs.append(bndcan)
        elif current_type==2:
            nor.append(bndcan)
        else:
            same_bnd[bndcan[0]].append(bndcan)
    invs.sort(key=lambda x:(x[1],x[2]))
    
    same_bnd[0].sort(key=lambda x:(x[2],x[3]))
    same_bnd[1].sort(key=lambda x:(x[2],x[3]))
    #匹配 invs和nor\same匹配
    ##    same_bnd和自身匹配
    for bnd0 in same_bnd[0]:
        chrom=bnd0[2]
        start=bnd0[3]
        length=bnd0[-1]
        for bnd1 in same_bnd[1]:
            chrom1=bnd1[2]
            start1=bnd1[3]#模糊的开始
            length1=bnd1[-1]#重点
            if chrom1==chrom:
                if abs(length1-length)<=5:
                    if bnd1[0]!=-1:
                        bnd_can.append([1,chrom,min(start,bnd1[6]),chrom,max(start,bnd1[6]),bnd0[1],max(length,length1),False,False])
                        bnd1[0]=-1
    nor_chr={}#对不同chr的bnd匹·配
    nor.sort(key=lambda x:(x[1],x[2]))
    for norcan in  nor:
        if norcan[1] not in nor_chr:
            nor_chr[norcan[1]]=[]
        nor_chr[norcan[1]].append(norcan)
    for key in nor_chr.keys():
        nor_chr[key].sort(key=lambda x:(x[1],x[2]))
    for norcan in nor:
        if norcan[0]!=-1:
            norcan[0]=-1
            chrom=norcan[1]
            anchorm=norcan[4]
            length=norcan[-1]
            if anchorm in nor_chr:
                for can in nor_chr[anchorm]:
                    if can[0]!=-1:
                    
                        length2=can[-1]
                        if abs(length2-length)<=5:
                            can[0]=-1
                            bnd_can.append([0,chrom,norcan[2],anchorm,can[2],norcan[-3],max(length,length2),False,False])
        
    for norcan  in  bnd_can:
        chrom0=norcan[1]
        start0=norcan[2]
        chrom1=norcan[3]
        start1=norcan[4]
        length0=norcan[-3]
        for invcan in invs:
            if invcan[0]!=-1:
                chrom=invcan[1]
                start=invcan[2]
                length=invcan[-1]
                if chrom==chrom0:
                    if abs(length0-length)<=5:
                        if abs(start0-start)<=5:
                            norcan[-2]=True
                            invcan[0]=-1
                elif chrom1==chrom:
                        if abs(length0-length)<=5:
                            if abs(start-start1)<=5:
                                norcan[-1]=True
                                invcan[0]=-1
    return bnd_can
def cluster_bnd(bnd_can):#去重[0,chrom1,start1,chorm2,start2,readname,length,False,False]
    seen = set()            #[5,ele1[al_chr_name],ele1[al_r_start],ele2[al_chr_name],ele2[al_r_end],readname,length,False,False]
    result = []
    for item in bnd_can:
        tpl = tuple(item[1:5])  # 把除了最后两个元素之外的部分转为元组，用于去重判断
        if tpl not in seen:
            seen.add(tpl)
            result.append(item)
    return result
def cluster_SV(svlist):#cuteSV的重实现
    svtype=svlist[0][0]
    final_cluster=[]
    if svtype=='DEL':
        return cluster_ins_in_hap(svlist,final_cluster)
    # elif svtype=='INS':
    #     return SR_del_in_hap(svlist,final_cluster)
    elif svtype=='DUP':
        return cluster_dup(svlist)
    else:
        return svlist
    pass

def cluster_ins_in_hap(ins_list,final_cluster):
    ins_list.sort(key=lambda x:(x[2],x[3]))
    temp_cluster=list()
    temp_cluster.append([])
    for sig in ins_list:
        if len(temp_cluster[-1])==0:
            temp_cluster[-1].append(sig)
        else:
            if not SR_ins(temp_cluster[-1],sig):
                final_cluster.append(temp_cluster)
                temp_cluster=list()
                temp_cluster[-1].append(sig)
            else:
                temp_cluster.append(sig)
    final_cluster.append(temp_cluster)
    return final_cluster
def SR_ins(sig1,sig2):
    if sig1[1]==sig2[1]:
        if abs(sig1[2]-sig2[2])<50:
            if min(sig1[3],sig2[3])/max(sig1[3],sig2[3])>0.8:
                return True
    return False
def SR_deldupinv(sig1,sig2):
    if sig1[1]==sig2[1]:
        # diffr=1 if abs(0.5*(sig1[1]+sig1[2])-) else 0
        if abs(sig1[2]-sig2[2])<50:
            if min(sig1[3],sig2[3])/max(sig1[3],sig2[3])>0.8:
                return True
    return False
def SR_ins(sig1,sig2):
    if sig1[1]==sig2[1]:
        if abs(sig1[2]-sig2[2])<50:
            if min(sig1[3],sig2[3])/max(sig1[3],sig2[3])>0.8:
                return True
    return False

if __name__ == "__main__":
    

    # 定义示例信号列表
    sig_list = [
        ['chr1', 100, 150, 50, 'info1'],  # 信号 1
        ['chr1', 120, 170, 50, 'info2'],  # 信号 2 (与信号 1 有重叠)
        ['chr1', 200, 250, 50, 'info3'],  # 信号 3
        ['chr1', 300, 350, 50, 'info4'],  # 信号 4
        ['chr1', 400, 450, 50, 'info5'],  # 信号 5
    ]

    # 调用 cluster_del 函数
    final_signals = cluster_del(sig_list, max_shift=50, min_overlap_ratio=0.5, min_size_similarity=0.5)

    # 输出结果
    print("最终选择的信号:")
    for sig in final_signals:
        print(sig)