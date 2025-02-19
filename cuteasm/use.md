# del ins 
# 对cigar里相近的del/ins信号聚类 合并成一个
**https://github.com/maiziezhoulab/FocalSV/blob/main/focalsv/4_sv_calling/Dippav/extract_contig_signature.py#L196**
def calculate_overlap_ratio(start1,end1,start2,end2):
        len1 = end1-start1
        len2 = end2-start2
        minlen = min(len1,len2)
        return (min(end1,end2)-max(start1,start2))/minlen
    
    while (np.array(cluster)==-1).sum()!=0:
        for i in range(len(sig_list)):
            if cluster[i]==-1:
                sig1 = sig_list[i]
                cluster[i] = i
                for j in range(len(sig_list)):
                    if cluster[j]==-1:
                        sig2 = sig_list[j]
                        start1 = sig1[2]
                        start2 = sig2[2]
                        end1 = start1+sig1[3]
                        end2 = start2+sig2[3]
                        overlap_ratio = calculate_overlap_ratio(start1,end1,start2,end2)
                        size_similarity = min(sig1[3],sig2[3])/max(sig1[3],sig2[3])
                        shift = abs(sig1[2]-sig2[2])
                        
                        if (shift<=max_shift) and\
                        (overlap_ratio >= min_overlap_ratio) and\
                        (size_similarity >= min_size_similarity):
                            cluster[j]=cluster[i]
#    ## 并将同一个类中的最长的信号作为最好信号
        final_sig_list = []
    for cluster_idx in valid_cluster:
        sig_idxs = np.where(cluster==cluster_idx)

        best_sig = sig_list[sig_idxs[0][0]]
        for idx in sig_idxs[0]:
            new_sig = sig_list[idx]
            best_sig_len = best_sig[3]
            new_sig_len = new_sig[3]
            # update best sig
            if new_sig_len > best_sig_len:
                best_sig = new_sig
        final_sig_list.append(best_sig)
# split
# 只考虑两个方向相同
        Diffdis = (Ref2s-Ref1e)-(Read2s-Read1e)
        Diffolp = Ref1e - Ref2s
        
        if reverse1:
            direction = '-'
        else:
            direction = '+'
            


        if abs(Diffdis)<=max_svlen:
            if (Diffolp<30)  and (Diffdis >= 30 ):
                sigdel = [read1.reference_name,'DEL', Ref1e, Diffdis,read1.qname, Read1e,Read2s,direction,'split-alignment',"%d-%d"%(read1.mapq,read2.mapq)]
                del_list.append(sigdel)
            elif (Diffolp<3000)  and (Diffdis >= 30 ):
                sigdel = [read1.reference_name,'DEL', Ref1e-Diffdis, Diffdis,read1.qname, Read1e-Diffdis,Read2s-Diffdis,direction,'split-alignment',"%d-%d"%(read1.mapq,read2.mapq)]
                del_list.append(sigdel)
            elif (Diffolp<3000) and (Diffdis<= -30) :
                svlen = abs(Read2s-Read1e+Diffolp)
                if abs(Diffolp) > 400:
                    pos_ref = int((Ref1e+Ref2s)/2)
                else:
                    pos_ref = Ref2s
                sigins = [read1.reference_name,'INS', pos_ref, svlen , read1.qname,Read1e-Diffolp, Read2s,direction,'split-alignment',"%d-%d"%(read1.mapq,read2.mapq)]

                ins_list.append(sigins)
#   #再对信号进行聚类 上边的split 和 cigar 里的信号 ；一样的方法

#   pair del和ins 阈值设置
def pair_sig(sig_hp1,sig_hp2,max_compare_dist,max_shift, min_overlap_ratio,min_size_similarity):
    pair_status_hp1 = [-1]*len(sig_hp1)
    pair_status_hp2 = [-1]*len(sig_hp2)
    
    for i in range(len(sig_hp1)):
        sig1 = sig_hp1[i]
        for j in range(len(sig_hp2)):
            sig2 = sig_hp2[j]
            dist = sig2[2]-sig1[2]
            if dist>max_compare_dist:
                break
            elif (sig1[:2]==sig2[:2]) and (pair_status_hp2[j]==-1):
                if sig1[1]=='DEL':
                    result = pair_del(sig1,sig2,200,0.5,0.5)
                else:
                    result = pair_ins(sig1,sig2,200,0.5)
                if result==1:
                    pair_status_hp1[i] = j
                    pair_status_hp2[j] = i
                    break
#   1到k的聚类 去重
def get_size_sim(svlen1, svlen2):
    """计算两个长度的相似度"""
    return min(abs(svlen1), abs(svlen2)) / max(abs(svlen1), abs(svlen2))

def calculate_svlen(sig):
    """计算信号的长度差"""
    return abs(len(sig[3]) - len(sig[4]))

def match_one_pair(sig1, sig2, dist_thresh, size_sim_thresh, overlap_thresh=None, seq_sim_thresh=None):
    """通用匹配函数，支持插入和删除匹配"""
    dist_ref = abs(sig2[1] - sig1[1])
    svlen1, svlen2 = calculate_svlen(sig1), calculate_svlen(sig2)
    size_sim = get_size_sim(svlen1, svlen2)

    if dist_ref <= dist_thresh and size_sim >= size_sim_thresh:
        if overlap_thresh is not None:
            overlap = get_reciprocal_overlap(sig1, sig2)
            return overlap >= overlap_thresh
        elif seq_sim_thresh is not None:
            seq_sim = edit_sim(sig1[4], sig2[4])
            return seq_sim >= seq_sim_thresh
    return False

def get_reciprocal_overlap(sig1, sig2):
    """计算两个信号的重叠比例"""
    svlen1, svlen2 = calculate_svlen(sig1), calculate_svlen(sig2)
    start1, start2 = sig1[1], sig2[1]
    end1, end2 = start1 + svlen1, start2 + svlen2
    overlap = (min(end1, end2) - max(start1, start2)) / max(svlen1, svlen2)
    return overlap

def match_del_chr(sig_list, dist_thresh, size_sim_thresh, overlap_thresh):
    """在染色体内匹配信号"""
    links = []
    for sig1 in sig_list:
        pos1 = sig1[1]
        window = (pos1 - dist_thresh, pos1 + dist_thresh)
        comp_sig_list = [sig2 for sig2 in sig_list if sig2[1] <= window[1 and sig2 != sig1 and window[0] <= sig2[1]]

        for sig2 in comp_sig_list:
            if match_one_pair(sig1, sig2, dist_thresh, size_sim_thresh, overlap_thresh=overlap_thresh):
                links.append((sig1[2], sig2[2]))
    return links

def match_del(sig_list, dist_thresh, size_sim_thresh, overlap_thresh):
    """匹配所有染色体的信号"""
    links = []
    for i in tqdm(range(1, 23)):
        chr_name = f'chr{i}'
        sig_list_chr = [sig for sig in sig_list if sig[0] == chr_name]
        links.extend(match_del_chr(sig_list_chr, dist_thresh, size_sim_thresh, overlap_thresh))

    G = nx.Graph()
    G.add_edges_from(links)
    components = nx.connected_components(G)
    return [nodes for nodes in components] 
def pick_best_sv_one_cluster(vcf_dc, index_list):
    """从给定的索引列表中选择最佳结构变异（SV）索引"""
    sig_list = [vcf_dc[idx] for idx in index_list]  # 提取信号
    length_diffs = [abs(len(sig[3]) - len(sig[4])) for sig in sig_list]  # 计算长度差
    best_index = index_list[length_diffs.index(max(length_diffs))]  # 找到最大长度差的索引
    return best_index

def pick_best_sv(vcf_dc, nodes_list):
    """选择每个聚类中的最佳结构变异（SV）并返回保留和移除的索引"""
    retain_index = {}
    remove_index = {}

    for i, node_indices in enumerate(nodes_list):  # 遍历每个聚类
        best_index = pick_best_sv_one_cluster(vcf_dc, node_indices)  # 找到最佳索引
        retain_index[best_index] = i  # 记录保留的最佳索引

        # 记录移除的索引
        for idx in node_indices:
            if idx != best_index:
                remove_index[idx] = i

    return retain_index, remove_index  # 返回保留和移除的索引字典
def remove_redundancy(vcf_path,output_dir,dist_thresh  =500,
                      dist_thresh_del = 3000,
                      overlap_thresh = 0,
                      size_sim_thresh = 0.5,
                      size_sim_thresh_del = 0.1,
                      seq_sim_thresh = 0.5,
                      ):
    
    os.system("mkdir -p "+output_dir)
    ## set logger
    logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
    global logger
    logger = logging.getLogger(" ")
            
                
    # vcf_path = "dippav_variant_filtered.vcf"
    # dist_thresh, size_sim_thresh, seq_sim_thresh = 500,0.5,0.5
    # overlap_thresh = 0 
    # dist_thresh_del = 1000
    # output_dir = './'
    prefix = 'dippav_variant'
    global vcf_dc
    del_sig,ins_sig,vcf_dc,header = vcf_to_sig(vcf_path)
    links_del = match_del(del_sig,dist_thresh_del, size_sim_thresh_del, overlap_thresh)
    links_ins = match_ins(ins_sig,dist_thresh, size_sim_thresh, seq_sim_thresh)
    retain_index_del,remove_index_del = pick_best_sv(vcf_dc,links_del)
    retain_index_ins,remove_index_ins = pick_best_sv(vcf_dc,links_ins)
    write_vcf(output_dir,prefix,header,retain_index_del,remove_index_del,retain_index_ins,remove_index_ins)
    return      



# inter 提信号 my
#   ## v1.0
def analysis_split_read1(split_read,pri_chr,is_reverse,readname,query_length,query,bam):
    sv_candidate=[]
    same_chr_reads={}
    #按照start排序并编号
    split_read.sort(key=lambda x:x[1])
    count=0
    for read in split_read:
        read[-1]=count
        count+=1
        if read[7]  not in same_chr_reads:
            same_chr_reads[read[7] ]=[read]
        else:
            same_chr_reads[read[7]].append(read)
        
    for key in same_chr_reads:#按照每个染色体的起始位置排序
        same_chr_reads[key].sort(key=lambda x:(x[4],x[1]))
        #print_pos_sort(same_chr_reads[key])
    # for read in split_read:
    #     if read[-2]:
    #         read_s=query_length-read[2]
    #         read_e=query_length-read[1]
    #         read[1]=read_s
    #         read[2]=read_e
    #分析
    #主pri
    list_chr=same_chr_reads[pri_chr]
    reverse_num=[item[6] for item in split_read]
    reverse_flag=detect_reverse(reverse_num,is_reverse)
    sequnece_num= [item[-1] for item in list_chr]
    backbone=max_antichain(sequnece_num)
    same_chr_reads = [i for i in sequnece_num if i not in backbone]
    for i in range(len(split_read)):
        ai=split_read[i]
        #[aln1,query_start(aln1),query_start(aln1)+aln1.query_alignment_length,mapq,aln1.reference_start,r_end,is_reverse,0]
        if ai[-1] in backbone:
            if i+1<len(split_read):
                if split_read[i+1][-1] in backbone:
                    #相邻 暗示del/ins或是dup
                   sv_candidate.extend( analaysi_del_ins_dup(ai,split_read[i+1],readname,query,bam))
            
        elif ai[-1] in same_chr_reads:
            k=i-1
            while k>=0 and k not in backbone:
                k-=1
            if k>=0:
                sv_candidate.extend(analysis_dup_sametrans(split_read[k],ai,readname,bam))
        else:
            k=i-1
            while k>=0 and k not in backbone:
                k-=1
            if k>=0:
                sv_candidate.extend(analysis_trans(split_read[k],ai,readname,bam))
#inv
    active_inv=[]
    for i in reverse_flag:
       
        if i and len(active_inv):
            active_inv.append(i)
        else:
            if not i:
                a2=split_read[i]
                a1=split_read[i-1] 
                ele_1=[a1[1],a1[2],a1[4],a1[5],a1[7],a1[6],a1[0]]
                ele_2=[a2[1],a2[2],a2[4],a2[5],a2[7],a2[6],a2[0]]
                analysis_inv(ele_1,ele_2,[readname],sv_candidate,30,bam)
                

            
#对candidate聚类 复杂结构变异
   
    return sv_candidate
#   ##v2.0
def analysis_split_read2(split_read,pri_chr,is_reverse,readname,query_length,query,bam,candidates):
  
    same_chr_reads={}
    #按照start排序并编号
    split_read.sort(key=lambda x:x[1])
    #分配每个染色体上的aln
    count=0
    for read in split_read:
        read[8]=count
        read.append(-1)
        count+=1
        if read[7]  not in same_chr_reads:
            same_chr_reads[read[7] ]=[read]
        else:
            same_chr_reads[read[7]].append(read)
        
    # for key in same_chr_reads:#按照每个染色体的起始位置排序
    #     same_chr_reads[key].sort(key=lambda x:(x[4],x[1]))
    #对主要的染色体区域进行分块
    ref_blocks=[]
    new_ref_blocks=[]
    # pri_reads_num=[]
    try:
        for read in same_chr_reads[pri_chr]:
            ref_blocks.append([read[4],read[5],0,True,[read[8]]])
    except:
        max_value=0
        for key, value_list in same_chr_reads.items():
            # 找到当前列表的最大值
            current_max = len(value_list)
            
            # 如果当前最大值大于已知最大值，则更新
            if current_max > max_value:
                max_value = current_max
                max_key = key
        pri_chr=max_key
        for read in same_chr_reads[pri_chr]:
            ref_blocks.append([read[4],read[5],0,True,[read[8]]])
        # pri_reads_num.append(read[8])
    #对ref_blocks聚类 因为存在相近的blocks
    #ref_start  ref_end ref_num has_aln aln_num_list
    #   0           1       2       3       (4)
    ref_blocks.sort(key=lambda x :(x[1],x[2]))
    new_ref_blocks.append(ref_blocks[0])
    for i in range(len(ref_blocks)-1):#尽可能短的blocks 但如果前后overlap过大，则认为是一个blocks
        if ref_blocks[i+1][0]-new_ref_blocks[-1][1]<-Max_overlap_tolerance:#overlap足够大
            new_ref_blocks[-1]=[new_ref_blocks[-1][0],max(ref_blocks[i+1][1],new_ref_blocks[-1][1]),0,True,new_ref_blocks[-1][-1]+ref_blocks[i+1][-1]]
        else:
            new_ref_blocks.append(ref_blocks[i+1])
    # for i in range(len(ref_blocks)-1):
    #     if ref_blocks[i+1][0]-ref_blocks[i][1]>(Max_gap_tolerance+Max_gap_tolerance):#足够远的block之间添加新block
    #         new_ref_blocks.append([ref_blocks[i][1],ref_blocks[i+1][0],0,False])
    del ref_blocks
    #对 ref_block 进行编号并对对应的read编号
    num=1
    inter_dup=[]
    tandem_dup=[]
    dup_region=[]
    for block in new_ref_blocks:
        block[2]=num
        
        if block[3]:
            for aln in block[4]:
                split_read[aln][9]=num
        # if len(block[4])>=2:#极可能代表里边存在dup
        #     reads=[split_read[i] for i in block[4]]
        #     plot(reads)
        #     resolve_dup(reads,dup_region)
        num+=1
    #窗口移动
    #按照在aln顺序 三个reads为一个窗口进行移动 每次移动加2
    #ref窗口自适应
    # True_list=[]
    # double_list=[]
    # backbone_list=[]
    # block_scan=[]
    # for i in split_read:
    #     print(str(i[al_ref_num])+','+str(i[al_is_reverse]),end='\t')
    # print('\n')
    
    inversions=[]
    # plot(split_read)
    if count>=2:
        a=0
        while a < count-1:
            flag=False
            ele1=split_read[a][:]
            ele2=split_read[a+1][:]
            ref_num1=ele1[al_ref_num]
            ref_num2=ele2[al_ref_num]
            if (ref_num1+1 == ref_num2) or ref_num1== ref_num2:
                flag=True
                
            if ele2[al_chr_name]==ele1[al_chr_name]:#chr
                if ele1[al_is_reverse]==ele2[al_is_reverse]:#strand
                    
                    if ele2[al_is_reverse]:#反转 修改坐标

                        query1=query
                        if ele1[al_r_start] > ele2[al_r_start]:
                        
                            q_start=query_length-ele1[al_q_end]
                            q_end=query_length-ele1[al_q_start]
                            ele1[1],ele1[2]=q_start,q_end

                            q_start=query_length-ele2[2]
                            q_end=query_length-ele2[1]
                            ele2[1],ele2[2]=q_start,q_end
                            e=ele1[:]
                            ele1=ele2[:]
                            ele2=e[:]
                            query1=query[::-1]
                    else:
                        query1=query
                        #split_read read    q_start q_end   mapq    r_start r_end   is_reversed chr_name    aln_num ref_num
                        #             0         1     2       3         4     5         6           7           8     9
                                            
                    dis_ref = ele2[4] - ele1[5]
                    dis_read = ele2[1] - ele1[2]
                    delta_length = dis_read-dis_ref
                    if flag:#表明是相邻s

                    
                        if dis_ref >= - max(SV_size, delta_length/5) and delta_length >= SV_size: #no overlap in ref
                            if dis_ref <= max(100, delta_length/5) and (delta_length <= MaxSize or MaxSize == -1):#no gap in ref
                                candidates[types_to_output[ins_num]].append([ele2[al_chr_name],(ele2[al_r_start]+ele1[al_r_end])/2,delta_length,readname,str(query1[ele1[2]+int(dis_ref/2):ele2[1]-int((dis_ref)/2)])])
                        elif dis_ref < -SV_size and delta_length>=SV_size: #overlap on ref
                            if dis_read >= - max(SV_size, delta_length/2) and  dis_read <=max(100, delta_length/2):#no overlap or gap in read
                                
                                if delta_length <=MaxSize or MaxSize == -1:#正常dup   没有分reverse处理
                                    candidates[types_to_output[duptan_num]].append([ele2[al_chr_name],ele2[al_r_start],ele1[al_r_end],readname,bam])
                            elif ele1[al_r_start]==ele2[al_r_start] or ele1[al_r_end]==ele2[al_r_end]:#完全一样的
                                if ele2[al_r_end]<=ele1[al_r_end]:
                                    if ele1[al_q_end]-ele2[al_q_end]>=0:
                                        candidates[types_to_output[duptan_num]].append([ele2[al_chr_name],ele2[al_r_end]-(ele1[al_q_end]-ele2[al_q_end]),ele2[al_r_end],readname])
                                    else:
                                        candidates[types_to_output[duptan_num]].append([ele2[al_chr_name],ele2[al_r_start],ele2[al_r_end],readname])
                                else:
                                    candidates[types_to_output[duptan_num]].append([ele2[al_chr_name],ele2[al_r_start],ele1[al_r_end],readname])
                            elif ele2[al_r_end]<ele1[al_r_start]:#trans
                                if dis_read >- max(SV_size, delta_length/5) :#无overlap
                                    candidates[types_to_output[bnd_num]].append([ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_start],"Same trans",readname])
                            # elif dis_read <-SV_size:#read重叠的一块del
                            #     if delta_length>0:
                            #         print("WHy")
                            #     candidates.append(CandidateDeletion(ele2[7],ele1[al_r_end],-delta_length,[readname]))
                            else:
                                print("hh")

                        if -dis_ref < max(SV_size, -delta_length/2) and -delta_length >= SV_size:#no overlap in ref
                                    if dis_read <= max(100, -delta_length/2) and (-delta_length <= MaxSize or MaxSize == -1):#no gap in read
                                        if ele1[al_ref_num]!=-1 or (ele1[al_ref_num]==-1 and -delta_length <10000) :
                                            candidates[types_to_output[del_num]].append([ele2[7],ele1[5],-delta_length,readname])
                                        else:
                                            candidates[types_to_output[bnd_num]].append([ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_start],"Same trans",readname])
                    else:#same trans
                        if dis_read >= - max(SV_size, delta_length/5) and  dis_read <=max(100, delta_length/5):#no overlap or gap in read
                            if dis_ref>max(100, delta_length/5) :#has gap in ref
                                candidates[types_to_output[bnd_num]].extend(analysis_trans(ele1,ele2,readname,bam))
                        else:
                            #read之间有点远 所以可能有ins+trans   #因为重复区域乱比对导致的
                            try:
                                if not ele1[al_is_reverse]:
                                    new_dis_ref=new_ref_blocks[ele1[al_ref_num]+1][0]-new_ref_blocks[ele1[al_ref_num]][1]
                                    new_ele1=[ele1[0],ele1[1],new_ref_blocks[ele1[al_ref_num]+1][0]]+ele1[3:]
                                else:
                                    new_dis_ref=new_ref_blocks[ele1[al_ref_num]-1][1]-new_ref_blocks[ele1[al_ref_num]][0]
                                    new_ele1=[ele1[0],new_ref_blocks[ele1[al_ref_num]-1][1],ele1[2]]+ele1[3:]
                                new_detalength=dis_read-new_dis_ref
                                if new_dis_ref >= - max(SV_size, new_detalength/5) and new_dis_ref<= max(100, new_detalength/5):#no overlap or gap in ref
                                    if new_detalength >= SV_size and (abs(new_detalength) <= MaxSize or MaxSize == -1) :
                                        candidates[types_to_output[ins_num]].append([ele2[al_chr_name],(ele2[al_r_start]+ele1[al_r_end])/2,new_detalength,[readname],str(query1[ele1[2]+int(new_dis_ref/2):ele2[1]-int((new_dis_ref)/2)])])
                                if -new_dis_ref < max(SV_size, -new_detalength/5) and -new_detalength >= SV_size:#no overlap in ref
                                    if dis_read <= max(100, -new_detalength/5) and (abs(new_detalength) <= MaxSize or MaxSize == -1):#no gap in read
                                        candidates[types_to_output[del_num]].append([ele2[7],ele1[5],-delta_length,readname])
                                   
                                candidates[types_to_output[bnd_num]].extend(analysis_trans(new_ele1,ele2,readname,bam))

                            except:
                                candidates[types_to_output[bnd_num]].extend(analysis_trans(ele1,ele2,readname,bam))
                                # print("else")


                else:#dif strand
                    if ele1[6]:#-+
                        
                        if a+2<count:
                            ele3=split_read[a+2][:]
                            if  ele3[al_is_reverse]:
                                 
                                if ele2[al_q_start] + 0.5 * (ele1[al_r_start] - ele3[al_r_end]) >= ele3[al_q_end] and ele1[al_q_start] + 0.5 * (ele1[al_r_start] - ele3[al_r_end]) >= ele2[al_q_end]:#正常在中间的iver
                                    # No overlaps in split reads

                                    if ele2[al_r_start]-ele3[al_r_end]>= -SV_size and ele1[al_r_start] - ele2[al_r_end]>=-SV_size:#这个和下边的不太一样
                                        candidates[types_to_output[inv_num]].append([ele1[al_chr_name],ele3[al_r_end],ele2[al_r_end],"++",readname])
                                        # head-to-head
                                        # 5'->5'
                                        candidates[types_to_output[inv_num]].append([ele1[al_chr_name],ele2[al_r_start],ele1[al_r_start],"--",readname])
                                        
                                        # tail-to-tail
                                    else:
                                        if not flag:#有same trans
                                            candidates[types_to_output[bnd_num]].append([ele1[al_chr_name],ele3[al_r_end],ele2[al_chr_name],ele2[al_r_start],"-+- same trans",readname])
                                        else:
                                            print("has 遗漏")
                                    a+=1
                                        
                            else:
                                ele_1=[ele1[1],ele1[2],ele1[4],ele1[5],ele1[7],ele1[6],ele1[0]]
                                ele_2=[ele2[1],ele2[2],ele2[4],ele2[5],ele2[7],ele2[6],ele2[0]]
                                analysis_inv(ele_1, 
                                            ele_2, 
                                            readname, 
                                            candidates, 
                                            SV_size,bam)
                                inversions.append([ele1[al_chr_name],ele2[al_r_start],ele2[al_r_end],"-++"])
                    #Normal to reverse
                    #这些处理什么意思
                    elif ele2[6]:#+-
                        if a+2<count:
                            ele3=split_read[a+2][:]
                            if not ele3[al_is_reverse]:
                                if ele2[al_q_start] + 0.5 * (ele3[al_r_start] - ele1[al_r_end]) >= ele1[al_q_end] and ele3[al_q_start] + 0.5 * (ele3[al_r_start] - ele1[al_r_end]) >= ele2[al_q_end]:#两边没有太大的gap
                                    # No overlaps in split reads

                                    if ele2[al_r_start]-ele1[al_r_end]>= -SV_size and ele3[al_r_start] - ele2[al_r_end]>=-SV_size:#保证比对也在中间
                                        candidates[types_to_output[inv_num]].append([ele1[al_chr_name],ele1[al_r_end],ele2[al_r_end],"++",readname])
                                        # head-to-head
                                        # 5'->5'
                                        candidates[types_to_output[inv_num]].append([ele1[al_chr_name],ele2[al_r_start],ele3[al_r_start],"--",readname])
                                        inversions.append([ele1[al_chr_name],ele1[al_r_end],ele2[al_r_end],"++"])
                                    else:
                                        if not flag:#有same trans
                                            candidates[types_to_output[bnd_num]].append([ele1[al_chr_name],ele1[al_r_end],ele2[al_chr_name],ele2[al_r_start],"+-+ same trans",readname])
                                        else:
                                            print("has 遗漏")
                                    a+=1
                                       
                                        # tail-to-tail
                                        # 3'->3'
                            #print("Overlapping read segments in read", read_name)
                            else:
                                ele_1=[ele1[1],ele1[2],ele1[4],ele1[5],ele1[7],ele1[6],ele1[0]]
                                ele_2=[ele2[1],ele2[2],ele2[4],ele2[5],ele2[7],ele2[6],ele2[0]]
                                analysis_inv(ele_1, 
                                            ele_2, 
                                            readname, 
                                            candidates, 
                                            SV_size,bam)
                                inversions.append([ele1[al_chr_name],ele2[al_r_start],ele2[al_r_end],"+--"])
                   
            else:#dif chr trans
                if ele2[al_ref_num]==-1:
                    if ele1[al_ref_num]!=-1:
                        
                        candidates[types_to_output[bnd_num]].extend(analysis_trans(ele1,ele2,readname,bam))
                                                    
            a+=1
    else:
        print("1")
    # if len(inversions)>0:
    #     print(len(inversions))
    return candidates    
# 对于第一次的版本
     # last_ref=current_aln[0]
                
    # for index in range(1,len(alignment_list)):
        
    #     dirty=False
    #     if alignment_list[index][0] not in refpos_dic:
    #         refpos_dic[alignment_list[index][0]]=[]#存储起始，终止，index,chr
    #         refpos_dic[alignment_list[index][0]].append((alignment_list[index][3],alignment_list[index][4],index,alignment_list[index][0]))
    #     else:
    #         max_refpos=refpos_dic[alignment_list[index][0]][-1][1]
    #         lastindex=refpos_dic[alignment_list[index][0]][-1][2]
    #         if  alignment_list[index][0]==last_ref:#一样染色体
    #             if alignment_list[index][al_r_end]>max_refpos:#并且ref正常增大
                    
    #                 if alignment_list[index][al_q_end]>max_alnpos+thred:# aln 正常增大
                        
    #                     # if index+1<len(alignment_list): 想要利用下一个来判断 但好像不行
    #                         # if alignment_list[index][al_q_end]<alignment_list[index+1][al_q_start]:
    #                         #     print('error')
    #                         #     continue
                        
    #                     #处理相邻染色体和aln之间距离
    #                     # print(alignment_list[index][al_q_start]-alignment_list[lastindex][al_q_end])
                       
                       
                        
    #                     if alignment_list[index][al_r_start]-max_refpos-(alignment_list[index][al_q_start]-max_alnpos)>max_thred :#ref异常增大 或者大于下一个的结束位置!!!!参数需要调整
    #                         #aln+ ref++ 易位
    #                         ele=alignment_list[lastindex]
    #                         print('pasetindex',index,lastindex)
    #                         detlength=alignment_list[index][al_q_start]-ele[al_q_end]
    #                         paste.append([0,ele[al_chr_name],ele[al_is_reverse],ele[al_r_start],ele[al_r_end],alignment_list[index][al_chr_name],flag,alignment_list[index][al_r_start], alignment_list[index][al_r_end],detlength,readname])
    #                         dirty=True
    #                     max_alnpos=alignment_list[index][al_q_end]
    #                     if not dirty:
    #                         nearby_collect(alignment_list[lastindex][:],alignment_list[index][:],query,query_length,readname,ins_candidates,del_candidates,duptans,dupinters)
    #                         refpos_dic[last_ref].append((alignment_list[index][3],alignment_list[index][4],index,alignment_list[index][0]))   
    #                         # max_refpos=alignment_list[index][al_r_end]
                
    #                 else:#aln 变小 
    #                     #case 1 内部反转信号 inv
    #                     #case 2 bnd 可能是3种tra 比对到其他地方？
    #                     if alignment_list[index][al_r_start]-max_refpos>max_thred :#
    #                         ele=alignment_list[lastindex]
    #                         bnd.append([1,lastindex,ele[al_is_reverse],ele[al_r_start],ele[al_r_end],index,flag,alignment_list[index][al_r_start], alignment_list[index][al_r_end]])
    #                         dirty=True
    #                     else:
    #                         print("?",index)
    #                     pass
                
    #             else:#ref 变小
    #                 if alignment_list[index][al_q_end]>max_alnpos+thred:#aln 正常变大
    #                     #说明是多了一部分 
    #                     # case 1:也可能是cutpaste（另一个部分有正好少一样的部分）和
    #                     # case 2:copy dup（另一部分没有少 可以多） 通过另一个位置的远近来判断 tandup or interdup
    #                     # case 3:#先不考虑同一个染色体上的易位 他会突然变很大和很小（另一部分也会有相似的段位信号，但长度不一定）
                       
    #                         #如果都放到后边判断是哪种类型
    #                         #易位
                        
    #                     ele=alignment_list[lastindex]
    #                     print('pasetindex',index,lastindex)
    #                     detlength=alignment_list[index][al_q_start]-ele[al_q_end]
    #                     paste.append([1,alignment_list[index][al_chr_name],flag,alignment_list[index][al_r_start], alignment_list[index][al_r_end],ele[al_chr_name],ele[al_is_reverse],ele[al_r_start],ele[al_r_end],detlength,readname])
    #                     max_alnpos=alignment_list[index][al_q_end]
    #                 else: #aln 变小
    #                     print("- -???",index)
    #                     pass
    #         else:#不一样染色体
    #             #bnd
    #             if edalignment_list[index][al_q_end]>max_alnpos+thr:
    #                 max_alnpos=alignment_list[index][al_q_end]
    #             pass


#
# 输入: svcandidates, refcandidates
# 输出: svcandidates_no_fp

# 初始化去除假阳性后的候选集
svcandidates_no_fp = {svtype: [] for svtype in types_to_output}

ref_buckets = dict() # 初始化参考桶
base = 100 #初始化base阈值

# 构建参考桶
FOR svtype, refsv_list in refcandidates.items() #对于每一种变异类型和对应的变异信号列表
    FOR refsv in refsv_list
        chrom,start = refsv[chrom],refsv[start] #获取变异信号的染色体和起始位置
        bucket = int(start / base)#确定桶号

        IF chrom not in ref_buckets THEN
            ref_buckets[chrom] = dict()
        ENDIF
        IF svtype not in ref_buckets[chrom] THEN
            ref_buckets[chrom][svtype] = dict()
        ENDIF
        IF bucket not in ref_buckets[chrom][svtype] THEN
            ref_buckets[chrom][svtype][bucket] = list()
        ENDIF

        ref_buckets[chrom][svtype][bucket].append(refsv)
ENDFOR
# 筛选待检查候选数据
FOR svtype, svcand_list in svcandidates.items():
    FOR svcand in svcand_list
        chrom,start = svcand[chrom],svcand[start]
        bucket = int(start / base)
        flag = 0#是否匹配标志

        if chrom in ref_buckets and svtype in ref_buckets[chrom]:
            relevant_buckets = [b for b in ref_buckets[chrom][svtype].keys() if abs(b - bucket) < 1]#找到关联桶
            for bucket_start in relevant_buckets:
                for refsv in ref_buckets[chrom][svtype][bucket_start]:#对于可能的每个变异信号
                    IF Same(refsv,svcand) THEN
                        flag=1
                    ENDIF
                IF flag == 1:#找到则退出循环
                    break
                ENDIF

        IF flag == 0 THEN
            svcandidates_no_fp[svtype].append(svcand)#添加到过滤后集合
        ENDIF
    ENDFOR
ENDFOR

RETURN svcandidates_no_fp


# 自匹配
# 初始化两个数据集合 A 和 B，对应不同类型的bnd，以及结果集 result
A = bnd[0]
B = bnd[1]
result = list()

foreach elementA in A:
    chromosomeA,startA,lengthA = elementA[chromosome],elementA[start],elementA[length]
    foreach elementB in B:
        chromosomeB,startB,lengthB = elementB[chromosome],elementB[start],elementB[length]
        if chromosomeA == chromosomeB:
            if abs(lengthA - lengthB) <= threshold:
                if elementB[status]!= marked:
                    new_entry = create_new_entry(elementA, elementB)
                    result.append(new_entry)
                    elementB[status] = marked
# invs添加方向

# 假设 bnd_can 和 invs 是已有的数据集合
foreach element1 in bnd_can do #对每一个bnd
    chrom1_1 = element1[chromosome1]
    start1_1 = element1[start1]
    chrom1_2 = element1[chromosome2]
    start1_2 = element1[start2]
    length1 = element1[length]

    foreach element2 in invs do #遍历反转信号
        IF element2[status]!= marked THEN
            chrom2,start2,length2 = element2[chromosome],element2[start],element2[length]

            IF chrom2 == chrom1_1 THEN
                IF abs(length1 - length2) <= length_threshold and abs(start1_1 - start2) <= start_threshold THEN
                    element1[flag1] = True #添加方向信息
                    element2[status] = marked
                ENDIF
            ELIF chrom2 == chrom1_2 THEN
                IF abs(length1 - length2) <= length_threshold and abs(start2 - start1_2) <= start_threshold THEN
                    element1[flag2] = True #添加方向信息
                    element2[status] = marked                    
                ENDIF
            ENDIF
        ENDIF

# paste _ del

# 假设已有以下抽象函数定义：
# del_bucket: 对删除列表进行分桶处理，返回分桶后的结构
# parse_candidate: 解析候选数据，返回关键属性
# find_del_in_bucket: 在分桶数据中查找相关删除数据
# simi_k: 计算两个序列相似度，返回相似度值和布尔值
# Merge:合并匹配的del和ins/paste，返回融合成功的bnd信号
# Counter: 统计元素出现次数的工具
# 输入：粘贴数据paste，删除列表del_list，参考基因组reference，重复相关数据dup_tan、dup_int，结果列表rm_cut
# 输出：更新后的删除列表

function paste_del_match(paste, del_list, reference, dup_tan, dup_int):
    cut_paste = []
    IF is_empty(paste) THEN
        return cut_paste,dup_int，dup_tan
    ENDIF
    del_budget = del_bucket(del_list)
    FOREACH candidate IN paste DO
        flag = False
        readname = candidate[readname]
        detalen = candidate[detalen]
        chrom, start, end, last_start, last_end, ins_flag, last_flag = parse_candidate(candidate)

        del_result = find_del_in_bucket(del_budget, chrom, start, end, end - start)
        IF length(del_result) < MAX_PAIR THEN #当候选信号少的时候执行
            query_seq = reference.fetch(chrom, start, end).to_upper() #获取该paste序列
            k_and_bool = [(simi_k(query_seq, del_can, reference)) FOREACH del_can IN del_result]#调用simi_k存储所有的相似度结果
            IF not is_empty(k_and_bool) THEN
                k_values, bool_values = unzip(k_and_bool)
                min_k = min(k_values)
                IF min_k < MIN_DISTANCE THEN #当小于距离相似度阈值
                    min_k_index = index_of(min_k, k_values)
                    is_true = bool_values[min_k_index]
                    cut_paste.append(Merge(del_result[min_k_index],candidate))
                    flag = True
                ENDIF
            ENDIF
        ENDIF

        IF not flag THEN
            handle_dup(last_end, dup_tan, dup_int, chrom, start, end, detalen, ins_flag, readname) #没匹配成功的进一步判断是dupint/duptan
        ENDIF
RETURN cut_paste,dup_int，dup_tan

# semi_k

# 假设已有以下抽象函数定义：
# align: 执行序列比对操作，返回包含比对结果的字典
# Seq: 用于处理序列的类，这里使用它生成反向互补序列
# reference.fetch: 从参考基因组获取指定区域序列

# 输入：查询序列queryseq，待比对的目标数据del_can，参考基因组reference
# 输出：编辑距离edi_dis与布尔标志flag组成的元组

function simi_k(queryseq, del_can, reference):
    flag = False
    chrom = del_can[chromosome_index]
    region_start = del_can[start_index]
    region_end = del_can[length_index] + region_start
    targetseq = reference.fetch(chrom, region_start, region_end).to_upper()

    align_ans = align(queryseq, targetseq, mode="NW", task="distance", k=20)
    edi_dis = align_ans[edit_distance_key]
    if edi_dis == -1:
        query_rc = reverse_complement(queryseq)
        align_ans1 = align(query_rc, targetseq, mode="SHW", task="distance", k=20)
        if align_ans1[edit_distance_key] == -1:
            return (-1, flag)
        edi_dis = align_ans1[edit_distance_key]
        flag = True
    if edi_dis == -1:
        edi_dis = 100

    return (edi_dis, flag)