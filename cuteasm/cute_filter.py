#将alignment reshape 使头到尾
#操纵cigar字符串 和参考序列
#输入bam 输出bam 修正好的bam
import pysam
import logging
import matplotlib.pyplot as plt
import numpy as np

###找到在某条染色体上的所有同名字的read集合，good_dic是指含有主比对的，然后按照query排序
# #获得特定read名的堆
# read_dic={}
# while True:
#     try:
#         current_alignment=next(alignment_it)
#         if current_alignment.query_name not in read_dic:
#             read_dic[current_alignment.query_name]=[current_alignment]
#         else:
#             read_dic[current_alignment.query_name].append(current_alignment)
#     except StopIteration:
#         break
# #筛选主比对
# good_readdic={}
# for read in read_dic:
#     flag=0
#     for alignment in read_dic[read]:
#         if  not alignment.is_supplementary and not alignment.is_secondary:
#             good_readdic[read]=read_dic[read]
#             break


# #并按照在query相对顺序排序
# for read in good_readdic:
#     good_readdic[read]=sorted(good_readdic[read],key=lambda x:x.query_alignment_start)
#     #对每一个read重构
#     for alignment in good_readdic[read]:
#         print(alignment.cigarstring)
#         print(alignment.is_unmapped)
#         print(alignment.is_secondary)
#         print(alignment.mapping_quality)
        

# #写到文件


###
#找到所有的主比对的SA
#将SA可视化
##首位连接 不同质量比对不同颜色  竖条纹   如果可以，在ref的相对位置


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

        other_alignments.append([a,0])
    return other_alignments
def get_length(cigartuple):#返回query序列涉及的长度
    qlen=0
    rlen=0
    for operation, length in cigartuple:
        if operation == 0 or operation ==7 or operation == 8:
            qlen+=length
            rlen+=length
        elif operation==2:
            rlen+=length
        elif operation == 1 :
            qlen+=length
    if cigartuple[0][0]!=4:
        start=0
    return qlen,rlen

# 创建一个颜色映射
mapq_colors = {
    10: 'red',
    30: 'orange',
    40: 'yellow',
    50: 'green',
    60: 'blue'
}

def draw_alns(reads_info,num,mode=1):
    # 绘图
    plt.figure(figsize=(12, 6))

    # 遍历读段信息并绘制
    for read in reads_info:
        qstart, qend, chr_name, rstart, rend, mapQ = read
        colors = plt.cm.viridis(mapQ / 60)
        plt.plot([rstart, rend], [qstart, qend], color=colors, linewidth=2, label=f'MAPQ: {mapQ}')
        # plt.plot([0, rstart], [qstart, qstart], color='red', linestyle='--', linewidth=1)
        # plt.plot([0, rend], [qend, qend], color='green', linestyle='--', linewidth=1)

    # 设置图例
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), title="MAPQ")

    # 设置坐标轴标签和标题
    plt.xlabel("Reference Position")
    plt.ylabel("Query Position")
    plt.title("Visualization of Read Positions on Chromosome")
    if mode==2:
        mane="new_read_positions"+str(num)+".png"
    else:
        mane="read_positions"+str(num)+".png"
    plt.savefig(mane, dpi=300, bbox_inches='tight') 

    # 显示图形
    plt.grid()
    plt.show()     
#找到对应的read
def fetch_read(sa_name,sa_chr,sa_pos,start_pos,bam_file):
    sa_length=200#默认200bp
    for sa_read in bam_file.fetch( sa_pos, sa_pos + sa_length,tid=sa_chr):
        if sa_read.reference_start == sa_pos and sa_read.query_name== sa_name and query_start(sa_read)==start_pos:
                #print(f"    Found Secondary Alignment: {sa_read.query_name}, Start: {sa_read.reference_start}")
            return sa_read
def fetch_read_dic(sa_name,sa_chr,sa_pos,bam_file):
    sa_length=200#默认200bp
    sa_dict=[]
    for sa_read in bam_file.fetch( sa_pos, sa_pos + sa_length,tid=sa_chr):
        if sa_read.reference_start == sa_pos and sa_read.query_name== sa_name:
                #print(f"    Found Secondary Alignment: {sa_read.query_name}, Start: {sa_read.reference_start}")
                sa_dict.append(sa_read)
    return sa_dict
#重叠部分的编辑距离
def overlap_nm(read1,overlap,string=1):
    nm=0
    if string==2:
        for item in read1.cigartuples[::-1]:
            if item[0]==4 or item[0]==5:
                continue
            elif item[0]==1:
                nm+=item[1]
                overlap-=item[1]
                if overlap<=0:
                    nm=nm+overlap
                    break
            elif item[0]==2:
                nm+=item[1]
            else:
                overlap-=item[1]
                if overlap<=0:
                    break
    else:
        for item in read1.cigartuples:
            if item[0]==4 or item[0]==5:
                continue
            elif item[0]==1:
                nm+=item[1]
                overlap-=item[1]
                if overlap<=0:
                    nm=nm+overlap
                    break
            elif item[0]==2:
                nm+=item[1]
            else:
                overlap-=item[1]
                if overlap<=0:
                    break

    return nm

def modify_cigar_for_overlap(alignment, overlap_length,string=1):
    """
    Modify the CIGAR string of the given alignment to change the last `overlap_length` bases
    to either hard clipping (H) or soft clipping (S).

    Parameters:
    - alignment: The alignment object to modify.
    - overlap_length: The number of overlapping bases to modify.
    string==1 表示修改前一个
    """
    # 获取当前的 CIGAR 字符串
    cigar = alignment.cigar
    mode='S'
    
    new_cigar = []

    # 计算当前的查询长度
    query_length = alignment.query_alignment_length+query_start(alignment)
    cigar_length=0

    # 处理 CIGAR 字符串
    total_length = 0
    flag_1=0
    ref_length=0
    if string==1:
        if alignment.cigartuples[-1][0]==5:
            mode='H'
        for  operation,length in cigar:
            if operation ==2:
                continue
            total_length += length
            if total_length < query_length - overlap_length:
                new_cigar.append((operation,length))  # 保留原有的 CIGAR
            else:
                if flag_1==0:
                    flag_1=1
                    # 修改重叠部分
                    new_cigar.append(( operation,query_length - total_length + length-overlap_length))
                    cigar_length+=-query_length +total_length + overlap_length
                else:
                    cigar_length+=length
                    if operation==5 or operation ==4:
                        new_cigar.append(( operation,cigar_length))


                
            
    else:
        if alignment.cigartuples[0][0]==5:
            mode='H'
        for operation,length in cigar:
            if (operation==4 or operation==5 )and flag_1==0:
                cigar_length+=length
                continue
            if operation ==2:
                ref_length+=length
                continue
            total_length += length
            if total_length <= overlap_length:
                # 修改重叠部分
                if operation==0 or operation ==7 or operation ==8:
                    ref_length+=length
                    cigar_length+=length                        
            else:
                if flag_1==0:
                    flag_1=1
                    if operation==0 or operation ==7 or operation ==8:
                            ref_length+=length-total_length+overlap_length
                    if mode == 'H':
                        new_cigar.append(( 5,length-total_length+overlap_length+cigar_length))  # 5 = H
                        
                    elif mode == 'S':
                        new_cigar.append((4,length-total_length+overlap_length+cigar_length))  # 4 = S
                    new_cigar.append(( operation,total_length-overlap_length))
                else:
                    new_cigar.append(( operation,length))
                
    # 更新 CIGAR
    alignment.cigar = new_cigar
    if string==2:
        alignment.reference_start = alignment.reference_start+ref_length

def query_start(alignment):
    total_length=0
    for  operation,length in alignment.cigartuples:
        if operation in (4, 5):  # 4 = S (soft clip), 5 = H (hard clip)
            total_length += length
        else:
            break  # 一旦遇到非剪切操作，停止计算

    return total_length
def query_start_cigar(cigar_string):
    start_position=0
      # 解析 CIGAR 字符串
    i = 0
    while i < len(cigar_string):
        # 读取数字部分
        num_str = ''
        while i < len(cigar_string) and cigar_string[i].isdigit():
            num_str += cigar_string[i]
            i += 1
        
        if num_str:
            length = int(num_str)  # 转换为整数
            if i < len(cigar_string):
                operation = cigar_string[i]  # 获取操作符
                if operation == 'S' or operation == 'H'  :
                    start_position += length  # 匹配会影响起始位置
                else:
                    break
                i += 1  # 移动到下一个字符
    return start_position

def is_same_alignment(read1,read2):
    if query_start(read1)==query_start(read2):
        return 1
    return 0


if __name__ == "__main__":
    ref_d="/data/1/wuxiaomia/data/GRCh38_chromosomes.fa"
    hap1_d="/data/1/wuxiaomia/data/aln/Assembly/NA24385_shiJie/NA24385.hap1.bam"
    hap2_d="/data/1/wuxiaomia/data/aln/Assembly/NA24385_shiJie/NA24385.hap2.bam"
    reference=pysam.FastaFile(ref_d)
    bam=pysam.AlignmentFile(hap1_d)
    hap2bam=pysam.AlignmentFile(hap2_d)


    chromosomes=bam.references
    current_chromosome=chromosomes[0]
    alignment_it=bam.fetch(contig=current_chromosome)
    print("Processing chromosome {0}...".format(current_chromosome))
    minmapQ=20
    minioverlap=200
    #将SA 按照mapQ填补
    primary_reads=[]
    unmaped=[]
    num=0
    #找到主比对
    while True:
        try:
            current_alignment=next(alignment_it)
            num+=1
            if current_alignment.is_unmapped :
                unmaped.append(current_alignment)
            elif  not current_alignment.is_supplementary and not  current_alignment.is_secondary:
                primary_reads.append(current_alignment)
            else:
                continue

        except StopIteration:
            break
    #Filter&Reshape

    #alns_ [[aln,0/-1/1]]
    #alns [[aln,mapping_quality,query_start,0/-1/1]
    #new_alns[[query_start,query_end,ref_start,ref_end,chr,strand,aln]]
    num_i=0#计数alignment
    repeat_region_primary=[]#记录repeat区域的primary
    for pri_aln in primary_reads:
        if num_i==10:
            pass
        if pri_aln.mapping_quality <10:#说明这很可能是重复区域
            #需要对重复区域聚类
            repeat_region_primary.append(pri_aln)
            continue
        suppl_alns=retrieve_other_alignments(pri_aln,bam)
        alns_=[[pri_aln,1]]
        alns_.extend(suppl_alns)
        #所有的同一个contig的alns按照mapQ和起始位置（cigar起始位置第一个大小）排序
        alns=[]
        for aln in alns_:
            alns.append([aln[0],aln[0].mapping_quality,aln[0].cigartuples[0][1],aln[1]])
        alns.sort(key=lambda x:(-x[1],x[2]))
        print("hello")
        new_alns=[]#记录新的不重叠的alns 即alns进行过滤和reshape后的new alns
        alni=0#遍历索引
        alnj=1
        flag_pri=0
        while(alni<len(alns)-1 and alnj<len(alns)):
            while alns[alni][-1]==-1 :
                alni+=1
            aln1=alns[alni][0]
            while alns[alnj][-1]==-1 :
                alnj+=1
            aln2=alns[alnj][0]

            # print(query_start(aln1))
            # print(query_start(aln2))
            if alns[alni][1] ==1:
                flag_pri=1
            if alns[alnj][1]==1:
                flag_pri=2
        
            read1=fetch_read(aln1.query_name,aln1.reference_id ,aln1.reference_start,query_start(aln1),bam)
            read2=fetch_read(aln2.query_name,aln2.reference_id ,aln2.reference_start,query_start(aln2),bam)
            
            overlap=aln1.query_alignment_length+query_start(read1)-query_start(read2)
            #对前后两个alignment
            if aln1.mapping_quality == aln2.mapping_quality:#第一波最好的质量
                ####这里是处理后一个比对包含在前面比对的地方 具体该怎么处理？？？
                if (read2.query_alignment_length+query_start(read2))< (read1.query_alignment_length+query_start(read1)):
                    #new_alns.append([read1,query_start(read1),read1.query_alignment_length+query_start(read1)])
                    #new_alns.append([read2,query_start(read2),read2.query_alignment_length+query_start(read2)])
                    alns[alnj][-1]=-1
                    alnj+=1
                    continue
                #print("mapq:"+str(aln1.mapping_quality))
                if overlap>=200:#有重叠
                    if flag_pri==1:
                        modify_cigar_for_overlap(read2,overlap,2)
                    elif flag_pri==2:
                        modify_cigar_for_overlap(read1,overlap,1)
                    else:
                        read1_nm=overlap_nm(read1,overlap,2)
                        read2_nm=overlap_nm(read2,overlap,1)
                        if read1_nm-read2_nm>100:#如果差别很大就修改read
                            #修改read1的cigar
                            modify_cigar_for_overlap(read1,overlap,1)
                        elif read2_nm-read1_nm>100:
                            modify_cigar_for_overlap(read2,overlap,2)
                new_alns.append([read1,query_start(read1),read1.query_alignment_length+query_start(read1)])
                alni+=1
                alnj+=1
                if alni==(len(alns)-1):
                    if alns[alni][-1]!=-1:
                        new_alns.append([read2,query_start(read2),read2.query_alignment_length+query_start(read2)])
            else:
                new_alns.append([read1,query_start(read1),read1.query_alignment_length+query_start(read1)])
                alni+=1
                alnj+=1
                
                break
            
        if alni<len(alns)-1:#当有次低比对时
            #对于剩下的插入空闲区域
            while alni<len(alns)-1 :
                aln11=alns[alni][0]
                aln1=fetch_read(aln11.query_name,aln11.reference_id ,aln11.reference_start,query_start(aln11),bam)
                if aln1.mapping_quality < minmapQ:
                    alni+=1
                    break
                start=query_start(aln1)
                end=start+aln1.query_alignment_length
                for newalni in range(len(new_alns)):
                
                    if new_alns[newalni][2]<end:
                        continue
                    overlap_1= end-new_alns[newalni][1]
                    
                    if newalni==0:
                        if overlap_1<=0:
                            new_alns.insert(0,[aln1,start,end])
                        else:
                            modify_cigar_for_overlap(aln1,overlap_1)
                            new_alns.insert(0,[aln1,query_start(aln1),aln1.query_alignment_length+query_start(aln1)])
                        alni+=1
                        break
                    overlap_2= start-new_alns[newalni-1][2]
                    if overlap_2>0:
                        if overlap_1<=0:
                            new_alns.insert(newalni,[aln1,start,end])
                        else:
                            if end-start-overlap_1>200:#多的话就插入
                                modify_cigar_for_overlap(aln1,overlap_1)
                                new_alns.insert(newalni,[aln1,query_start(aln1),aln1.query_alignment_length+query_start(aln1)])
                        alni+=1
                        break
                    else:
                        if overlap_1<=0:
                            modify_cigar_for_overlap(aln1,overlap_2,2)
                            new_alns.insert(newalni,[aln1,query_start(aln1),aln1.query_alignment_length+query_start(aln1)])
                        else:
                            if end-start+overlap_2-overlap_1>200:
                                modify_cigar_for_overlap(aln1,overlap_1,1)
                                modify_cigar_for_overlap(aln1,-overlap_2,2)
                                new_alns.insert(newalni,[aln1,query_start(aln1),aln1.query_alignment_length+query_start(aln1)])
                        alni+=1
                        break
                if newalni==len(new_alns)-1:
                    new_alns.insert(newalni,[aln1,query_start(aln1),aln1.query_alignment_length+query_start(aln1)])

    





        




        #可视化相对位置 【qstart,qend,chr,rstart,rend,mapQ]
        pos_info=[]
        flag=0
        pos_info2=[]
        for a_ in alns:
            a=a_[0]
            start=query_start(a)
            q_len,r_len=get_length(a.cigartuples)
            pos_info.append([start,start+q_len,a.reference_id,a.reference_start,a.reference_start+r_len,a.mapping_quality])
        for ax in new_alns:
            a=ax[0]
            start=query_start(a)
            q_len,r_len=get_length(a.cigartuples)
            pos_info2.append([start,start+q_len,a.reference_id,a.reference_start,a.reference_start+r_len,a.mapping_quality])
        pos_info.sort(key=lambda x:x[0])
        print("\t")
        print(str(len(pos_info))+":")
        print(pos_info)
        print(str(len(pos_info2))+":")
        print(pos_info2)
        for i in range(0,len(pos_info2)-1):
            #有重叠大于200bp的绘图
            if pos_info2[i][1]-pos_info2[i+1][0]>200:
                print(i)
                flag=1
                break
        if flag==1:
            print("***%d"%(num_i))
            if num_i==5:
                pass

            draw_alns(pos_info,num_i,1)
            draw_alns(pos_info2,num_i,2)
        num_i+=1

    


def filter_alignments(alns_,bam,minmapQ=20):
    #所有的同一个contig的alns按照mapQ和起始位置（cigar起始位置第一个大小）排序
    alns=[]
    for aln in alns_:
        alns.append([aln,aln.mapping_quality,aln.cigartuples[0][1],0])
    alns.sort(key=lambda x:(-x[1],x[2]))
    new_alns=[]#记录新的不重叠的alns 即alns进行过滤和reshape后的new alns
    alni=0#遍历索引
    alnj=1
    while(alni<len(alns)-1 and alnj<len(alns)):
        while alns[alni][-1]==-1 :
            alni+=1
        aln1=alns[alni][0]
        while alns[alnj][-1]==-1 :
            alnj+=1
        aln2=alns[alnj][0]        
        read1=fetch_read(aln1.query_name,aln1.reference_id ,aln1.reference_start,query_start(aln1),bam)
        read2=fetch_read(aln2.query_name,aln2.reference_id ,aln2.reference_start,query_start(aln2),bam)
        overlap=aln1.query_alignment_length+query_start(read1)-query_start(read2)
        #对前后两个alignment
        if aln1.mapping_quality == aln2.mapping_quality:#第一波最好的质量
            ####这里是处理后一个比对包含在前面比对的地方 具体该怎么处理？？？
            if (read2.query_alignment_length+query_start(read2))< (read1.query_alignment_length+query_start(read1)):
                #new_alns.append([read1,query_start(read1),read1.query_alignment_length+query_start(read1)])
                #new_alns.append([read2,query_start(read2),read2.query_alignment_length+query_start(read2)])
                alns[alnj][-1]=-1
                alnj+=1
                continue
            #print("mapq:"+str(aln1.mapping_quality))
            if overlap>=200:#有重叠
                
                read1_nm=overlap_nm(read1,overlap,2)
                read2_nm=overlap_nm(read2,overlap,1)
                if read1_nm-read2_nm>100:#如果差别很大就修改read
                    #修改read1的cigar
                    modify_cigar_for_overlap(read1,overlap,1)
                elif read2_nm-read1_nm>100:
                    modify_cigar_for_overlap(read2,overlap,2)
            new_alns.append([query_start(read1),read1.query_alignment_length+query_start(read1),read1.reference_start,read1.reference_end,read1.reference_name,read1.is_reverse,read1])
            alni+=1
            alnj+=1
            if alni==(len(alns)-1):
                if alns[alni][-1]!=-1:
                    new_alns.append([query_start(read2),read2.query_alignment_length+query_start(read2),read2.reference_start,read2.reference_end,read2.reference_name,read2.is_reverse,read2])
        else:
            new_alns.append([query_start(read1),read1.query_alignment_length+query_start(read1),read1.reference_start,read1.reference_end,read1.reference_name,read1.is_reverse,read1])
            alni+=1
            alnj+=1
            
            break
        
    if alni<len(alns)-1:#当有次低比对时
        #对于剩下的插入空闲区域
        while alni<len(alns)-1 :
            aln11=alns[alni][0]
            aln1=fetch_read(aln11.query_name,aln11.reference_id ,aln11.reference_start,query_start(aln11),bam)
            if aln1.mapping_quality < minmapQ:
                alni+=1
                break
            start=query_start(aln1)
            end=start+aln1.query_alignment_length
            for newalni in range(len(new_alns)):
            
                if new_alns[newalni][2]<end:
                    continue
                overlap_1= end-new_alns[newalni][1]
                
                if newalni==0:
                    if overlap_1<=0:
                        new_alns.insert(0,[start,end,aln1.reference_start,aln1.reference_end,aln1.reference_name,aln1.is_reverse,aln1])
                    else:
                        modify_cigar_for_overlap(aln1,overlap_1)
                        new_alns.insert(0,[query_start(aln1),aln1.query_alignment_length+query_start(aln1),aln1.reference_start,aln1.reference_end,aln1.reference_name,aln1.is_reverse,aln1])
                    alni+=1
                    break
                overlap_2= start-new_alns[newalni-1][2]
                if overlap_2>0:
                    if overlap_1<=0:
                        new_alns.insert(newalni,[start,end,aln1.reference_start,aln1.reference_end,aln1.reference_name,aln1.is_reverse,aln1])
                    else:
                        if end-start-overlap_1>200:#多的话就插入
                            modify_cigar_for_overlap(aln1,overlap_1)
                            new_alns.insert(newalni,[query_start(read1),read1.query_alignment_length+query_start(read1),read1.reference_start,read1.reference_end,read1.reference_name,read1.is_reverse,read1])
                    alni+=1
                    break
                else:
                    if overlap_1<=0:
                        modify_cigar_for_overlap(aln1,overlap_2,2)
                        new_alns.insert(newalni,[query_start(read1),read1.query_alignment_length+query_start(read1),read1.reference_start,read1.reference_end,read1.reference_name,read1.is_reverse,read1])
                    else:
                        if end-start+overlap_2-overlap_1>200:
                            modify_cigar_for_overlap(aln1,overlap_1,1)
                            modify_cigar_for_overlap(aln1,-overlap_2,2)
                            new_alns.insert(newalni,[query_start(read1),read1.query_alignment_length+query_start(read1),read1.reference_start,read1.reference_end,read1.reference_name,read1.is_reverse,read1])
                    alni+=1
                    break
            if newalni==len(new_alns)-1:
                 new_alns.insert(newalni,[query_start(read1),read1.query_alignment_length+query_start(read1),read1.reference_start,read1.reference_end,read1.reference_name,read1.is_reverse,read1])
#[query_start(read1),read1.query_alignment_length+query_start(read1),read1.reference_start,read1.reference_end,read1.reference_name,read1.is_reverse,read1]
#[aln1,query_start(aln1),query_start(aln1)+aln1.query_alignment_length,mapq,aln1.reference_start,0]
    #read_start	read_end	ref_start	ref_end	chr	strand aln
    #0			#1			#2			#3		#4	#5     #6
    return new_alns

#为构建SA中的完整序列而做
#输入：SA标签
#输出：完整的contig经过筛选后的reads

def retrieve_other_alignments1(main_alignment, bam,header):
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
        query_start=query_start_cigar(cigar)
        query_length,rlength=get_length(a.cigartuples)
        try:
            refename=a.reference_name
        except:
            refename= header['SQ'][a.reference_id]['SN']

        other_alignments.append([a,query_start,query_start+query_length,mapq,pos,a.reference_end,a.is_reverse,refename,0])
    
    return other_alignments
def filter_alignments_in_contig(alns,bam):
    #形成前三桶 构成 contig 并且将其中的
    alns_.sort(key=lambda x:(-x[3],x[1]))
    new_alns=[]
    buget=0
    current_mapq=90
    for aln in alns:
        mapq=aln[3]
        a=aln[0]
        if current_mapq != mapq:
            current_mapq=mapq
            buget+=1
        if buget>=4:
            break
        if buget==1 :
            new_alns.append(aln)
        else:
            start=aln[1]
            end=aln[2]
            if mapq>=59:#避免重复提取信号
                aln1=aln[0]
            else:
                aln1=fetch_read(a.query_name,a.reference_id ,a.reference_start,query_start(a),bam)
            for newalni in range(len(new_alns)):
            
                if new_alns[newalni][2]<end:
                    continue
                overlap_1= end-new_alns[newalni][1]
                
                if newalni==0:
                    if overlap_1<=0:
                        aln[0]=aln1
                        new_alns.insert(0,aln)
                    else:
                        modify_cigar_for_overlap(aln1,overlap_1)
                        if aln1.query_alignment_length>200:
                            new_alns.insert(0,[aln1,query_start(aln1),query_start(aln1)+aln1.query_alignment_length,mapq,aln1.reference_start,aln1.reference_end,aln1.is_reverse,0])
                    break
                overlap_2= start-new_alns[newalni-1][2]
                if overlap_2>0:
                    if overlap_1<=0:
                        aln[0]=aln1
                        new_alns.insert(newalni,aln)
                    else:
                        if end-start-overlap_1>200:#多的话就插入
                            aln1=fetch_read(a.query_name,a.reference_id ,a.reference_start,query_start(a),bam)
                            modify_cigar_for_overlap(aln1,overlap_1)
                            if aln1.query_alignment_length>200:
                                new_alns.insert(newalni,[aln1,query_start(aln1),query_start(aln1)+aln1.query_alignment_length,mapq,aln1.reference_start,aln1.reference_end,aln1.is_reverse,0])
                    break
                else:
                    
                    if overlap_1<=0:
                        modify_cigar_for_overlap(aln1,overlap_2,2)
                        if aln1.query_alignment_length>200:
                            new_alns.insert(newalni,[aln1,query_start(aln1),query_start(aln1)+aln1.query_alignment_length,mapq,aln1.reference_start,aln1.reference_end,aln1.is_reverse,0])
                    else:
                        if end-start+overlap_2-overlap_1>200:
                            modify_cigar_for_overlap(aln1,overlap_1,1)
                            modify_cigar_for_overlap(aln1,-overlap_2,2)
                            if aln1.query_alignment_length>200:
                                new_alns.insert(newalni,[aln1,query_start(aln1),query_start(aln1)+aln1.query_alignment_length,mapq,aln1.reference_start,aln1.reference_end,aln1.is_reverse,0])
                    break
            if newalni==len(new_alns):
                 aln[0]=aln1
                 new_alns.append(aln)
    return new_alns#按start顺序排好的完整contig 并且返回了mapq《59的aln
min_length=500
def filter(alns,bam):
    #形成前三桶 构成 contig 并且将其中的
    alns.sort(key=lambda x:(-x[3],x[1]))
    new_alns=[]
    buget=0
    current_mapq=90
    for aln in alns:
        if aln[0].is_unmapped:
            continue
        len_=aln[2]-aln[1]
        if len_<min_length:
            continue
        mapq=aln[3]
        if current_mapq != mapq:
            current_mapq=mapq
            if mapq<20:
                break
            buget+=1
        if buget>=4:
            break
        if buget==1 :
            try:
                bef=new_alns[-1]
                if bef[2]>=aln[2]:#过滤掉重合部分
                    continue
            except:
                pass
            new_alns.append(aln)
        else:
            start=aln[1]
            end=aln[2]
            # if mapq>=59:#避免重复提取信号
            #     aln1=aln[0]
            # else:
            #     aln1=fetch_read(a.query_name,a.reference_id ,a.reference_start,query_start(a),bam)
            for newalni in range(len(new_alns)):
            
                if new_alns[newalni][2]<end:
                    continue
                if start>new_alns[newalni][1]:
                    continue
                overlap_1= end-new_alns[newalni][1]
                
                if newalni==0:
                    if overlap_1<=0:
                        new_alns.insert(0,aln)
                    else:
                        # modify_cigar_for_overlap(aln1,overlap_1)
                        if end-start>500:
                            new_alns.insert(0,aln)
                    break
                overlap_2= start-new_alns[newalni-1][2]
                if overlap_2>0:
                    if overlap_1<=0:
                        new_alns.insert(newalni,aln)
                    else:
                        if end-start-overlap_1>200:#多的话就插入
                            # aln1=fetch_read(a.query_name,a.reference_id ,a.reference_start,query_start(a),bam)
                            # modify_cigar_for_overlap(aln1,overlap_1)
                            if end-start>500:
                                new_alns.insert(newalni,aln)
                    break
                else:
                    
                    if overlap_1<=0:
                        # modify_cigar_for_overlap(aln1,overlap_2,2)
                        if end-start>500:
                            new_alns.insert(newalni,aln)
                    else:
                        if end-start+overlap_2-overlap_1>200:
                            # modify_cigar_for_overlap(aln1,overlap_1,1)
                            # modify_cigar_for_overlap(aln1,-overlap_2,2)
                            if end-start>500:
                                new_alns.insert(newalni,aln)
                    break
            if newalni==len(new_alns):
                 new_alns.append(aln)
    newalns=[]
    alni=0
    while alni <len(new_alns)-2:
        if (-new_alns[alni][2]+new_alns[alni+1][2]<50):#减去包含在内的
            newalns.append(new_alns[alni])
            alni+=2
        else:
            newalns.append(new_alns[alni])
            alni+=1
    while alni<len(new_alns):
        newalns.append(new_alns[alni])
        alni+=1

    
    return newalns#按start顺序排好的完整contig 并且返回了mapq《59的al