import ast
from collections import defaultdict
"""先找到染色体上重复的区域，再把该区域对应的变异信号给过滤"""


def find_dup_region(dup_candidates, Max_gap_tolerance=500, threshold=50000):
    same_chr_reads = {}
    same_chr_blocks = {}
    ref_start = 2
    ref_end = 3
    chr_num = 1

    # 将重复候选按染色体分组
    for read in dup_candidates:
        same_chr_reads.setdefault(read[chr_num], []).append([read[ref_start], read[ref_end],0])

    # 对每个染色体的起始位置排序
    for key in same_chr_reads:
        same_chr_reads[key].sort(key=lambda x: (x[0], x[1]))

        # 对主要的染色体区域进行分块
        ref_blocks = same_chr_reads[key]
        new_ref_blocks = [ref_blocks[0]]

        for i in range(len(ref_blocks) - 1):
            if ref_blocks[i + 1][0] - new_ref_blocks[-1][1] <= Max_gap_tolerance:
                new_ref_blocks[-1][1] = max(ref_blocks[i + 1][1], new_ref_blocks[-1][1])
                new_ref_blocks[-1][-1] = new_ref_blocks[-1][-1]+1
            else:
                new_ref_blocks.append(ref_blocks[i + 1])

        same_chr_blocks[key] = [block for block in new_ref_blocks if block[-1]>=2]
    # 扩展块的范围
    for key in same_chr_blocks:
        for block in same_chr_blocks[key]:
            block[0] = max(0,block[0]-threshold)
            block[1] += threshold
    for key in same_chr_blocks:
        for block in same_chr_blocks[key]:
            print(key,block)
    return same_chr_blocks



def remove_can_in_dupregion(candidates, chrom_dup_blocks,mode=0):

    # 过滤候选
    fil_can = []
    chr_num = 1
    start_num = 2
    len_num = 3
    can_list=[]
    for can in candidates:
        flag=0
        if can[chr_num] in chrom_dup_blocks:
            for block in chrom_dup_blocks[can[chr_num]]:
                if mode==0:

                    if block[1] >= can[start_num] >= block[0] or block[0] <= can[start_num] + can[len_num] <= block[1]:
                        fil_can.append(can)
                        flag=1
                        break  # 找到匹配后跳出内层循环
                else:#ins
                    if block[1] >= can[start_num] >= block[0] :
                        fil_can.append(can)
                        flag=1
                        break  # 找到匹配后跳出内层循环
        if not flag:
            can_list.append(can)
    # 输出结果
    print(len(fil_can))
    # for can in can_list:
    #     print(can)   # 输出未匹配的候选
    return can_list

# 调用函数
if __name__=='__main__':
    data_dup = """
    [4, 'chr1', 73129297, 73939836, 58356713, 'h1tg000009l', 3]
    [4, 'chr1', 121745679, 121746845, 107327522, 'h1tg000009l', 3]
    """
    data = """
    [0, 'chr1', 1605539, 40, 'h1tg000189l']
    [0, 'chr1', 4632566, 63, 'h1tg000146l']
    [0, 'chr1', 6004963, 238, 'h1tg000146l']
    [0, 'chr1', 31595303, 44, 'h1tg000009l']"""
    dup_candidates= [eval(line.strip()) for line in data_dup.strip().split('\n')]
    candidates = [ast.literal_eval(line.strip()) for line in data.strip().split('\n')]
    chrom_dup_blocks = find_dup_region(dup_candidates)
    can_list = remove_can_in_dupregion(candidates, chrom_dup_blocks)