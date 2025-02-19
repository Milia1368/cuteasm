def find_permutation_cycles(original, permuted):
    n = len(original)
    visited = [False] * n  # 用于跟踪访问过的元素
    cycles = []

    for i in range(n):
        if not visited[i]:
            cycle = []
            current = i
            
            # 追踪循环
            while not visited[current]:
                visited[current] = True
                cycle.append(original[current])
                current = permuted.index(original[current])
            
            # 如果循环不为空，添加到结果中
            if cycle:
                cycles.append(cycle)

    return cycles

# original = [1, 2, 3, 4, 5, 6, 7]
# permuted = [1, 4, 3, 5, 2, 6, 7]

# cycles = find_permutation_cycles(original, permuted)

# # 输出结果
# for cycle in cycles:
#     print('(' + ' '.join(map(str, cycle)) + ')', end='')


import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
def plot(reads,c=None):
    # 示例数据：假设有3条reads
    
    # 配置图形
    plt.figure(figsize=(10, 6))
   
    # 为每条read绘制斜线
    
    first_positions = [item[0] for item in reads]

    # 使用 Counter 统计出现次数
    counter = Counter(first_positions)
    for key, count in counter.items():
        num=0
        output="./plot_"+str(key)+'_'+str(c)+".png"
        new_read=[]
        if count > 3:  # 只处理出现次数大于3的元素
            for read in reads:
                if read[0] == key:  # 判断读取信息中的对应元素是否匹配当前处理的元素
                    new_read.append(read)
                    plt.plot([read[1], read[2]], 
                            [read[3], read[4]], 
                            marker='o', linewidth=2, label=num)
                    num += 1
            if new_read:  # 如果有符合条件的读取信息，进行后续图形设置和保存操作
                min_y = min([read[3] for read in new_read])
                max_y = max([read[4] for read in new_read])
                min_x = min([read[1] for read in new_read])
                max_x = max([read[2] for read in new_read])

                # 设置图形属性
                plt.xlabel('Contig Position')
                plt.ylabel('Reference Position')
                plt.title('Reads Relative Position on Contig vs Reference')
                plt.xlim(min_x - 100, max_x + 100)
                plt.ylim(min_y - 100, max_y + 100)
                plt.grid(True)
                plt.legend()

                # 显示图形
                plt.savefig(output)
                plt.close()
    # # 找到出现次数最多的值
    # most_common = counter.most_common(1)  # 返回出现次数最多的元素及其计数
    # print(most_common)
    # # 输出结果
    # if most_common:
    #     value, count = most_common[0]
    # print(value)
    # for read in reads:
    #     if read[0]!=value:
    #         num+=1
    #         continue
        
    #     new_read.append(read)
    #     plt.plot([read[1], read[2]], 
    #             [read[3], read[4]], 
    #             marker='o', linewidth=2, label=num)
    #     num+=1
    # min_y=min([read[3] for read in new_read])
    # max_y=max([read[4] for read in new_read])
    # min_x=min([read[1] for read in new_read])
    # max_x=max([read[2] for read in new_read])


    # # 设置图形属性
    # plt.xlabel('Contig Position')
    # plt.ylabel('Reference Position')
    # plt.title('Reads Relative Position on Contig vs Reference')
    # plt.xlim(min_x-100, max_x+100)
    # plt.ylim(min_y-100, max_y+100)
    # plt.grid(True)
    # plt.legend()

    # # 显示图形
    # plt.savefig(output)
    # plt.close()

