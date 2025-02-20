"""变异及长度分布图 svim-asm"""
from collections import Counter
import matplotlib
import logging
import random
import math
matplotlib.use('Agg')
import matplotlib.pyplot as plt
types = ['DEL', 'INS', 'INV', 'DUP']
def plot_sv_lengths(deletion_candidates, inversion_candidates, tan_dup_candidates, novel_insertion_candidates, options):
    # len_dict = dict()
    # len_dict["DEL"] = [v.get_source()[2] - v.get_source()[1] for v in deletion_candidates]
    # len_dict["INV"] = [v.get_source()[2] - v.get_source()[1] for v in inversion_candidates]
    # # len_dict["DUP_INT"] = [v.get_destination()[2] - v.get_destination()[1] for v in int_duplication_candidates]
    # len_dict["DUP"] = [v.get_source()[2] - v.get_source()[1] for v in tan_dup_candidates]
    # len_dict["INS"] = [v.get_source()[2] - v.get_source()[1]for v in novel_insertion_candidates]
    len_dict = dict()
    len_dict["DEL"] = [v[3] for v in deletion_candidates]
    len_dict["INV"] = [v[3]-v[2] for v in inversion_candidates]
    # len_dict["DUP_INT"] = [v.get_destination()[2] - v.get_destination()[1] for v in int_duplication_candidates]
    len_dict["DUP"] = [v[3]-v[2] for v in tan_dup_candidates]
    len_dict["INS"] = [v[3] for v in novel_insertion_candidates]
    draw_sv_length_plot(dict_of_lengths=len_dict, output=options.working_dir + "/sv-lengths.png")


def draw_sv_length_plot(dict_of_lengths, output):
    """Makes two stacked bar charts
    Plotting two bar charts of number of SVs by length split by SV type
    Use a consistent colouring scheme for those in "standard_order" to
    make comparison reasonable

    First bar chart is up to 2kb with bins of 10bp
    Second bar chart is up to 20kb, with bins of 100bp
     and uses log scaling on the y-axis
    """
    standard_order = [ 'DEL','INS', 'INV', 'DUP']
    names, lengths = zip(
            *sorted([(svtype, lengths) for svtype, lengths in dict_of_lengths.items()],
                    key=lambda x: standard_order.index(x[0])))
    plt.subplot(2, 1, 1)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 2000, 10)],
             stacked=True,
             histtype='bar',
             label=names)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")

    plt.subplot(2, 1, 2)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 20000, 100)],
             stacked=True,
             histtype='bar',
             label=names,
             log=True)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.tight_layout()
    plt.savefig(output)
    plt.clf()

def plot(reads,c=None):#绘制横坐标为contig位置，纵坐标为reference位置的图形，同一个ref，并去掉极端值
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
import pysam
if __name__ == "__main__":
    # 打开 VCF 文件
    vcf_file = '/data/1/wuxiaomia/Cuteasm/NA24385/ans.vcf'  # 替换为您的 VCF 文件路径
    vcf = pysam.VariantFile(vcf_file,"r")

    # 遍历 VCF 文件中的每一行
    candidates={}
    for type in types:
        candidates[type]=[]
    for record in vcf:
        # 获取每一行的信息
        
        svtype = record.info.get('SVTYPE', None)
        svlength=record.info.get('SVLEN', None)
        candidates[svtype].append(abs(svlength))
    draw_sv_length_plot(candidates,"/data/1/wuxiaomia/Cuteasm/NA24385/sv_length_ans.png")
    # 关闭 VCF 文件
    vcf.close()
