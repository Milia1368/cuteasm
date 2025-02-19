"""统计read比对段前后重叠分布等信息"""
#输入列表的supl 
#输出差距分布柱形图
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.ticker import MaxNLocator

def plot_deta1(detal_lengths):
    
    # 计算差值的范围和数量
    max_distance = max(detal_lengths)
    min_distance = min(detal_lengths)
    bins = np.arange(min_distance, max_distance + 1, 5000)  # 设置差值块的范围和大小

    # 计算每个差值块中的数量
    counts, _ = np.histogram(detal_lengths, bins=bins)

    # 计算统计值
    mean_distance = np.mean(detal_lengths)
    std_distance = np.std(detal_lengths)

    # 绘制柱形图
    plt.figure(figsize=(10, 6))
    plt.bar(bins[:-1], counts, width=np.diff(bins), align='edge', color='skyblue', edgecolor='black', alpha=0.7)

    # 绘制分布曲线
    # 计算分布曲线的 x 和 y 值
    x = np.linspace(min_distance, max_distance, 100)
    y = norm.pdf(x, mean_distance, std_distance) * sum(counts) * np.diff(bins)[0]  # 归一化

    plt.plot(x, y, color='red', linewidth=2, label='Distribution Curve')

    # 设置简约风格
    plt.xlabel('Difference Blocks', fontsize=12)
    plt.ylabel('Count of Reads in Each Block', fontsize=12)
    plt.title('Count of Reads in Difference Blocks with Distribution Curve', fontsize=14)
    plt.xticks(bins, fontsize=10)  # 设置 x 轴刻度
    plt.yticks(fontsize=10)  # 设置 y 轴刻度为整数
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # 标注统计值
    plt.axvline(mean_distance, color='green', linestyle='--', label=f'Mean: {mean_distance:.2f}')
    plt.axvline(mean_distance + std_distance, color='orange', linestyle='--', label=f'Std Dev: {std_distance:.2f}')
    plt.axvline(mean_distance - std_distance, color='orange', linestyle='--')

    # 添加图例
    plt.legend()

    # 确保 y 轴为整数
    plt.gca().yaxis.get_major_locator().set_params(integer=True)

    plt.tight_layout()  # 自动调整子图参数
    plt.savefig('test.png', dpi=300)  # 保存图片
    # plt.show()  # 显示图片
def collect(group,detal_lengths,readname=None,d_list=None):
    
    # 提取所有 read 的结束位置并计算相邻 read 之间的差值
    deta_list=[]
    end_pos = group[0][2]  # 初始化第一个 read 的结束位置
    for read in group[1:]:  # 从第二个 read 开始
        a=read[1] - end_pos
        deta_list.append(a)  # 计算差值
        # if abs(a)>41000:
        #     print(a)
        end_pos = read[2]  # 更新结束位置
    i=0
    for i in range(len(deta_list)-1):
        if deta_list[i]<0 and deta_list[i+1]<0:
            a=group[i+1][2]-group[i+1][1]
            b=deta_list[i]
            c=deta_list[i+1]
            d=abs(b+c)/a*100
            d_list.append(d)
            detal_lengths.append([a,b,c,d,group[i+1][0],group[i+1][1],group[i+1][3],readname])
def plot_deta(detal_lengths):


   # 计算均值和标准差
    mean_distance = np.mean(detal_lengths)
    std_distance = np.std(detal_lengths)

    # 计算差值的范围和数量
    max_distance = max(detal_lengths)
    min_distance = min(detal_lengths)
    bins = np.arange(min_distance, max_distance + 1, 20000)  # 设置差值块的范围和大小

    # 计算每个差值块中的数量
    counts, _ = np.histogram(detal_lengths, bins=bins)

    # 计算总数并归一化为百分比
    total_counts = sum(counts)
    percentages = (counts / total_counts) * 100  # 计算百分比

    # 绘制完整分布的柱形图
    plt.figure(figsize=(10, 6))
    plt.bar(bins[:-1], percentages, width=np.diff(bins), align='edge', color='skyblue', edgecolor='black', alpha=0.7)

    # 绘制完整分布的分布曲线
    x = np.linspace(min_distance, max_distance, 100)
    y = norm.pdf(x, mean_distance, std_distance) * total_counts * np.diff(bins)[0] / total_counts * 100  # 归一化为百分比
    plt.plot(x, y, color='red', linewidth=2, label='Full Distribution Curve')

    # 设置简约风格
    plt.xlabel('Difference Blocks', fontsize=12)
    plt.ylabel('Percentage of Reads in Each Block (%)', fontsize=12)
    plt.title('Percentage of Reads in Difference Blocks with Full Distribution Curve', fontsize=14)

    # 自动调整横坐标刻度的显示个数
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True, nbins='auto'))  # 自动选择刻度数量
    plt.xticks(fontsize=10)  # 设置 x 轴刻度
    plt.yticks(fontsize=10)  # 设置 y 轴刻度
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # 标注统计值
    plt.axvline(mean_distance, color='green', linestyle='--', label=f'Mean: {mean_distance:.2f}')
    plt.axvline(mean_distance + std_distance, color='orange', linestyle='--', label=f'Std Dev: {std_distance:.2f}')
    plt.axvline(mean_distance - std_distance, color='orange', linestyle='--')

    # 添加图例
    plt.legend()

    # 设置坐标轴范围以包含负值
    plt.xlim(min_distance - 10, max_distance + 10)  # 根据需要调整范围

    # 确保 y 轴为整数
    plt.gca().yaxis.get_major_locator().set_params(integer=True)

    plt.tight_layout()  # 自动调整子图参数
    plt.savefig('full_distribution.png', dpi=300)  # 保存完整分布图
    # plt.show()  # 显示完整分布图

    #  确定极端值的阈值（均值 ± 2 标准差）
    lower_bound = mean_distance - 2 * std_distance
    upper_bound = mean_distance + 2 * std_distance

    # 过滤掉极端值
    filtered_lengths = [x for x in detal_lengths if lower_bound <= x <= upper_bound]

    # 绘制过滤后的数据的分布图
    plt.figure(figsize=(10, 6))

    # 设置更小的分桶间隔
    bins = np.linspace(min(filtered_lengths), max(filtered_lengths), 100)  

    # 绘制过滤后的直方图，归一化为百分比
    counts, _ = np.histogram(filtered_lengths, bins=bins)
    percentages = (counts / sum(counts)) * 100  # 计算百分比

    # 绘制直方图
    plt.bar(bins[:-1], percentages, width=np.diff(bins), align='edge', alpha=0.6, color='skyblue', edgecolor='black')

    # 绘制过滤后的分布曲线
    filtered_mean = np.mean(filtered_lengths)
    filtered_std = np.std(filtered_lengths)
    x_filtered = np.linspace(min(filtered_lengths), max(filtered_lengths), 100)
    y_filtered = norm.pdf(x_filtered, filtered_mean, filtered_std) * sum(counts) * np.diff(bins)[0] / sum(counts) * 100  # 归一化为百分比

    

    # 设置图形属性
    plt.xlabel('Difference Blocks', fontsize=12)
    plt.ylabel('Percentage (%)', fontsize=12)
    plt.title('Filtered Distribution of Difference Blocks', fontsize=14)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    # 标注统计值
    plt.axvline(filtered_mean, color='green', linestyle='--', label=f'Mean: {filtered_mean:.2f}')
    plt.axvline(filtered_mean + filtered_std, color='orange', linestyle='--', label=f'Std Dev: {filtered_std:.2f}')
    plt.axvline(filtered_mean - filtered_std, color='orange', linestyle='--')
    plt.plot(x_filtered, y_filtered, color='red', linewidth=2, label='Filtered Distribution Curve')
    # 添加图例
    plt.legend()

    plt.tight_layout()  # 自动调整子图参数
    plt.savefig('filtered_distribution.png', dpi=300)  # 保存过滤后的分布图
    plt.show()  # 显示过滤后的分布图
def plot_d(percentages):
    # 绘制数据分布的直方图
    plt.figure(figsize=(10, 6))
    plt.hist(percentages, bins=50, color='skyblue', edgecolor='black', alpha=0.7)

    # 设置图形属性
    plt.xlabel('Percentage (%)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title('Distribution of Percentages', fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()  # 自动调整子图参数
    plt.savefig('d.png')  # 显示图片