# cuteasm 检测拼接的基因组结构变异
## 文件结构
--cuteasm/anocode 无用的code，一些没有正式完成的想法
--cuteasm/overlap_collect.py、cuteasm/plot.py 画图，用于调试和统计
--cuteasm/cuteasm 主函数
    --cuteasm/cute_input_parsing.py、cuteasm/cute_out_format.py 分别为参数部分和格式化输出模块
    --cuteasm/cute_candidate.py 变异信号对象定义，参考svim-asm进行了微调
    --cuteasm/cute_collect.py 调用intra和inter，处理染色体上的比对信息
        --cuteasm/cute_intra.py、cuteasm/cute_inter.py分别处理cigar和split信号
    --cuteasm/cute_genotyping.py 基因分型模块，参考svim-asm，并进行优化
    --cuteasm/realign.py  重比对模块，对应复杂结构变异部分
    --cuteasm/ref_filter.py  参考基因组假阳性过滤部分
## 创新点
1.比对信息过滤
2.复杂结构变异重比对
3.假阳性过滤
4.基因分型
5、聚类和split信号检测变异方面有优化和改进
## 结果
与cuteSV、SVIM-asm、dipcall、SVanalyzer在hg002数据集上
recall、GT 领先6%~10%
基因分型 95%~99%
