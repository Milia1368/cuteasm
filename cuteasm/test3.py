# import networkx as nx

# # 创建一个空图（有向图或无向图都可以，根据需求）
# G = nx.Graph()

# # 变异信号列表
# sv_signals = {
#     'DEL': [0, 'chr1', 1076341, 238, 'h1tg000217l'],
#     'INS': [1, 'chr1', 58518307, 47, 'h1tg000009l', 'AAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAG'],
#     'INV': [2, 'chr1', 16590490, 16677785, '-+-', 'h1tg000009l'],
#     'DUP_TAN': ['chr1', 124910049, 124914447, 2, True, 'h1tg000090l'],
#     'DUP_INT': [False, False, 'chr1', 143190235, 'chr1', 143184700, 'h1tg000066l', 'dup_int_before'],
#     'BND': [5, 'chr1', 2324343, 'chr11', 29517493, 'h1tg000189l', 'Nor']
# }

# # 创建图并添加节点
# for sv_type, sv_data in sv_signals.items():
#     node_name = f"{sv_type}_{sv_data[1]}_{sv_data[2]}_{sv_data[3]}"  # 可以自定义节点名称
#     G.add_node(node_name, sv_type=sv_type, details=sv_data)

# # 添加变异之间的边（例如，如果它们靠得很近或者存在支持关系）
# # 假设我们使用一些简单的规则来连接变异，如位置接近的变异
# # 例如，DEL和INS如果在相同染色体上并且接近，则连接它们

# for sv_type1, sv_data1 in sv_signals.items():
#     for sv_type2, sv_data2 in sv_signals.items():
#         if sv_type1 != sv_type2 and sv_data1[1] == sv_data2[1]:  # 只连接同一染色体的变异
#             # 变异的位置接近（这里设定阈值为5000bp）
#             if abs(sv_data1[2] - sv_data2[2]) < 5000:  # 起始位置相差小于5000bp
#                 node1 = f"{sv_type1}_{sv_data1[1]}_{sv_data1[2]}_{sv_data1[3]}"
#                 node2 = f"{sv_type2}_{sv_data2[1]}_{sv_data2[2]}_{sv_data2[3]}"
#                 G.add_edge(node1, node2)

# # 查看图中的所有节点和边
# print(f"Graph Nodes: {G.nodes(data=True)}")
# print(f"Graph Edges: {G.edges()}")



import networkx as nx

def build_tandem_duplication_graph(dup_signals, max_gap=5000):
    """
    构建串联重复信号的图结构
    :param dup_signals: Tandem Duplication 信号列表，每个信号是一个字典，包含 chr、start、end、copies 等信息。
    :param max_gap: 最大间隔，若两个信号之间的距离小于 max_gap，则认为它们相关联。
    :return: 一个 NetworkX 图结构
    """
    G = nx.Graph()
    
    # 添加节点
    for i, dup in enumerate(dup_signals):
        G.add_node(i, **dup)  # 每个信号作为一个节点
    
    # 添加边：根据信号的距离和位置关系
    for i in range(len(dup_signals)):
        for j in range(i + 1, len(dup_signals)):
            if dup_signals[i]['chrom'] == dup_signals[j]['chrom']:  # 同染色体
                distance = abs(dup_signals[i]['end'] - dup_signals[j]['start'])
                if distance <= max_gap:  # 距离在允许范围内
                    G.add_edge(i, j)
    
    return G
dup_signals = [
    {'chrom': 'chr1', 'start': 1000, 'end': 2000, 'copies': 2},
    {'chrom': 'chr1', 'start': 2100, 'end': 3100, 'copies': 3},
    {'chrom': 'chr1', 'start': 4000, 'end': 5000, 'copies': 1},
    {'chrom': 'chr1', 'start': 5100, 'end': 6000, 'copies': 4},
]
G = build_tandem_duplication_graph(dup_signals)
# 找到所有连通分量
connected_components = list(nx.connected_components(G))

# 输出连通分量
for component in connected_components:
    print(f"Connected component: {component}")
    # 对每个连通分量的信号进行分析
    signals_in_component = [G.nodes[node] for node in component]
    print(f"Signals in component: {signals_in_component}")
def merge_signals_in_component(component, G):
    """
    合并连通分量中的信号
    :param component: 连通分量（节点集合）
    :param G: 图结构
    :return: 合并后的信号
    """
    signals = [G.nodes[node] for node in component]
    merged_signal = {
        'chrom': signals[0]['chrom'],  # 假设同一连通分量中的信号在同一染色体上
        'start': min(signal['start'] for signal in signals),
        'end': max(signal['end'] for signal in signals),
        'copies': sum(signal['copies'] for signal in signals) / len(signals)  # 平均拷贝数
    }
    return merged_signal
final_results = []
for component in connected_components:
    merged_signal = merge_signals_in_component(component, G)
    final_results.append(merged_signal)

# 输出最终结果
for result in final_results:
    print(result)
