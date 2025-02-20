#负责将combine csv
#按染色体和位置遍历，按照距离聚类SV 将大于两次出现的作为CSV
#candidates=[ins,del,isi]没有道理啊
#按照位置将candidate中的变异聚类 sorted
#多了csv的记录和
#
def merge(i,j,candidate):
    
    pass
def combine_svs(cluster):
    pass
def recombine(candidate,bam):
    min_sv_size=200
    sorted_candidates = sorted(candidate, key=lambda evi: evi.get_key())
    new_candidates=[]
    cluster=[]
    csvs=[]#复杂结构变异记录
    for i in range(len(candidate)-1):
        if sorted_candidates[i].get_key()[0]==sorted_candidates[i+1].get_key()[0]:
            if sorted_candidates[i+1].get_key()[1]-sorted_candidates[i].get_key()[1]<min_sv_size:
                if sorted_candidates[i].get_key()[2]==sorted_candidates[i+1].get_key()[2]:
                    candidate_i=merge(i,i+1,sorted_candidates)
                else:
                    if len(cluster)==0:
                        cluster.append(sorted_candidates[i])
                    else:
                        if (sorted_candidates[i].get_key()[1]-cluster[-1].get_key()[1]<min_sv_size )and ( sorted_candidates[i].get_key()[0]==cluster[-1].get_key()[0]) :
                            cluster.append(sorted_candidates[i])
                        else:
                            csvs.append(combine_svs(cluster))
                            cluster=[]
        else:
            new_candidates.append(sorted_candidates[i])
    if len(cluster) > 0:
        csvs.append(combine_svs(cluster))


    return new_candidates,csvs
    