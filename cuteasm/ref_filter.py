#去除变异中与参考基因组本身带来的FP
#输入：检测的变异候选和参考基因组的预输入候选
#输出：不交叠的变异候选
types_to_output = ['DEL', 'INS', 'INV', 'DUP_TAN', 'DUP_INT', 'BND']

def remove_fp(svcandidates, refcandidates):
    svcandidates_no_fp = {svtype: [] for svtype in types_to_output}
    ref_buckets = {}
    base=100

    # 建立桶
    for svtype, refsv_list in refcandidates.items():
        for refsv in refsv_list:
            chrom = refsv[1]
            start = refsv[2]
            bucket = int(int(start) / base)

            if chrom not in ref_buckets:
                ref_buckets[chrom] = {}
            if svtype not in ref_buckets[chrom]:
                ref_buckets[chrom][svtype] = {}
            if bucket not in ref_buckets[chrom][svtype]:
                ref_buckets[chrom][svtype][bucket] = []

            ref_buckets[chrom][svtype][bucket].append(refsv)

    # 检查并筛选
    for svtype, svcand_list in svcandidates.items():
        for svcand in svcand_list:
            chrom = svcand[1]
            start = svcand[2]
            bucket = int(start / base)
            flag = 0

            if chrom in ref_buckets and svtype in ref_buckets[chrom]:
                relevant_buckets = [b for b in ref_buckets[chrom][svtype].keys() if abs(b - bucket) < 1]
                for bucket_start in relevant_buckets:
                    for refsv in ref_buckets[chrom][svtype][bucket_start]:
                        if svtype == 'BND' or svtype == 'DUP_INT':
                            if svcand[3] == refsv[3] and abs(svcand[2]-refsv[2])<50 and abs(svcand[4]-refsv[4])<50:
                                flag = 1
                                break
                        elif svtype == 'INV' or svtype == 'DUP_TAN':
                            if abs(start - refsv[2]) < 50 and abs(svcand[3] - refsv[3]) < 50:
                                flag = 1
                                break
                        else:
                            if abs(start-refsv[2])<50:#如果距离相近
                                flag = 1
                                break
                    if flag == 1:
                        print(svtype, start,svcand[3])
                        break

            if flag == 0:
                svcandidates_no_fp[svtype].append(svcand)

    return svcandidates_no_fp
if __name__ == '__main__':
    svcandidates={'DEL':[['chr1',1234,23]],'INS':[['chr1',12343,78,'ACTTTACTGATCAGTAC']],'BND':[['chr1',12343,'chr2',2342]],'DUP_TAN':[['chr1',12434,12467]],'INV':[['chr1',1234,1245]],'DUP_INT':[]}
    refcandidates={'DEL':[['chr1',123,23]],'INS':[['chr1',12343,78,'ACTTTACTGATCAGTAC']],'BND':[['chr1',12343,'chr2',2342]],'DUP_TAN':[['chr1',12434,12467]],'INV':[['chr1',1234,1245]],'DUP_INT':[]}  

    print(remove_fp(svcandidates, refcandidates))