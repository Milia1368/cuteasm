import os
import logging
import re

from collections import defaultdict
from math import pow, sqrt
import time
from statistics import mean, stdev
from edlib import align
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

from cute_candidate import *

def form_partitions(sv_candidates_with_haplotype, max_distance):
    """Form partitions of signatures using mean distance."""
    sorted_candidates_with_haplotype = sorted(sv_candidates_with_haplotype, key=lambda evi: evi[1].get_key())
    partitions = []
    current_partition = []#可能是按照参考染色体、位置排，将可以放到一块的的hap1、2划分到一个patition
    for haplotype, candidate in sorted_candidates_with_haplotype:
        if len(current_partition) > 0:
            candidate_key = candidate.get_key()
            len1=candidate_key[2]-candidate_key[1]
            last_candidate_key = current_partition[-1][1].get_key()
            len2=last_candidate_key[2]-last_candidate_key[1]
            thred=min(int(min(len1,len2)*0.2),max_distance)
            if last_candidate_key[0] != candidate_key[0] or \
               abs(last_candidate_key[1] - candidate_key[1])>thred or \
               abs(last_candidate_key[2] - candidate_key[2]) > thred:
                partitions.append(current_partition[:])
                current_partition = []
        current_partition.append((haplotype, candidate))
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    # par_len=[len(p) for p in partitions]
    # print(par_len)
    return partitions
def form_partitions_ins(sv_candidates_with_haplotype, max_distance):
    """Form partitions of signatures using mean distance."""
    sorted_candidates_with_haplotype = sorted(sv_candidates_with_haplotype, key=lambda evi: evi[1].get_key())
    partitions = []
    current_partition = []#可能是按照参考染色体、位置排，将可以放到一块的的hap1、2划分到一个patition
    for haplotype, candidate in sorted_candidates_with_haplotype:
        if len(current_partition) > 0:
            candidate_key = candidate.get_key()
            len1=candidate_key[2]
            last_candidate_key = current_partition[-1][1].get_key()
            len2=last_candidate_key[2]
            thred=min(int(min(len1,len2)*0.2),max_distance)
            if last_candidate_key[0] != candidate_key[0] or \
               abs(last_candidate_key[1] - candidate_key[1])>thred or \
               abs(last_candidate_key[2] - candidate_key[2]) > thred:
                partitions.append(current_partition[:])
                current_partition = []
        current_partition.append((haplotype, candidate))
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    # par_len=[len(p) for p in partitions]
    # print(par_len)
    return partitions
def form_partitions_bnd(sv_candidates_with_haplotype, max_distance):
    """Form partitions of signatures using mean distance."""
    sorted_candidates_with_haplotype = sorted(sv_candidates_with_haplotype, key=lambda evi: evi[1].get_key())
    partitions = []
    current_partition = []#可能是按照参考染色体、位置排，将可以放到一块的的hap1、2划分到一个patition
    for haplotype, candidate in sorted_candidates_with_haplotype:
        if len(current_partition) > 0:
            candidate_key = candidate.get_key()
            last_candidate_key = current_partition[-1][1].get_key()
            if last_candidate_key[0] != candidate_key[0] or \
               abs(last_candidate_key[1] - candidate_key[1])>200 or \
               abs(last_candidate_key[2] - candidate_key[2]) > 200:
                partitions.append(current_partition[:])
                current_partition = []
        current_partition.append((haplotype, candidate))
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    # par_len=[len(p) for p in partitions]
    # print(par_len)
    return partitions
def compute_distance(candidate_with_haplotype1, candidate_with_haplotype2, reference):
    haplotype1, candidate1 = candidate_with_haplotype1
    haplotype2, candidate2 = candidate_with_haplotype2
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    if haplotype1 == haplotype2:
        return 1000000000

    if candidate1.type == "DEL":
        region_chr = candidate1.source_contig
        chr_length = reference.get_reference_length(region_chr)
        region_start = max(0, min(candidate1.source_start, candidate2.source_start) - 100)
        region_end = min(chr_length, max(candidate1.source_end, candidate2.source_end) + 100)
        # if region_start==248637623 or candidate1.source_end==248637623:#248387328
        #     pass
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.source_start).upper() + reference.fetch(region_chr, candidate1.source_end, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.source_start).upper() + reference.fetch(region_chr, candidate2.source_end, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "INV":
        region_chr = candidate1.source_contig
        chr_length = reference.get_reference_length(region_chr)
        region_start = max(0, min(candidate1.source_start, candidate2.source_start) - 100)
        region_end = min(chr_length, max(candidate1.source_end, candidate2.source_end) + 100)
        inverted_seq1 = "".join(complement.get(base.upper(), base.upper()) for base in reversed(reference.fetch(region_chr, candidate1.source_start, candidate1.source_end).upper()))
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.source_start).upper() + \
                     inverted_seq1 + \
                     reference.fetch(region_chr, candidate1.source_end, region_end).upper()
        inverted_seq2 = "".join(complement.get(base.upper(), base.upper()) for base in reversed(reference.fetch(region_chr, candidate2.source_start, candidate2.source_end).upper()))
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.source_start).upper() + \
                     inverted_seq2 + \
                     reference.fetch(region_chr, candidate2.source_end, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "INS":
        region_chr = candidate1.source_contig
        chr_length = reference.get_reference_length(region_chr)
        region_start = max(0, min(candidate1.source_start, candidate2.source_start) - 100)
        region_end = min(chr_length, max(candidate1.source_start, candidate2.source_start) + 100)
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.source_start).upper() + \
                     candidate1.sequence + \
                     reference.fetch(region_chr, candidate1.source_start, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.source_start).upper() + \
                     candidate2.sequence + \
                     reference.fetch(region_chr, candidate2.source_start, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "DUP_TAN":
        region_chr = candidate1.source_contig
        chr_length = reference.get_reference_length(region_chr)
        region_start = max(0, min(candidate1.source_start, candidate2.source_start) - 100)
        region_end = min(chr_length, max(candidate1.source_end, candidate2.source_end) + 100)
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.source_start).upper() + \
                     reference.fetch(region_chr, candidate1.source_start, candidate1.source_end).upper() * (candidate1.copies + 1) + \
                     reference.fetch(region_chr, candidate1.source_end, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.source_start).upper() + \
                     reference.fetch(region_chr, candidate2.source_start, candidate2.source_end).upper() * (candidate2.copies + 1) + \
                     reference.fetch(region_chr, candidate2.source_end, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "DUP_INT":
        if candidate1.source_contig != candidate2.source_contig or candidate1.reverse != candidate2.reverse:
            return 1000000000

        return abs((candidate1.get_source()[2]- candidate1.get_source()[1])-(candidate2.get_source()[2]- candidate2.get_source()[1]))

    return editDistance


def span_position_distance_breakends(candidate1, candidate2):
    candidate1_hap, candidate1_pos1, candidate1_dir1, candidate1_pos2, candidate1_dir2 = candidate1
    candidate2_hap, candidate2_pos1, candidate2_dir1, candidate2_pos2, candidate2_dir2 = candidate2
    if candidate1_hap != candidate2_hap:
        if candidate1_dir1 == candidate2_dir1 and candidate1_dir2 == candidate2_dir2:
            dist1 = abs(candidate1_pos1 - candidate2_pos1)
            dist2 = abs(candidate1_pos2 - candidate2_pos2)
            position_distance = (dist1 + dist2) / 3000
        else:
            position_distance = 99999
    else:
        position_distance = 99999
    return position_distance


def pair_haplotypes(partitions, reference, edit_distance_threshold = 10):
    clusters_final = []
    for partition in partitions:
        if len(partition) < 2:
            new_clusters = [partition]
        #Ignore very large partitions because they tend to be in difficult regions
        elif len(partition) > 10:
            logging.debug("Ignored partition of size {0} and type {1}: {2}".format(len(partition), partition[0][1].get_key()[0], ",".join(["{0}:{1}".format(sig[1].get_key()[1], sig[1].get_key()[2]) for sig in partition])))
            continue
        else:
            distances = []
            for i in range(len(partition)-1):
                for j in range(i+1, len(partition)):
                    distances.append(compute_distance(partition[i], partition[j], reference))
            Z = linkage(np.array(distances), method = "complete")#！！！？？？应该是使用编辑距离合并hap1/2
            cluster_indices = list(fcluster(Z, edit_distance_threshold, criterion='distance'))
            new_clusters = [[] for i in range(max(cluster_indices))]
            for candidate_index, cluster_index in enumerate(cluster_indices):
                new_clusters[cluster_index-1].append(partition[candidate_index])
        clusters_final.extend(new_clusters)
    return clusters_final


def pair_haplotypes_breakends(partitions, span_position_distance_threshold = 0.3):
    """Finds clusters in partitions using span-position distance and hierarchical clustering. 
    Assumes that all signatures in the given partition are of the same type and on the same contig"""
    clusters_final = []
    for partition in partitions:
        if len(partition) < 2:
            new_clusters = [partition]
        #Ignore very large partitions because they tend to be in difficult regions
        elif len(partition) > 10:
            continue
        else:
            data = np.array( [[haplotype, candidate.get_source()[1], 1 , candidate.get_destination()[1], 1 ] for (haplotype, candidate) in partition])
            Z = linkage(data, method = "complete", metric = span_position_distance_breakends)
            cluster_indices = list(fcluster(Z, span_position_distance_threshold, criterion='distance'))
            new_clusters = [[] for i in range(max(cluster_indices))]
            for candidate_index, cluster_index in enumerate(cluster_indices):
                new_clusters[cluster_index-1].append(partition[candidate_index])
        clusters_final.extend(new_clusters)
    return clusters_final


def pair_candidates(sv_candidates1, sv_candidates2, reference, bam, options):
    #sv_can1\2分别是两条bam信息
    #分别对将信息提取出来
    deletion_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "DEL"]
    insertion_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "INS"]
    inversion_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "INV"]
    tandem_duplication_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "DUP_TAN"]
    breakend_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "BND"]
    interspersed_duplication_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "DUP_INT"]

    deletion_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "DEL"]
    insertion_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "INS"]
    inversion_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "INV"]
    tandem_duplication_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "DUP_TAN"]
    breakend_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "BND"]
    interspersed_duplication_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "DUP_INT"]

    paired_candidates = []
    #DELETIONS
    logging.info("Pairing {0} deletions...".format(len(deletion_candidates1) + len(deletion_candidates2)))
    partitions = form_partitions(deletion_candidates1 + deletion_candidates2,options.partition_max_distance )#
    # with open("./partion_del1.txt",'w') as file_out:
    #     for part in partitions:
    #         part_list=[]
    #         for p in part:
    #             part_list.append([p[0],p[1].get_key()])
    #         print(part_list,file=file_out)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateDeletion(candidate.source_contig, 
                                                       candidate.source_start, 
                                                       candidate.length, 
                                                       candidate.reads,
                                                       genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateDeletion(candidate.source_contig, 
                                                       candidate.source_start, 
                                                       candidate.length,
                                                       reads,
                                                       genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))
    
    #INVERSIONS
    logging.info("Pairing {0} inversions...".format(len(inversion_candidates1) + len(inversion_candidates2)))
    partitions = form_partitions(inversion_candidates1 + inversion_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateInversion(candidate.source_contig, 
                                                        candidate.source_start, 
                                                        candidate.source_end, 
                                                        candidate.layer,
                                                        candidate.reads, 
                                                        bam, 
                                                        genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            # complete = cluster[0][1].complete or cluster[1][1].complete
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateInversion(candidate.source_contig, 
                                                        candidate.source_start, 
                                                        candidate.source_end, 
                                                        candidate.layer,
                                                        reads,  
                                                        bam, 
                                                        genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    #INSERTIONS
    logging.info("Pairing {0} insertions...".format(len(insertion_candidates1) + len(insertion_candidates2)))
    partitions = form_partitions_ins(insertion_candidates1 + insertion_candidates2, options.partition_max_distance)
    # with open("./partion_ins2_modins.txt",'w') as file_out:
    #     for part in partitions:
    #         part_list=[]
    #         for p in part:
    #             part_list.append([p[0],p[1].get_key()])
    #         print(part_list,file=file_out)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    # with open("./partion_cluster_ins2_modins.txt",'w') as file_out:
    #     for part in clusters:
    #         part_list=[]
    #         for p in part:
    #             part_list.append([p[0],p[1].get_key()])
    #         print(part_list,file=file_out)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateInsertion(candidate.source_contig,candidate.source_start,candidate.length,candidate.reads,candidate.sequence,genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateInsertion(candidate.source_contig,candidate.source_start,candidate.length,reads,candidate.sequence,genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    #TANDEM DUPLICATIONS
    logging.info("Pairing {0} tandem duplications...".format(len(tandem_duplication_candidates1) + len(tandem_duplication_candidates2)))
    partitions = form_partitions(tandem_duplication_candidates1 + tandem_duplication_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateDuplicationTandem(candidate.source_contig, 
                                                                candidate.source_start, 
                                                                candidate.source_end, 
                                                                candidate.copies, 
                                                                candidate.fully_covered, 
                                                                candidate.reads,
                                                                bam, 
                                                                genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            fully_covered = cluster[0][1].fully_covered or cluster[1][1].fully_covered
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateDuplicationTandem(candidate.source_contig, 
                                                                candidate.source_start, 
                                                                candidate.source_end, 
                                                                round(mean([cluster[0][1].copies, cluster[1][1].copies])),
                                                                fully_covered,
                                                                reads, 
                                                                bam, 
                                                                genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))
    

     #INTERSPERSED DUPLICATIONS
    logging.info("Pairing {0} interspersed duplications...".format(len(interspersed_duplication_candidates1) + len(interspersed_duplication_candidates2)))
    partitions = form_partitions(interspersed_duplication_candidates1 + interspersed_duplication_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateDuplicationInterspersed(candidate.source_contig, 
                                                                      candidate.source_start, 
                                                                      candidate.source_end, 
                                                                      candidate.dest_contig, 
                                                                      candidate.dest_start, 
                                                                      candidate.dest_end,
                                                                      candidate.reads,
                                                                      bam, 
                                                                      candidate.cutpaste,
                                                                      genotype,candidate.reverse,candidate.support))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1].reads + cluster[1][1].reads
            cutpaste = cluster[0][1].cutpaste or cluster[1][1].cutpaste
            genotype = "1/1"
            paired_candidates.append(CandidateDuplicationInterspersed(candidate.source_contig, 
                                                                      candidate.source_start, 
                                                                      candidate.source_end, 
                                                                      candidate.dest_contig, 
                                                                      candidate.dest_start, 
                                                                      candidate.dest_end,
                                                                      reads,
                                                                      bam, 
                                                                      cutpaste,
                                                                      genotype,candidate.reverse,candidate.support))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    # #INTERSPERSED DUPLICATIONS
    # logging.info("Pairing {0} interspersed duplications...".format(len(interspersed_duplication_candidates1) + len(interspersed_duplication_candidates2)))
    # partitions = form_partitions_bnd(interspersed_duplication_candidates1 + interspersed_duplication_candidates2, options.partition_max_distance)
    # clusters = pair_haplotypes_breakends(partitions)
    # for cluster in clusters:
    #     if len(cluster) == 1:
    #         candidate = cluster[0][1]
    #         genotype = "1/0" if cluster[0][0] == 1 else "0/1"
    #         paired_candidates.append(CandidateDuplicationInterspersed(candidate.source_contig, 
    #                                                     candidate.source_start, 
    #                                                     candidate.dest_contig, 
    #                                                     candidate.dest_start, 
    #                                                     candidate.layer,
    #                                                     candidate.reads,
    #                                                     bam,
    #                                                     genotype))
    #     elif len(cluster) == 2:
    #         candidate = cluster[0][1]
    #         reads = cluster[0][1].reads + cluster[1][1].reads
    #         genotype = "1/1"
    #         paired_candidates.append(CandidateDuplicationInterspersed(candidate.source_contig, 
    #                                                     candidate.source_start, 
    #                                                     candidate.dest_contig, 
    #                                                     candidate.dest_start, 
    #                                                     candidate.layer,
    #                                                     reads,
    #                                                     bam,
    #                                                     genotype))
    #     else:
    #         logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    #BREAKENDS
    logging.info("Pairing {0} breakends...".format(len(breakend_candidates1) + len(breakend_candidates2)))
    partitions = form_partitions_bnd(breakend_candidates1 + breakend_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes_breakends(partitions)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateBreakend(candidate.source_contig, 
                                                        candidate.source_start, 
                                                        candidate.dest_contig, 
                                                        candidate.dest_start, 
                                                        candidate.layer,
                                                        candidate.reads,
                                                        bam,
                                                        candidate.length,
                                                        candidate.forword,
                                                        candidate.backword,
                                                        genotype))

        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateBreakend(candidate.source_contig, 
                                                        candidate.source_start, 
                                                        candidate.dest_contig, 
                                                        candidate.dest_start, 
                                                        candidate.layer,
                                                        reads,
                                                        bam,
                                                        candidate.length,
                                                        candidate.forword,
                                                        candidate.backword,
                                                        genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))
    return paired_candidates


# def form_partitions(sv_candidates_with_haplotype, max_distance):
#     """Form partitions of signatures using mean distance."""
#     sorted_candidates_with_haplotype = sorted(sv_candidates_with_haplotype, key=lambda evi: evi[1].get_key())
#     partitions = []
#     current_partition = []#可能是按照参考染色体、位置排，将可以放到一块的的hap1、2划分到一个patition
#     _prtions=[]
#     for haplotype, candidate in sorted_candidates_with_haplotype:
#         if len(current_partition) > 0:
#             candidate_key = candidate.get_key()
#             len1=candidate_key[2]-candidate_key[1]
#             last_candidate_key = current_partition[-1][1].get_key()
#             len2=last_candidate_key[2]-last_candidate_key[1]
#             thred=min(int(min(len1,len2)*0.2),max_distance)
#             if last_candidate_key[0] != candidate_key[0] or \
#                abs(last_candidate_key[1] - candidate_key[1])>thred or \
#                abs(last_candidate_key[2] - candidate_key[2]) > thred:
#                 partitions.append(current_partition[:])
#                 cluster_prtions(current_partition)
#                 current_partition = []
#         current_partition.append((haplotype, candidate))
#     if len(current_partition) > 0:
#         partitions.append(current_partition[:])
#     return partitions

def pair_haplotypes1(partitions, reference, edit_distance_threshold = 10):
    clusters_final = []
    for partition in partitions:
        if len(partition) < 2:
            new_clusters = [partition]
        #Ignore very large partitions because they tend to be in difficult regions
        elif len(partition) > 10:
            logging.debug("Ignored partition of size {0} and type {1}: {2}".format(len(partition), partition[0][1][0], ",".join(["{0}:{1}".format(sig[1][1], sig[1][2]) for sig in partition])))
            continue
        else:
            distances = []
            for i in range(len(partition)-1):
                for j in range(i+1, len(partition)):
                    distances.append(compute_distance(partition[i], partition[j], reference))
            Z = linkage(np.array(distances), method = "complete")
            cluster_indices = list(fcluster(Z, edit_distance_threshold, criterion='distance'))
            new_clusters = [[] for i in range(max(cluster_indices))]
            for candidate_index, cluster_index in enumerate(cluster_indices):
                new_clusters[cluster_index-1].append(partition[candidate_index])
        clusters_final.extend(new_clusters)
    return clusters_final


#实现基因分型 基于距离聚类基于序列相似度分类
            # 0    1       2     3         4            5           6               7       8               9
# del:      [0, 'chr1', 1076341,    238,        'h1tg000217l']
# ins:      [1, 'chr1', 58518307,    47,       'h1tg000009l', 'AAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAG']
# INV：     [2, 'chr1', 16590490,  16677785,    '-+-',        'h1tg000009l']
# DuP_tan:  [3，'chr1', 124910049, 124914447,    2,             True,           'h1tg000090l']
# dup_int:  [4， False,  False,     'chr1',      143190235,      'chr1',        143184700, 'h1tg000066l', 'dup_int_before']
# BND：     [5, 'chr1', 2324343,   'chr11',     29517493,       'h1tg000189l',       'Nor']


def form_partitions1(sv_candidates_with_haplotype, max_distance):
    """Form partitions of signatures using mean distance."""
    sorted_candidates_with_haplotype = sorted(sv_candidates_with_haplotype, key=lambda evi: evi[1].get_key())
    partitions = []
    current_partition = []#可能是按照参考染色体、位置排，将可以放到一块的的hap1、2划分到一个patition
    for haplotype, candidate in sorted_candidates_with_haplotype:
        if len(current_partition) > 0:
            candidate_key = candidate.get_key()
            len1=candidate_key[2]-candidate_key[1]
            last_candidate_key = current_partition[-1][1].get_key()
            len2=last_candidate_key[2]-last_candidate_key[1]
            thred=min(int(min(len1,len2)*0.2),max_distance)
            if last_candidate_key[0] != candidate_key[0] or \
               abs(last_candidate_key[1] - candidate_key[1])>thred or \
               abs(last_candidate_key[2] - candidate_key[2]) > thred:
                partitions.append(current_partition[:])
                current_partition = []
        current_partition.append((haplotype, candidate))
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    par_len=[len(p) for p in partitions]
    print(par_len)
    return partitions


def pair_candidates1(sv_candidates1, sv_candidates2, reference, bam, options):
    deletion_candidates1 = [(1, cand) for cand in sv_candidates1 if cand[0] == "DEL"]
    insertion_candidates1 = [(1, cand) for cand in sv_candidates1 if cand[0] == "INS"]
    inversion_candidates1 = [(1, cand) for cand in sv_candidates1 if cand[0] == "INV"]
    tandem_duplication_candidates1 = [(1, cand) for cand in sv_candidates1 if cand[0] == "DUP_TAN"]
    breakend_candidates1 = [(1, cand) for cand in sv_candidates1 if cand[0] == "BND"]
    interspersed_duplication_candidates1 = [(1, cand) for cand in sv_candidates1 if cand[0] == "DUP_INT"]

    deletion_candidates2 = [(2, cand) for cand in sv_candidates2 if cand[0] == "DEL"]
    insertion_candidates2 = [(2, cand) for cand in sv_candidates2 if cand[0] == "INS"]
    inversion_candidates2 = [(2, cand) for cand in sv_candidates2 if cand[0] == "INV"]
    tandem_duplication_candidates2 = [(2, cand) for cand in sv_candidates2 if cand[0] == "DUP_TAN"]
    breakend_candidates2 = [(2, cand) for cand in sv_candidates2 if cand[0] == "BND"]
    interspersed_duplication_candidates2 = [(2, cand) for cand in sv_candidates2 if cand[0] == "DUP_INT"]

    paired_candidates = []
    # DELETIONS
    logging.info("Pairing {0} deletions...".format(len(deletion_candidates1) + len(deletion_candidates2)))
    partitions = form_partitions(deletion_candidates1 + deletion_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(["DEL", candidate[1], candidate[2], candidate[3], genotype])
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1][4] + cluster[1][1][4]
            genotype = "1/1"
            paired_candidates.append(["DEL", candidate[1], candidate[2], candidate[3], reads, genotype])
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    # INVERSIONS
    logging.info("Pairing {0} inversions...".format(len(inversion_candidates1) + len(inversion_candidates2)))
    partitions = form_partitions(inversion_candidates1 + inversion_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(["INV", candidate[1], candidate[2], candidate[3], candidate[4], genotype])
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1][5] + cluster[1][1][5]
            genotype = "1/1"
            paired_candidates.append(["INV", candidate[1], candidate[2], candidate[3], candidate[4], reads, genotype])
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    # INSERTIONS
    logging.info("Pairing {0} insertions...".format(len(insertion_candidates1) + len(insertion_candidates2)))
    partitions = form_partitions_ins(insertion_candidates1 + insertion_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(["INS", candidate[1], candidate[2], candidate[3], candidate[4], candidate[5], genotype])
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1][4] + cluster[1][1][4]
            genotype = "1/1"
            paired_candidates.append(["INS", candidate[1], candidate[2], candidate[3], reads, candidate[5], genotype])
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    # TANDEM DUPLICATIONS
    logging.info("Pairing {0} tandem duplications...".format(len(tandem_duplication_candidates1) + len(tandem_duplication_candidates2)))
    partitions = form_partitions(tandem_duplication_candidates1 + tandem_duplication_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(["DUP_TAN", candidate[1], candidate[2], candidate[3], candidate[4], candidate[5], genotype])
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1][5] + cluster[1][1][5]
            genotype = "1/1"
            paired_candidates.append(["DUP_TAN", candidate[1], candidate[2], candidate[3], candidate[4], reads, genotype])
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    # INTERSPERSED DUPLICATIONS
    logging.info("Pairing {0} interspersed duplications...".format(len(interspersed_duplication_candidates1) + len(interspersed_duplication_candidates2)))
    partitions = form_partitions_bnd(interspersed_duplication_candidates1 + interspersed_duplication_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes_breakends(partitions)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(["DUP_INT", candidate[1], candidate[2], candidate[3], candidate[4], candidate[5], candidate[6], genotype])
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1][6] + cluster[1][1][6]
            genotype = "1/1"
            paired_candidates.append(["DUP_INT", candidate[1], candidate[2], candidate[3], candidate[4], candidate[5], reads, genotype])
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    # BREAKENDS
    logging.info("Pairing {0} breakends...".format(len(breakend_candidates1) + len(breakend_candidates2)))
    partitions = form_partitions_bnd(breakend_candidates1 + breakend_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes_breakends(partitions)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(["BND", candidate[1], candidate[2], candidate[3], candidate[4], candidate[5], candidate[6], genotype])
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1][5] + cluster[1][1][5]
            genotype = "1/1"
            paired_candidates.append(["BND", candidate[1], candidate[2], candidate[3], candidate[4], reads, candidate[6], genotype])
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    return mergeAtrans1(paired_candidates, [], bam)

def mergeAtrans1(inter_can, intra_can, bam, mode=0):
    gen = '1/0' if mode == 0 else '0/1'
    sv_candidates = []
    for svtype in ['DEL', 'INS', 'INV', 'DUP_TAN', 'DUP_INT', 'BND']:
        if svtype == 'DEL':
            type_candidates = inter_can[svtype] + intra_can[svtype]
            type_candidates.sort(key=lambda x: (x[1], x[2], x[3]))
            for can in type_candidates:
                sv_candidates.append(CandidateDeletion(can[1], can[2], can[3], [can[4]], gen))
        elif svtype == 'INS':
            type_candidates = inter_can[svtype] + intra_can[svtype]
            type_candidates.sort(key=lambda x: (x[1], x[2], x[3]))
            for can in type_candidates:
                sv_candidates.append(CandidateInsertion(can[1], can[2], can[3], [can[4]], can[5], gen))
        elif svtype == 'BND':
            for can in inter_can[svtype]:
                sv_candidates.append(CandidateBreakend(can[1], can[2], can[3], can[4], can[6], [can[5]], bam, gen))
        elif svtype == 'INV':
            for can in inter_can[svtype]:
                sv_candidates.append(CandidateInversion(can[1], can[2], can[3], can[4], [can[5]], bam, gen))
        elif svtype == 'DUP_TAN':
            for can in inter_can[svtype] + intra_can[svtype]:
                sv_candidates.append(CandidateDuplicationTandem(can[0], can[1], can[2], can[3], can[4], [can[5]], bam, gen))
        elif svtype == 'DUP_INT':
            for can in inter_can[svtype]:
                sv_candidates.append(CandidateDuplicationInterspersed(can[2], can[3], can[4], can[5], can[7], [can[6]], bam, gen))
    return sv_candidates