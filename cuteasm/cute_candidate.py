class Candidate:
    """Candidate class for structural variant candidates. Candidates reflect the final SV types and can be merged from signatures of several reads.
    """
    def __init__(self, source_contig, source_start, source_end,genotype = "./."):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end

        self.type = None
        self.genotype = genotype


    # def get_key(self):
    #     return (self.source_contig, self.source_start, self.type)
    def get_key(self):
        return (self.source_contig, self.source_start, self.source_end)

    def set_genotype(self,genotype="0/1"):
        self.genotype=genotype

    def get_vcf_entry(self):
        raise NotImplementedError

class CandidateDeletion(Candidate):
    def __init__(self, source_contig, source_start, length, reads,  genotype = "./."):
       # assert source_end >= source_start, "Deletion end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(source_contig, source_end, source_start, reads)
        self.source_contig = source_contig
        #contig_length = bam.get_reference_length(source_contig)
        #0-based start of the deletion (first deleted base)
        self.source_start = int(max(0, source_start))
        #0-based end of the deletion (one past the last deleted base)
        self.source_end = int(source_start)+abs(int(length))
        self.length=abs(int(length))

        self.type = "DEL"
        self.reads = reads
        self.genotype = genotype


    def get_vcf_entry(self, sequence_alleles = True, reference = None, read_names = True):
        contig, start, end = self.get_key()
        filters = []
        if sequence_alleles:
            try:
                ref_allele = reference.fetch(contig, max(0, start-1), end).upper()
            except:
                ref_allele='Not_present'
            try:
                alt_allele = reference.fetch(contig, max(0, start-1), start).upper()
            except:
                alt_allele='Not_present'
        else:
            ref_allele = "N"
            alt_allele = "<" + self.type + ">"
        info_template="SVTYPE={0};END={1};SVLEN={2}"
        info_string = info_template.format(self.type, 
                                           end, 
                                           start - end)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=max(1, start),
                    id="PLACEHOLDERFORID",
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


class CandidateInversion(Candidate):
    def __init__(self, source_contig, source_start, source_end,types, reads, ref, genotype = "./."):
        assert source_end >= source_start, "Inversion end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(source_contig, source_end, source_start, reads)
        
       
        self.source_contig = source_contig
        contig_length = ref.get_reference_length(self.source_contig)
        #0-based start of the inversion (first inverted base)
        self.source_start = max(0, source_start)
        #0-based end of the inversion (one past the last inverted base)
        self.source_end = min(contig_length, source_end)
        self.layer=types
        self.length=source_end-source_start
        self.type = "INV"
        self.reads = reads
        self.genotype = genotype

        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


    def get_vcf_entry(self, sequence_alleles = True, reference = None, read_names = True):
        contig, start, end = self.get_key()
        if sequence_alleles:
            ref_allele='Not_present'
            alt_allele='Not_present'
            # try:
            #     ref_allele = reference.fetch(contig, start, end).upper()
            # except:
            #     ref_allele='Not_present'
            # try:
            #     alt_allele = "".join(self.complement.get(base.upper(), base.upper()) for base in reversed(ref_allele))
            # except:
            #     alt_allele='Not_present'
        else:
            ref_allele = "N"
            alt_allele = "<" + self.type + ">"
        info_template="SVTYPE={0};END={1};SVLEN={2};STRAND={3}"
        info_string = info_template.format(self.type, 
                                            end,self.length,self.layer)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=start+1,
                    id="PLACEHOLDERFORID",
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=".",
                    filter="PASS",
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


class CandidateInsertion(Candidate):
    def __init__(self, source_contig, source_start, length, reads, sequence,  genotype = "./."):
        self.source_contig=source_contig
        self.source_start=int(source_start)
        self.length=length
        self.source_end=int(source_start)

        self.type = "INS"
        self.reads = reads
        self.sequence = sequence
        self.genotype = genotype
    def get_source(self):
        return (self.source_contig, self.source_start, self.source_end)

    def get_key(self):
        return (self.source_contig, self.source_start, self.length)
    def get_vcf_entry(self, sequence_alleles = True, reference = None, read_names = True):
        contig, start, end = self.get_source()
        filters = []
        if sequence_alleles:
            try:
                ref_allele = reference.fetch(contig, max(0, start-1), start).upper()
                alt_allele = ref_allele + self.sequence
            except:
                ref_allele = "N"
                alt_allele = "<" + self.type + ">"
        else:
            ref_allele = "N"
            alt_allele = "<" + self.type + ">"
        info_template="SVTYPE={0};END={1};SVLEN={2}"
        info_string = info_template.format(self.type, 
                                           start, 
                                           self.length) 
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=max(1, start),
                    id="PLACEHOLDERFORID",
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


# class CandidateDuplication(Candidate):
#     def __init__(self, source_contig, source_start, source_end, copy_num, reads, ref, genotype = "./."):
#         assert source_end >= source_start, "Tandem duplication end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(source_contig, source_end, source_start, reads)
#         self.source_contig = source_contig
#         contig_length = ref.get_reference_length(source_contig)
#         #0-based start of the region (first copied base)
#         self.source_start = int(max(0, source_start))
#         #0-based end of the region (one past the last copied base)
#         self.source_end = int(min(contig_length, source_end))
#         self.type = "DUP"
#         self.reads = reads
#         self.genotype = genotype
#         self.copy=copy_num

#     def get_vcf_entry_as_ins(self, sequence_alleles = True, reference = None, read_names = True):
#         contig = self.source_contig
#         start = self.source_start
#         end = self.source_end
#         svtype = "INS"
        
#         ref_allele = "N"
#         alt_allele = "<" + svtype + ">"
#         info_template="SVTYPE={0};END={1};SVLEN={2}"
#         info_string = info_template.format(svtype, 
#                                            end, 
#                                            (end - start))
        
#         if read_names:
#             info_string += ";READS={0}".format(",".join(self.reads))
#         return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
#                     chrom=contig,
#                     pos=start+1,
#                     id="PLACEHOLDERFORID",
#                     ref=ref_allele,
#                     alt=alt_allele,
#                     qual=".",
#                     filter="PASS",
#                     info=info_string,
#                     format="GT",
#                     samples="{gt}".format(gt=self.genotype))


#     def get_vcf_entry_as_dup(self,reference=None, read_names = False):
#         contig = self.source_contig
#         start = self.source_start
#         end = self.source_end
#         length = self.source_end - self.source_start
#         svtype = "DUP"
#         info_template="SVTYPE={0};END={1};SVLEN={2}"
#         info_string = info_template.format(svtype, 
#                                            end, 
#                                            length)
#         if read_names:
#             info_string += ";READS={0}".format(",".join(self.reads))
#         return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
#                     chrom=contig,
#                     pos=start+1,
#                     id="PLACEHOLDERFORID",
#                     ref="N",
#                     alt="<" + svtype + ">",
#                     qual=".",
#                     filter="PASS" ,
#                     info=info_string,
#                     format="GT:CN",
#                     samples="{gt}:{cn}".format(gt=self.genotype,cn=self.copy))



class CandidateBreakend(Candidate):
    def __init__(self, source_contig, source_start,dest_contig, dest_start,types, reads, ref,length=0,forword=False,backword=False, genotype = "./."):
        self.source_contig = source_contig
        source_contig_length = ref.get_reference_length(source_contig)
        #0-based source of the translocation (first base before the translocation)
        self.source_start = min(source_contig_length, max(0, source_start))
        self.dest_contig = dest_contig
        try:
            dest_contig_length = ref.get_reference_length(dest_contig)
            #0-based destination of the translocation (first base after the translocation)
            self.dest_start = min(dest_contig_length, max(0, dest_start))
        except:
            self.dest_start=max(0, dest_start)
        self.layer=types
        self.length=length
        self.forword=forword
        self.backword=backword
        self.type = "BND"
        self.reads = reads
        self.genotype = genotype


    def get_source(self):
        return (self.source_contig, self.source_start)
    def get_key(self):
        return (self.source_contig, self.source_start, self.dest_start)

    def get_destination(self):
        return (self.dest_contig, self.dest_start)

    def get_vcf_entry(self,reference=None, read_names = True):
        source_contig, source_start = self.get_source()
        info_template="SVTYPE={0};Dest_chr={1};Dest_start={2};SVLEN={3}"
        info_string = info_template.format(self.type,self.dest_contig,self.dest_start,self.length)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=source_contig,
                    pos=source_start+1,
                    id="PLACEHOLDERFORID",
                    ref="N",
                    alt="<" + self.type + ">",
                    qual=".",
                    filter="PASS" ,
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))
class CandidateDuplicationTandem(Candidate):
    def __init__(self, source_contig, source_start, source_end, copies, fully_covered, reads, bam, genotype = "1/1"):
        assert source_end >= source_start, "Tandem duplication end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(source_contig, source_end, source_start, reads)
        self.source_contig = source_contig
        contig_length = bam.get_reference_length(source_contig)
        #0-based start of the region (first copied base)
        self.source_start = max(0, source_start)
        #0-based end of the region (one past the last copied base)
        self.source_end = min(contig_length, source_end)
        
        #number of additional copies
        self.copies = copies
        self.type = "DUP_TAN"
        self.reads = reads
        self.fully_covered = fully_covered
        self.genotype = genotype


    def get_destination(self):
        source_contig, source_start, source_end = self.get_key()
        return (source_contig, source_end, source_end + self.copies * (source_end - source_start))


    def get_vcf_entry_as_ins(self, sequence_alleles = False, reference = None, read_names = False):
        contig = self.source_contig
        start = self.source_start
        end = self.source_end
        svtype = "INS"
        filters = []
        if sequence_alleles:
            ref_allele = reference.fetch(contig, self.source_start, self.source_end).upper()
            alt_allele = ref_allele * (self.copies + 1)
        else:
            ref_allele = "N"
            alt_allele = "<" + svtype + ">"
        if not(self.fully_covered):
            filters.append("not_fully_covered")
        info_template="SVTYPE={0};END={1};SVLEN={2}"
        info_string = info_template.format(svtype, 
                                           end, 
                                           (end - start) * self.copies)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=start+1,
                    id="PLACEHOLDERFORID",
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


    def get_vcf_entry_as_dup(self, read_names = False):
        contig = self.source_contig
        start = self.source_start
        end = self.source_end
        length = self.source_end - self.source_start
        svtype = "DUP:TANDEM"
        filters = []
        if not(self.fully_covered):
            filters.append("not_fully_covered")
        info_template="SVTYPE={0};END={1};SVLEN={2}"
        info_string = info_template.format(svtype, 
                                           end, 
                                           length)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=start+1,
                    id="PLACEHOLDERFORID",
                    ref="N",
                    alt="<" + svtype + ">",
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT:CN",
                    samples="{gt}:{cn}".format(gt=self.genotype, cn=self.copies + 1))

# [-1,chrom,ins_start,ins_end,ins_flag,qs,len,readname,1,'pari']

class CandidateDuplicationInterspersed(Candidate):
    def __init__(self, source_contig, source_start, source_end, dest_contig, dest_start, dest_end, reads, bam, cutpaste=False, genotype = "1/1",reverse=False,repeats=1):
        # assert source_end >= source_start, "Interspersed duplication source end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(source_contig, source_end, source_start, reads)
        # assert dest_end >= dest_start, "Interspersed duplication destination end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(dest_contig, dest_end, dest_start, reads)
        self.source_contig = source_contig
        source_contig_length = bam.get_reference_length(source_contig)
        #0-based start of the region (first copied base)
        self.source_start = max(0, source_start)
        #0-based end of the region (one past the last copied base)
        self.source_end = min(source_contig_length, source_end)

        self.dest_contig = dest_contig
        dest_contig_length = bam.get_reference_length(dest_contig)
        #0-based start of the insertion (base after the insertion)
        self.dest_start = max(0, dest_start)
        #0-based end of the insertion (base after the insertion) + length of the insertion
        self.dest_end = min(dest_contig_length, dest_end)
        self.reverse=reverse
        self.cutpaste= cutpaste
        self.type = "DUP_INT"
        self.reads = reads
        self.genotype = genotype
        self.support=repeats


    def get_destination(self):
        return (self.dest_contig, self.dest_start, self.dest_end)

    def get_source(self):
        return (self.source_contig, self.source_start, self.source_end)
    def get_key(self):
        return (self.dest_contig, self.source_start, self.dest_start)


    def get_vcf_entry_as_ins(self, sequence_alleles = False, reference = None, read_names = False):
        contig, start, end = self.get_destination()
        svtype = "INS"
        filters = []
       
        ref_allele = "N"
        alt_allele = "<" + svtype + ">"
        info_template="SVTYPE={0};{1}END={2};SVLEN={3}"
        info_string = info_template.format(svtype, 
                                           "CUTPASTE;" if self.cutpaste else "", 
                                           start, 
                                           end - start)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=max(1, start),
                    id="PLACEHOLDERFORID",
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


    def get_vcf_entry_as_dup(self, read_names = True):
        contig, start, end = self.get_source()
        svtype = "DUP:INT"
        filters = []
        info_template="SVTYPE={0};{1}END={2};SVLEN={3};Dest_chr={4};Dest_start={5};Support={6};Reverse={7}"
        info_string = info_template.format(svtype, 
                                           "CUTPASTE;" if self.cutpaste else "", 
                                           end, 
                                           end - start,
                                           self.dest_contig,
                                            self.dest_start,
                                           self.support,"Ture" if self.reverse else "False")
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=start+1,
                    id="PLACEHOLDERFORID",
                    ref="N",
                    alt="<" + svtype + ">",
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


# class CandidateDuplicationInterspersed(Candidate):
#     def __init__(self, source_contig, source_start,dest_contig, dest_start,types, reads, ref, genotype = "./."):
#         self.source_contig = source_contig
#         source_contig_length = ref.get_reference_length(source_contig)
#         #0-based source of the translocation (first base before the translocation)
#         self.source_start = min(source_contig_length, max(0, source_start))
#         self.dest_contig = dest_contig
#         try:
#             dest_contig_length = ref.get_reference_length(dest_contig)
#             #0-based destination of the translocation (first base after the translocation)
#             self.dest_start = min(dest_contig_length, max(0, dest_start))
#         except:
#             self.dest_start=max(0, dest_start)
#         self.layer=types

#         self.type = "DUP_INT"
#         self.reads = reads
#         self.genotype = genotype


#     def get_source(self):
#         return (self.source_contig, self.source_start)
#     def get_key(self):
#         return (self.source_contig, self.source_start, self.dest_start)

#     def get_destination(self):
#         return (self.dest_contig, self.dest_start)

#     def get_vcf_entry(self,reference=None, read_names = True):
#         source_contig, source_start = self.get_source()
#         info_template="SVTYPE={0};Dest_chr={1};Dest_start={2}"
#         info_string = info_template.format(self.type,self.dest_contig,self.dest_start)
#         if read_names:
#             info_string += ";READS={0}".format(",".join(self.reads))
#         return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
#                     chrom=source_contig,
#                     pos=source_start+1,
#                     id="PLACEHOLDERFORID",
#                     ref="N",
#                     alt=self.layer,
#                     qual=".",
#                     filter="PASS" ,
#                     info=info_string,
#                     format="GT",
#                     samples="{gt}".format(gt=self.genotype))
