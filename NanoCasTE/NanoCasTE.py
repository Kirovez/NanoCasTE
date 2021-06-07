from Bio import SeqIO
import pysam
import os
import re

def BLASt_vs_Target(reads, target_fasta):
    outFile = reads + "blast_target.tab"
    os.system('makeblastdb -in {0} -out {0} -dbtype nucl'.format(target_fasta))
    blastCMD = 'blastn -query {0} -db {1} -outfmt 6 -num_threads 150 -word_size 11 -out {2}'.format(reads, target_fasta, outFile)
    print(blastCMD)
    os.system(blastCMD)
    clipped_read_length = {seq.id:len(seq.seq) for seq in SeqIO.parse(reads, 'fasta')}
    return([outFile, clipped_read_length])

def _mergeBlastIntervals(blastTab, min_distance = 20):
    with open(blastTab) as inFile:
        hsps_per_read = {}

        for lines in inFile:
            sp = lines.split('\t')
            if (int(sp[7]) - int(sp[6])) > 100:
                if sp[0] not in hsps_per_read:
                    hsps_per_read[sp[0]] = []
                strand = "+"
                if int(sp[9]) < int(sp[8]):
                    strand = "-"

                hsps_per_read[sp[0]].append([int(sp[6]), int(sp[7]), [strand]])

        coverage_per_read = {}
        for reads in hsps_per_read:
            hsps_per_read[reads].sort(key = lambda x: x[0])
            if len(hsps_per_read[reads]) == 1:
                coverage_per_read[reads] = [hsps_per_read[reads][0][1] - hsps_per_read[reads][0][0], hsps_per_read[reads][0][-1][0]]
            else:
                new_splits = []
                for hsps in hsps_per_read[reads]:
                    if not new_splits:
                        new_splits = [hsps]
                    else:
                        if abs(new_splits[-1][1] - hsps[0]) <= min_distance:
                            #Extend the interval
                            new_splits[-1][1] = hsps[1] # change end
                            new_splits[-1][-1].append(hsps[-1][0])  # add strand
                            #print(new_splits)
                        else:
                            new_splits.append(hsps)

                ##count sum of intervals
                cov = 0
                strand = []
                #print(reads)
                for hsps in new_splits:
                    cov += hsps[1] - hsps[0]
                    strand += hsps[-1]
                    #print(strand)
                #print(','.join(set(strand)))
                coverage_per_read[reads] = [cov, set(strand)]
        return(coverage_per_read) #read:[coverage len, strand)

def filterClippedPartsByBLAST(split_clip_pos, split_clip_neg, outFileBLAST, minBlastHitLen=120, min_clipped_part_coverage = 0.6, long_clipped_report = 8000):
    # select blast positive reads
    blast_out, clipped_read_length = outFileBLAST
    ids_target_reads = {} # read:set(strand)
    cover_clipped = _mergeBlastIntervals(blast_out, min_distance=20)
    for reads in cover_clipped:
        cover_length = cover_clipped[reads][0]
        if cover_length > min_clipped_part_coverage * clipped_read_length[reads]:
            ids_target_reads[reads] = cover_clipped[reads][1]
            if clipped_read_length[reads] > long_clipped_report:
                print("!!!Long clipping!!!:", reads, 'Length of the clipped part:', clipped_read_length[reads], 'Blast covered: ', cover_clipped[reads][0], round(cover_length/clipped_read_length[reads], 1))

    cnt_clipped_reads_accepted = 0

    print('Number of selected reads:', len(ids_target_reads))

    filtered_split_clip_pos, filtered_split_clip_neg = {}, {}
    ## collect all indexes (split positions) for each chromosome
    for readIds in ids_target_reads:
        if 'positive' in readIds:
            chromosome, start, idx = readIds.split('positive')
            start, idx = int(start), int(idx)
            # print([idx], clipped_length_pos[chromosome], chromosome)
            if chromosome not in filtered_split_clip_pos:
                filtered_split_clip_pos[chromosome] = []
            filtered_split_clip_pos[chromosome].append(split_clip_pos[chromosome][idx])
            filtered_split_clip_pos[chromosome][-1].append('|'.join(ids_target_reads[readIds]))
        if 'negative' in readIds:
            chromosome, start, idx = readIds.split('negative')
            start, idx = int(start), int(idx)
            if chromosome not in filtered_split_clip_neg:
                filtered_split_clip_neg[chromosome] = []
            filtered_split_clip_neg[chromosome].append(split_clip_neg[chromosome][idx])
            filtered_split_clip_neg[chromosome][-1].append('|'.join(ids_target_reads[readIds]))

    return (filtered_split_clip_pos, filtered_split_clip_neg)


def _areOverlapped(win_p, win_n, min_distance=20):
    #     win_p: [203013, 203013, 1],
    #    win_n: [311704, 311704, 1]
    # they are not overlapped
    if (win_p[0] > win_n[1] and win_p[1] > win_n[1]) or (win_p[0] < win_n[0] and win_p[1] < win_n[0]):
        if (abs(win_p[1] - win_n[0]) <= min_distance) or (abs(win_p[0] - win_n[1]) <= min_distance):
            return True
        else:
            return False
    # they are overlapped
    return True


def mergeBothStrandWindow(merged_split_clip_pos, merged_split_clip_neg, min_distance=20):
    cnt_both_covered = 0
    overlapped = {}
    merged_intervals = {}  # chromosome: interval_start: [[pos start, pos end, number of reads, clipped pos],[neg start,neg end, number of reads, clipped neg]]
    for chrom in merged_split_clip_pos:
        merged_intervals[chrom] = {}
        overlapped[chrom] = {}
        if chrom in merged_split_clip_neg:
            for pos_i, winds_pos in enumerate(merged_split_clip_pos[chrom]):
                overlapped_indic = False
                for neg_i, winds_neg in enumerate(merged_split_clip_neg[chrom]):
                    if _areOverlapped(winds_pos, winds_neg, min_distance=min_distance):
                        winds_pos[-1], winds_neg[-1] = ','.join([str(i) for i in winds_pos[-1]]), \
                                                       ','.join([str(i) for i in winds_neg[-1]]),

                        # print(chrom, "start:", winds_pos[0], 'reads pos:', winds_pos[2], 'clipped pos:', winds_pos[-1],
                        #       'reads neg:', winds_neg[2], 'clipped neg:', winds_neg[-1])

                        merged_intervals[chrom][winds_pos[0]] = [winds_pos, winds_neg]
                        overlapped[chrom][winds_neg[0]] = 0
                        overlapped_indic = True
                        break
                if not overlapped_indic:
                    winds_pos[-1] = ','.join([str(i) for i in winds_pos[-1]])
                    merged_intervals[chrom][winds_pos[0]] = [winds_pos, [0, 0, 0, '']]

    for chrom in merged_split_clip_neg:
        for neg_i, winds_neg in enumerate(merged_split_clip_neg[chrom]):
            if chrom in overlapped:
                if winds_neg[0] not in overlapped[chrom]:
                    winds_neg[-1] = ','.join([str(i) for i in winds_neg[-1]])
                    merged_intervals[chrom][winds_neg[0]] = [[0, 0, 0, ''], winds_neg]
            else:
                winds_neg[-1] = ','.join([str(i) for i in winds_neg[-1]])
                merged_intervals[chrom] = {}
                merged_intervals[chrom][winds_neg[0]] = [[0, 0, 0, ''], winds_neg]

    return (merged_intervals)



## merge neigboring split positions and count
def mergeSplitCount(split_clip, min_distance=20):
    merged_split_clip = {}
    cnt = 0
    for chromosome in split_clip:
        #print(chromosome)
        merged_split_clip[chromosome] = []
        new_splits = []  # [start, end, count reads, [clip legnths list]]
        split_clip[chromosome].sort(key=lambda x: x[0])

        for i, spl_pos in enumerate(split_clip[chromosome]):
            #print(spl_pos)
            if not new_splits:
                new_splits.append([spl_pos[0], spl_pos[0], 1, ["{0}:{1}".format(spl_pos[1], spl_pos[2])]])
            else:
                if abs(new_splits[-1][1] - spl_pos[0]) <= min_distance:
                    # Extend the interval
                    new_splits[-1][2] += 1
                    new_splits[-1][1] = spl_pos[0]
                    new_splits[-1][-1].append("{0}:{1}".format(spl_pos[1], spl_pos[2]))  # add clip length
                else:
                    cnt += 1
                    new_splits.append([spl_pos[0], spl_pos[0], 1, ["{0}:{1}".format(spl_pos[1], spl_pos[2])]])

        merged_split_clip[chromosome] = new_splits
    return (merged_split_clip)

def defineClippedSizes(target_fasta, gRNAs):
    forward = []
    reverse = []
    for seq in SeqIO.parse(target_fasta, 'fasta'):
        t_seq = str(seq.seq).upper()
        rev_seq = str(seq.seq.reverse_complement()).upper()
        for grna in gRNAs:
            if grna.upper() in t_seq:
                matches = re.finditer(grna.upper(), t_seq)
                matches_positions = [len(t_seq) - match.start() for match in matches]
                forward += matches_positions
            elif grna.upper() in rev_seq:
                matches = re.finditer(grna.upper(), rev_seq)
                matches_positions = [len(rev_seq) - match.start() for match in matches]
                reverse += matches_positions
                # print(matches_positions)

            else:
                print(grna, 'Not found!!')
    return (forward + reverse)


def clippedLengthIsOK(clipped_length, expected_clipped_length, dev_from_expected_clipped_length):
    for exp in expected_clipped_length:
        if (min(exp, clipped_length) / max(exp, clipped_length)) >= dev_from_expected_clipped_length:
            return True
    return False


### read alignments and collect split coordinates

### read alignments and collect split coordinates
def getSplitPositions(bam_file, target_fasta, reads, grnas, dev_from_expected_clipped_length=0.7,
                      min_len_clipped=120, min_read_length=750, mapping_quality=10):
    expected_clipped_length = defineClippedSizes(target_fasta, grnas)
    print("Expected clipped length:", expected_clipped_length)
    clipped_fasta = "clipped_parts.fasta"
    outsplitted = open(clipped_fasta, 'w')
    cnt = 0
    cnt_i = 0
    seq_idx = SeqIO.index(reads, 'fastq')
    splitted_seqs = {}  # read_id:{chromosome:[split_start, split_end]}
    split_clip_neg = {} #chromosome:[[split position, clip length],[],]
    split_clip_pos = {} #chromosome:[[split position, clip length],[],]
    for algns in pysam.AlignmentFile(bam_file, 'rb'):
        if algns.mapping_quality >= mapping_quality and not algns.is_supplementary and algns.infer_read_length() >= min_read_length:
            splitted_seqs[algns.query_name] = [0, 0]
            cnt_i += 1
            # add chromosome to the dictionary
            if algns.reference_name not in split_clip_pos:
                split_clip_pos[algns.reference_name] = []
            if algns.reference_name not in split_clip_neg:
                split_clip_neg[algns.reference_name] = []
            if algns.is_reverse:
                if (algns.cigartuples[-1][0] == 5 or algns.cigartuples[-1][0] == 4) and algns.cigartuples[-1][
                    1] > min_len_clipped and \
                        clippedLengthIsOK(algns.cigartuples[-1][1], expected_clipped_length,
                                          dev_from_expected_clipped_length):
                    split_clip_neg[algns.reference_name].append([algns.reference_end, algns.cigartuples[-1][1]])
                    cnt += 1
                    # write fasta file of clipped part. The id of each sequence is position of splitted end on the reference +
                    # idx + index of this position in the list
                    clipped = str(seq_idx[algns.query_name].seq.reverse_complement())[-algns.cigartuples[-1][1]:]
                    outsplitted.write('>{0}\n{1}\n'.format(algns.reference_name + "negative" + \
                                                           str(algns.reference_end) + "negative" + str(
                        len(split_clip_neg[algns.reference_name]) - 1), clipped))

            else:
                if (algns.cigartuples[0][0] == 5 or algns.cigartuples[0][0] == 4) and algns.cigartuples[0][
                    1] > min_len_clipped and \
                        clippedLengthIsOK(algns.cigartuples[0][1], expected_clipped_length,
                                          dev_from_expected_clipped_length):
                    split_clip_pos[algns.reference_name].append([algns.reference_start, algns.cigartuples[0][1]])
                    cnt += 1
                    # write fasta file of clipped part. The id of each sequence is position of splitted end on the reference +
                    # idx + index of this position in the list
                    clipped = str(seq_idx[algns.query_name].seq)[0:algns.cigartuples[0][1]]
                    outsplitted.write('>{0}\n{1}\n'.format(
                        algns.reference_name + "positive" + str(algns.reference_start) + "positive" + \
                        str(len(split_clip_pos[algns.reference_name]) - 1), clipped))

    print("Number of clipped reads added", cnt)
    outsplitted.close()
    return ([split_clip_pos, split_clip_neg, clipped_fasta])


def getBed(sorted_bam, merged_intervals, outBed):
    # 4. select asymetrically covered intervals using knowledge about gRNA positions
    asymmetrical_cov = []
    symmetrical_cov = []
    i = 0

    with open(outBed + 'bed', 'w') as outBED, open(outBed, 'w') as outTab, pysam.AlignmentFile(sorted_bam, 'rb') as aln_bam:
        outTab.write("pp\tChromosome:Start..End\tNo + reads\tClipped lengths of + reads\tNo - reads\tClipped lengths of - reads\tTotal number of reads in the region\n")
        for chrom in merged_intervals:
            for insertions in merged_intervals[chrom]:
                i += 1
                pos, neg = merged_intervals[chrom][insertions]
                if pos[0] != 0:
                    int_start, int_end = pos[0] - 20, pos[1] + 20
                else:
                    int_start, int_end = neg[0] - 20, neg[1] + 20
                if int_start < 0:
                    int_start = 0
                #print("{0}\t{1}\t{2}\n".format(chrom, int_start, int_end))
                outBED.write("{0}\t{1}\t{2}\t{0}:{1}..{2}\n".format(chrom, int_start, int_end))
                #print([i for i in pysam.AlignmentFile(sorted_bam, 'rb').fetch(chrom, int_start, int_end)])
                total_reads_in_region = len([i for i in aln_bam.fetch(chrom, int_start, int_end)])
                toWrite =  "{0}:{1}..{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(chrom, int_start, int_end, pos[2], pos[3], neg[2], neg[3],
                                                                          total_reads_in_region)
                if pos[2] != 0 and neg[2] != 0:
                    print(toWrite)
                    symmetrical_cov.append(toWrite)
                else:
                    asymmetrical_cov.append(toWrite)

        cnt = 0
        for lines in symmetrical_cov:
            cnt += 1
            outTab.write(str(cnt) + "\t" + lines + "\n")
        print("Number of insertion sites covered by + and - reads:", cnt)

        for lines in asymmetrical_cov:
            cnt += 1
            outTab.write(str(cnt) + "\t" + lines + "\n")
        print("Total number of insertion sites detected:", cnt)

def main(sorted_bam, reads, target_fasta, guides, outBed, mapping_quality=60, min_clipped_part_coverage = 0.6,min_read_length=750, long_clipped_report = 8000):
    ## 1. get split positions on the reference and length of clipped parts for each strand
    split_clip_pos, split_clip_neg, clipped_fasta = getSplitPositions(sorted_bam, target_fasta, reads, guides, dev_from_expected_clipped_length=0.65,
                                                                      mapping_quality=mapping_quality,  min_read_length=min_read_length)
    # 2. blast clipped parts vs target and remove one without similarities
    blast_out = BLASt_vs_Target(clipped_fasta, target_fasta)
    filtered_split_clip_pos, filtered_split_clip_neg = filterClippedPartsByBLAST(split_clip_pos, split_clip_neg, blast_out, min_clipped_part_coverage = min_clipped_part_coverage, long_clipped_report = long_clipped_report)

    # 3. merge neighbor split positions amd their clipped sizes
    merged_split_clip_pos = mergeSplitCount(filtered_split_clip_pos, min_distance=20)
    merged_split_clip_neg = mergeSplitCount(filtered_split_clip_neg, min_distance=20)

    # 4. intesect pos and neg split positions (intervals after merging) to find out symmetrically and asymetrically covered covered intervals
    merged_intervals = mergeBothStrandWindow(merged_split_clip_pos, merged_split_clip_neg)

    # 5. generate bed file with insertion coordinates
    getBed(sorted_bam, merged_intervals, outBed)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('sorted_bam', help='path to bam file after minimap2 mapping')
    parser.add_argument('reads', help='path to fastq file of reads')
    parser.add_argument('target_fasta', help='path to target sequence in fasta format')
    parser.add_argument('guides', help='path to quide RNA fasta file')
    parser.add_argument('outBed', help='path to bed file')

    parser.add_argument('-mrl', '--min_read_length', help='minimum read length', default=750, type=int)
    parser.add_argument('-q', '--map_q', help='minimum mapping quality', default=60,type=int)
    parser.add_argument('-mlc', '--min_len_clipped', help='minimum coverage of clipped part', default=0.6,type=float)
    parser.add_argument('-mlcr', '--long_clipped_report', help='minimum length of the clipped part to report', default=5000,type=float)

    args = parser.parse_args()
    guides = [seq for seq in args.guides.split(",")]
    main(args.sorted_bam, args.reads, args.target_fasta, guides, args.outBed,
         mapping_quality=args.map_q, min_clipped_part_coverage=args.min_len_clipped,min_read_length=args.min_read_length, long_clipped_report = args.long_clipped_report)
