import sys


def repeat_region():
    methylation_file = open(sys.argv[1], 'r')
    overlap_file = open(sys.argv[2], 'r')
    methylation_dict = dict()
    overlap_dict = dict()
    for line in overlap_file:
        seq = line.strip().split('\t')
        if seq[0] not in overlap_dict:
            overlap_dict[seq[0]] = list()
        # overlap_dict[seq[0]].append((int(seq[1]), int(seq[2])))
        overlap_dict[seq[0]].append((int(seq[3])-1000, int(seq[3])+1000))
    # print(len(overlap_dict))
    # for chr in overlap_dict:
    #     print('%s\t%d'%(chr, len(overlap_dict[chr])))
    print('finish overlapping file')
    # print('chr10' in overlap_dict)
    cnt = 0
    for line in methylation_file:
        seq = line.strip().split('\t')
        if seq[0] not in methylation_dict:
            methylation_dict[seq[0]] = list()
        methylation_dict[seq[0]].append(int(seq[1]))
    print('finish methylation file')
    for chrom in methylation_dict:
        # print(chrom)
        if chrom in overlap_dict:
            sort_list = list()
            idx = 0
            for i in methylation_dict[chrom]:
                sort_list.append([i, 1, idx])
                idx += 1
            idx = 0
            for i in overlap_dict[chrom]:
                sort_list.append([i[0], 0, idx])
                sort_list.append([i[1], 2, idx])
                idx += 1
            sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
            overlap_set = set()
            for node in sort_list:
                if node[1] == 1: # set2(methylation)
                    if len(overlap_set) > 0:
                        cnt += 1
                elif node[1] == 0: # set1(repeat region) left
                    overlap_set.add(node[2])
                elif node[1] == 2: # set1(repeat region) right
                    overlap_set.remove(node[2])
            print('%s\t%d'%(chrom, cnt))
    print('total count=%d'%(cnt))

def gtf():
    methylation_file = open(sys.argv[1], 'r')
    overlap_file = open(sys.argv[2], 'r')
    methylation_dict = dict()
    overlap_dict = dict()
    idx = 1
    for line in overlap_file:
        if line[0] == '#':
            continue
        seq = line.strip().split('\t')
        if seq[0] not in overlap_dict:
            overlap_dict[seq[0]] = list()
        overlap_dict[seq[0]].append((int(seq[3]), int(seq[4]), idx))
        idx += 1
    gene_tot = idx - 1 # gene number
    # print(len(overlap_dict))
    # for chr in overlap_dict:
    #     print('%s\t%d'%(chr, len(overlap_dict[chr])))
    print('finish overlapping file')
    cnt = dict()
    for line in methylation_file:
        seq = line.strip().split('\t')
        if seq[0] not in methylation_dict:
            methylation_dict[seq[0]] = list()
        methylation_dict[seq[0]].append(int(seq[1]))
    print('finish methylation file')
    for chrom in methylation_dict:
        if chrom in overlap_dict:
            sort_list = list()
            idx = 0
            for i in methylation_dict[chrom]:
                sort_list.append([i, 1, idx])
                idx += 1
            idx = 0
            for i in overlap_dict[chrom]:
                sort_list.append([i[0], 0, idx])
                sort_list.append([i[1], 2, idx])
                idx += 1
            sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
            overlap_set = set()
            for node in sort_list:
                if node[1] == 1: # set2(methylation)
                    if len(overlap_set) > 0:
                        # print(node)
                        # print(methylation_dict[chrom][node[2]])
                        # print(overlap_set)
                        info_set = set()
                        for i in overlap_set:
                            info = overlap_dict[chrom][i][2]
                            # print(overlap_dict[chrom][i])
                            info_set.add(info)
                        # print(info_set)
                        for info in info_set:
                            if info not in cnt:
                                cnt[info] = 0
                            cnt[info] += 1
                        #     print(overlap_dict[chrom][i])
                        # exit(0)
                elif node[1] == 0: # set1(repeat region) left
                    overlap_set.add(node[2])
                elif node[1] == 2: # set1(repeat region) right
                    overlap_set.remove(node[2])
            print('%s\t%d'%(chrom, len(cnt)))
    # print(cnt)
    for i in range(1, gene_tot + 1, 1):
        if i not in cnt:
            cnt[i] = 0
    sorted_cnt = sorted(cnt.items(), key=lambda x:x[1], reverse=True)
    print(len(sorted_cnt))
    for i in range(len(sorted_cnt)):
        print(sorted_cnt[i])

def bpcnt():
    methylation_dict = dict()
    with open(sys.argv[1], 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            if seq[0] not in methylation_dict:
                methylation_dict[seq[0]] = list()
            methylation_dict[seq[0]].append(int(seq[1]))
    for chrom in methylation_dict:
        methylation_dict[chrom].sort()
    cnt = 0
    tot = 0
    for chrom in methylation_dict:
        for i in range(len(methylation_dict[chrom])):
            tot += 1
            if i == 0:
                if methylation_dict[chrom][i+1] - methylation_dict[chrom][i] <= 5:
                    cnt += 1
            elif i == len(methylation_dict[chrom]) - 1:
                if methylation_dict[chrom][i] - methylation_dict[chrom][i-1] <= 5:
                    cnt += 1
            else:
                if methylation_dict[chrom][i] - methylation_dict[chrom][i-1] <= 5 or methylation_dict[chrom][i+1] - methylation_dict[chrom][i] <= 5:
                    cnt += 1
        print(cnt)
    print(tot)

def dnaline():
    methylation_file = open(sys.argv[1], 'r')
    overlap_file = open(sys.argv[2], 'r')
    methylation_dict = dict()
    overlap_dict = dict()
    for line in overlap_file:
        if line[0] == '#':
            continue
        seq = line.strip().split('\t')
        if seq[5] not in overlap_dict:
            overlap_dict[seq[5]] = list()
        overlap_dict[seq[5]].append((int(seq[6]), int(seq[7]), seq[11]))
    print(len(overlap_dict))
    # for chr in overlap_dict:
    #     print('%s\t%d'%(chr, len(overlap_dict[chr])))
    print('finish overlapping file')
    cnt = dict()
    for line in methylation_file:
        seq = line.strip().split('\t')
        if seq[0] not in methylation_dict:
            methylation_dict[seq[0]] = list()
        methylation_dict[seq[0]].append(int(seq[1]))
    print('finish methylation file')
    for chrom in methylation_dict:
        if chrom in overlap_dict:
            sort_list = list()
            idx = 0
            for i in methylation_dict[chrom]:
                sort_list.append([i, 1, idx])
                idx += 1
            idx = 0
            for i in overlap_dict[chrom]:
                sort_list.append([i[0], 0, idx])
                sort_list.append([i[1], 2, idx])
                idx += 1
            sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
            overlap_set = set()
            for node in sort_list:
                if node[1] == 1: # set2(methylation)
                    if len(overlap_set) > 0:
                        # print(node)
                        # print(methylation_dict[chrom][node[2]])
                        # print(overlap_set)
                        info_set = set()
                        for i in overlap_set:
                            info = overlap_dict[chrom][i][2]
                            # print(overlap_dict[chrom][i])
                            info_set.add(info)
                        # print(info_set)
                        for info in info_set:
                            if info not in cnt:
                                cnt[info] = 0
                            cnt[info] += 1
                        #     print(overlap_dict[chrom][i])
                        # exit(0)
                elif node[1] == 0: # set1(repeat region) left
                    overlap_set.add(node[2])
                elif node[1] == 2: # set1(repeat region) right
                    overlap_set.remove(node[2])
            print('%s\t%d'%(chrom, len(cnt)))
    print(cnt)

def intron():
    methylation_file = open(sys.argv[1], 'r')
    overlap_file = open(sys.argv[2], 'r')
    methylation_dict = dict()
    overlap_dict = dict()
    remain_gene_id = ''
    remain_end = -1
    remain_start = -1
    for line in overlap_file:
        if line[0] == '#':
            continue
        seq = line.strip().split('\t')
        if seq[2] != 'exon':
            continue
        if seq[0] not in overlap_dict:
            overlap_dict[seq[0]] = list()
        start = int(seq[3])
        end = int(seq[4])
        geneid = seq[8].split('"')[3]
        if geneid == remain_gene_id:
            if start > remain_start:
                overlap_dict[seq[0]].append((remain_end, start))
            else:
                overlap_dict[seq[0]].append((end, remain_start))
        remain_end = end
        remain_start = start
        remain_gene_id = geneid
        if start == 29554:
            print(overlap_dict)
        #     exit(0)
    # print(len(overlap_dict))
    # for chr in overlap_dict:
    #     print(chr)
    #     print(overlap_dict[chr][:5])
    print('finish overlapping file')
    for line in methylation_file:
        seq = line.strip().split('\t')
        if seq[0] not in methylation_dict:
            methylation_dict[seq[0]] = list()
        methylation_dict[seq[0]].append(int(seq[1]))
    print('finish methylation file')
    cnt = 0
    for chrom in methylation_dict:
        # print(chrom)
        if chrom in overlap_dict:
            sort_list = list()
            idx = 0
            for i in methylation_dict[chrom]:
                sort_list.append([i, 1, idx])
                idx += 1
            idx = 0
            for i in overlap_dict[chrom]:
                sort_list.append([i[0], 0, idx])
                sort_list.append([i[1], 2, idx])
                idx += 1
            sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
            overlap_set = set()
            for node in sort_list:
                if node[1] == 1: # set2(methylation)
                    if len(overlap_set) > 0:
                        cnt += 1
                elif node[1] == 0: # set1(repeat region) left
                    overlap_set.add(node[2])
                elif node[1] == 2: # set1(repeat region) right
                    overlap_set.remove(node[2])
            print('%s\t%d'%(chrom, cnt))
    print('total count=%d'%(cnt))

def UTR_region():
    methylation_file = open(sys.argv[1], 'r')
    overlap_file = open(sys.argv[2], 'r')
    methylation_dict = dict()
    overlap_dict = dict()
    for line in overlap_file:
        seq = line.strip().split('\t')
        if seq[0] not in overlap_dict:
            overlap_dict[seq[0]] = list()
        overlap_dict[seq[0]].append((int(seq[3])-1000, int(seq[3])+1000))
    # print(len(overlap_dict))
    # for chr in overlap_dict:
    #     print('%s\t%d'%(chr, len(overlap_dict[chr])))
    print('finish overlapping file')
    # print('chr10' in overlap_dict)
    cnt = 0
    for line in methylation_file:
        seq = line.strip().split('\t')
        if seq[0] not in methylation_dict:
            methylation_dict[seq[0]] = list()
        methylation_dict[seq[0]].append(int(seq[1]))
    print('finish methylation file')
    for chrom in methylation_dict:
        # print(chrom)
        if chrom in overlap_dict:
            sort_list = list()
            idx = 0
            for i in methylation_dict[chrom]:
                sort_list.append([i, 1, idx])
                idx += 1
            idx = 0
            for i in overlap_dict[chrom]:
                sort_list.append([i[0], 0, idx])
                sort_list.append([i[1], 2, idx])
                idx += 1
            sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
            overlap_set = set()
            for node in sort_list:
                if node[1] == 1: # set2(methylation)
                    if len(overlap_set) > 0:
                        cnt += 1
                elif node[1] == 0: # set1(repeat region) left
                    overlap_set.add(node[2])
                elif node[1] == 2: # set1(repeat region) right
                    overlap_set.remove(node[2])
            print('%s\t%d'%(chrom, cnt))
    print('total count=%d'%(cnt))

def transcript_region():
    methylation_file = open(sys.argv[1], 'r')
    overlap_file = open(sys.argv[2], 'r')
    methylation_dict = dict()
    overlap_dict = dict()
    for line in overlap_file:
        seq = line.strip().split('\t')
        if seq[0] not in overlap_dict:
            overlap_dict[seq[0]] = list()
        overlap_dict[seq[0]].append((int(seq[3]), int(seq[4])))
    # print(len(overlap_dict))
    # for chr in overlap_dict:
    #     print('%s\t%d'%(chr, len(overlap_dict[chr])))
    print('finish overlapping file')
    # print('chr10' in overlap_dict)
    cnt = 0
    for line in methylation_file:
        seq = line.strip().split('\t')
        if seq[0] not in methylation_dict:
            methylation_dict[seq[0]] = list()
        methylation_dict[seq[0]].append(int(seq[1]))
    print('finish methylation file')
    for chrom in methylation_dict:
        # print(chrom)
        if chrom in overlap_dict:
            sort_list = list()
            idx = 0
            for i in methylation_dict[chrom]:
                sort_list.append([i, 1, idx])
                idx += 1
            idx = 0
            for i in overlap_dict[chrom]:
                sort_list.append([i[0], 0, idx])
                sort_list.append([i[1], 2, idx])
                idx += 1
            sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
            overlap_set = set()
            for node in sort_list:
                if node[1] == 1: # set2(methylation)
                    if len(overlap_set) > 0:
                        cnt += 1
                elif node[1] == 0: # set1(repeat region) left
                    overlap_set.add(node[2])
                elif node[1] == 2: # set1(repeat region) right
                    overlap_set.remove(node[2])
            print('%s\t%d'%(chrom, cnt))
    print('total count=%d'%(cnt))

# python overlap.py /home/user/liuyadong/zhongyu/hifi/ONT_hg002/genome_ch3_ont.bed /home/user/liuyadong/zhongyu/test/repeat-file repeat
# python overlap.py /home/user/liuyadong/zhongyu/hifi/ONT_hg002/genome_ch3_ont.bed /home/user/liuyadong/zhongyu/test/gencode.v41.annotation.gtf gtf
if __name__ == '__main__':
    if sys.argv[3] == 'repeat':
        repeat_region()
    if sys.argv[3] == 'gtf':
        gtf()
    if sys.argv[3] == '5bp':
        bpcnt()
    if sys.argv[3] == 'dnaline':
        dnaline()
    if sys.argv[3] == 'intron':
        intron()
    if sys.argv[3] == 'UTR':
        UTR_region()
    if sys.argv[3] == 'transcript':
        transcript_region()