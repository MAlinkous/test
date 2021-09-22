#read MaxQuant Output File

import re

def read_lable_free(line1, fi):
    matrix = {}
    exp_site = {}
    lines1 = line1.strip('\n').split('\t')
    con = rev = inten = 0
    for i, value in enumerate(lines1):
        if re.search(r'Reverse', value):
            rev = i
        elif re.search(r'Potential contaminant', value):
            con = i
        elif re.search(r'^Intensity$', value):
            inten = i
    for i,exp in enumerate(lines1[inten+1: rev]):
        exp_name = exp[10:]
        matrix[exp_name] = {}
        exp_site[exp_name] = i + inten + 1
    #print(exp_site)
    for line in fi:
        li = line.strip('\n').split('\t')
        if line.startswith('Sequence') or li[con] == '+' or li[rev] == '+':
            continue
        seq = li[0]
        if set(li[inten+1: rev]) == {"0"}:
            continue
        for exp in exp_site:
            intensity = float(li[exp_site[exp]]) if li[exp_site[exp]] != "" else 0
            matrix[exp][seq] = intensity

    return matrix

def read_labled(line1, fi):
    matrix = {}
    exp_site = {}
    lines1 = line1.strip('\n').split('\t')
    con = rev = 0
    for i, value in enumerate(lines1):
        if re.search(r'Reverse', value):
            rev = i
        elif re.search(r'Contaminant', value):
            con = i
        elif re.search(r'Reporter intensity corrected \d+ .*', value):
            exp_tag = re.search(r'Reporter intensity corrected \d+ (.*)', value).group(1)
            channal = re.search(r'Reporter intensity corrected (\d+) .*', value).group(1)
            exp_name = exp_tag + "_" + channal
            exp_site[exp_name] = i
    for exp in exp_site:
        #print(exp)
        if exp not in matrix:
            matrix[exp] = {}
    #print(exp_site)
    exp_list = list(exp_site.values())
    for line in fi:
        li = line.strip('\n').split('\t')
        if line.startswith('Sequence') or li[con] == '+' or li[rev] == '+':
            continue
        seq = li[0]
        if set(li[exp_list[0]: exp_list[-1] + 1]) == {"0"}:
            continue
        for exp in exp_site:
            intensity = float(li[exp_site[exp]]) if li[exp_site[exp]] != "" else 0
            matrix[exp][seq] = intensity
    #print(matrix["tmt2b_10"])

    return matrix


def read_file(fi_name):
    # transfer file to intensity matrix
    #print('start reading evidence file...')
    #print('...')
    fi_trans = open(fi_name, 'r')
    matrix = {}

    line1 = fi_trans.readline()
    labled = 1 if re.search("Reporter intensity", line1) else 0

    if labled == 0:
        matrix = read_lable_free(line1, fi_trans)
    else:
        matrix = read_labled(line1, fi_trans)
    fi_trans.close()
    #print('file reading finished')
    return matrix

#read_file("../input/peptides_EXP.txt")
#read_file("C:/Users/12891/OneDrive - sjtu.edu.cn/桌面/毕业设计/rebuild/芯片系统rebuild/EPtrans_v2/输出格式/\MaxQuant_evidence&peptides/peptides_LFQ_paper.txt")