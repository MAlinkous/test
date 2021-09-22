import re
import xml.dom.minidom
from xml.dom.minidom import parse
import re
import os
import operator
import numpy as np
import chardet


# read experiment design file
def experiment_design(experiment_design_file):
    exp = {}
    f = open(experiment_design_file, 'rb')
    data = f.read()
    file_encoding = chardet.detect(data).get('encoding')
    f.close()
    with open(experiment_design_file, encoding = file_encoding) as fi:
        for line in fi:
            if not line.startswith("File"):
                li = line.strip('\n').split('\t')
                filex = li[0]
                expx = li[1]
                fracx = li[2]
                masicx = li[3]
                if expx not in exp:
                    exp[expx] = [[filex, masicx]]
                else:
                    exp[expx].append([filex, masicx])
    return exp


# get the lowest score of FDR filter
def get_i(li):
    j = 0.0
    for i in range(1, len(li)):
        s2 = li[i].strip().split('\t')
        if s2[5] == '1':
            j = j + 1
        fdr = j / i
        if fdr >= 0.01:
            x = float(s2[3])
            return x
    return 0

#####################################################
#######################################################
###############################################################

# MASIC ANALYSE
def masic_analyse_TMT(filey):
    fi = open(filey)
    line1 = fi.readline().strip('\n').split('\t')
    matrix = {}
    exp_num = []
    exp = []
    for i, title in enumerate(line1):
        if re.search('Ion_[\d\.]+$', title):
            exp_num.append(i)
            exp.append(title)

    for line in fi:
        li = line.strip('\n').split('\t')
        scan_num = li[1]
        Ion_inten = li[exp_num[0]: exp_num[-1] + 1]
        matrix[scan_num] = dict(zip(exp, Ion_inten))

    return matrix, exp


def masic_analyse_label_free(filey):
    fi = open(filey)
    line1 = fi.readline().strip('\n').split('\t')
    matrix = {}
    for i, title in enumerate(line1):
        if title == "FragScanNumber":
            scan = i
        elif title == "PeakArea":
            inten = i
    for line in fi:
        lines = line.strip('\n').split('\t')
        matrix[lines[scan]] = float(lines[inten])
    return matrix


# get specified intensity of one peptide in specified channle
def merge(masic_matrix, ids, is_labeled, exp="null"):
    inten = 0
    if is_labeled == 1:
        for idx in ids:
            inten = float(masic_matrix[str(idx)][exp]) + inten

    else:
        #for idx in ids:
            #inten = float(masic_matrix[str(idx)]) + inten
        inten = np.mean([float(masic_matrix[str(idx)]) for idx in ids])
    return inten

#####################################################
#######################################################
###############################################################

def tandem_analyse(filex, main_path):
    tree = xml.dom.minidom.parse(filex)
    root = tree.documentElement
    source = root.getAttribute('label')
    p = re.compile(r':reversed$')
    fo = open(main_path + "/temp/" + 'analyse.txt', 'w')
    fo.write('group_id' + '\t' + 'expect' + '\t' + 'pep_sumI' +
             '\t' + 'hyperscore' + '\t' + 'seq' + '\t' + 'isReversed' + '\t' + 'scan' + '\n')

    for node in root.childNodes:
        if node.nodeType != 3 and node.getAttribute('type') == 'model':
            pep_sumI = node.getAttribute('sumI')
            group_id = node.getAttribute('id')
            group_expect = node.getAttribute('expect')
            isReversed = 0
            protein = node.getElementsByTagName("protein")
            groups = node.getElementsByTagName("group")
            for pro in protein:
                peptide = pro.getElementsByTagName('domain')
                hyperscore = peptide[0].getAttribute('hyperscore')
                expect = peptide[0].getAttribute('expect')
                seq = peptide[0].getAttribute('seq')
                note = pro.getElementsByTagName('note')
                text = note[0].firstChild.data
                if p.search(text) == None:
                    pass
                else:
                    isReversed = 1
                    break
            for gro in groups:
                if gro.getAttribute('label') == "fragment ion mass spectrum":
                    note2 = gro.getElementsByTagName('note')
                    scan_info = note2[0].firstChild.data
                    scan_num = re.search("scan=(\d+).*", scan_info).group(1)

            fo.write(group_id + '\t' + expect + '\t' + pep_sumI +
                     '\t' + hyperscore + '\t' + seq + '\t' + str(isReversed) + '\t' + scan_num + '\n')

    fo.close()

    fi2 = open(main_path + "/temp/" + 'analyse.txt', 'r')
    fo2 = open(main_path + "/temp/" + 'res.txt', 'w')
    lines1 = fi2.readlines()
    dic = {}
    dic2 = {}
    for i in range(1, len(lines1)):
        s1 = lines1[i].strip().split('\t')
        m = s1[0]
        dic2[i] = m
        dic[m] = float(s1[3])
    sorted_dic = sorted(dic.items(), key=operator.itemgetter(1), reverse=True)
    fo2.write(lines1[0])
    for i in sorted_dic:
        for j in dic2:
            if i[0] == dic2[j]:
                fo2.write(lines1[j])
                break
    fi2.close()
    fo2.close()
    filename_in = "FDR_output.txt"
    fi3 = open(main_path + "/temp/" + 'res.txt', 'r')
    #fo3 = open(main_path + "/temp/" + 'test_' + filename_in, 'w')
    lines2 = tuple(fi3.readlines())
    score = get_i(lines2)
    fi3.close()
    matrix = {}

    for i in lines1:
        s3 = i.strip().split('\t')
        if i.strip().split('\t')[0] != 'group_id' and s3[5] == '0' and float(s3[3]) > score:
            seq = s3[4]
            scan = s3[6]
            if seq not in matrix:
                matrix[seq] = {"scan": [scan]}
            else:
                matrix[seq]["scan"].append(scan)

    return matrix

#####################################################
#######################################################
###############################################################

def tandem_matrix_analyse(filein, is_labeled, main_path):
    path_temp = main_path + "/temp/"
    if not os.path.exists(path_temp):
        os.mkdir(path_temp)
    matrix = {}
    exp_des = experiment_design(filein)
    if is_labeled == 1:
        for exp_l in exp_des:
            channel = []
            for frac in exp_des[exp_l]:
                tandem_file = frac[0]
                masic_file = frac[1]
                tandem_output_trans = tandem_analyse(tandem_file, main_path)
                masic_output_trans_all = masic_analyse_TMT(masic_file)
                masic_output_trans = masic_output_trans_all[0]
                if len(channel) == 0:
                    channel = masic_output_trans_all[1]
                for expx in channel:
                    exp_name = expx + "_" + exp_l
                    if exp_name not in matrix:
                        matrix[exp_name] = {}

                for peptide in tandem_output_trans:
                    s_num = tuple(tandem_output_trans[peptide]["scan"])
                    for expy in channel:
                        if peptide not in matrix[expy + "_" + exp_l]:
                            matrix[expy + "_" + exp_l][peptide] = merge(masic_output_trans, s_num, 1, expy)
                        else:
                            matrix[expy + "_" + exp_l][peptide] += merge(masic_output_trans, s_num, 1, expy)

    else:
        for exp_l in exp_des:
            if exp_l not in matrix:
                matrix[exp_l] = {}
            for frac in exp_des[exp_l]:
                tandem_file = frac[0]
                masic_file = frac[1]
                tandem_output_trans = tandem_analyse(tandem_file, main_path)
                masic_output_trans = masic_analyse_label_free(masic_file)
                for peptide in tandem_output_trans:
                    s_num = tuple(tandem_output_trans[peptide]["scan"])
                    if peptide not in matrix[exp_l]:
                        matrix[exp_l][peptide] = merge(masic_output_trans, s_num, 0)
                    else:
                        matrix[exp_l][peptide] += merge(masic_output_trans, s_num, 0)

    return matrix
