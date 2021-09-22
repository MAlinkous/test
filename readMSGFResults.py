import re
import numpy as np
import chardet

# read experiment design file
def experiment_design(experimrnt_design_file):
    exp = {}
    f = open(experimrnt_design_file, 'rb')
    data = f.read()
    file_encoding = chardet.detect(data).get('encoding')
    f.close()
    with open(experimrnt_design_file, encoding=file_encoding) as fi:
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


# MSGF output ANALYSE
def msgf_analyse(filex):
    fi = open(filex, "r")
    matrix = {}
    for line in fi:
        if line.startswith("#SpecFile"):
            line1s = line.strip("\n").split("\t")
            for i, value in enumerate(line1s):
                if value == "ScanNum":
                    scan_index = i
                elif value == "Peptide":
                    peptide_index = i
                elif value == "PepQValue":
                    fdr_index = i
            continue
        if re.search("REV", line):
            continue
        li = line.strip().split("\t")
        scan = int(li[scan_index])
        peptide_raw = li[peptide_index]
        fdr = float(li[fdr_index])
        peptide = ""
        for i in peptide_raw:
            if i.isalpha():
                peptide += i
        if fdr <= 0.01:
            if peptide not in matrix:
                matrix[peptide] = {"scan": [scan]}
            else:
                matrix[peptide]["scan"].append(scan)
    fi.close()
    return matrix


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
        #elif title == "StatMomentsArea":
            inten = i
    for line in fi:
        lines = line.strip('\n').split('\t')
        matrix[lines[scan]] = float(lines[inten])
    return matrix

# get specified intensity of one peptide in specified channle
def merge(masic_matrix, ids, is_labeled, exp = "null"):
    inten = 0
    if is_labeled == 1:
        for idx in ids:
            inten = float(masic_matrix[str(idx)][exp]) + inten
    else:
        #for idx in ids:
            #inten = float(masic_matrix[str(idx)]) + inten
        inten = np.mean([float(masic_matrix[str(idx)]) for idx in ids])
    return inten


# matrix analyse
def msgf_matrix_analyse(filein, is_labeled):
    matrix = {}
    exp_des = experiment_design(filein)

    if int(is_labeled) == 1:
        for exp_l in exp_des:
            channel = []
            for frac in exp_des[exp_l]:
                msgf_file = frac[0]
                masic_file = frac[1]
                msgf_output_trans = msgf_analyse(msgf_file)
                masic_output_trans_all = masic_analyse_TMT(masic_file)
                masic_output_trans = masic_output_trans_all[0]
                if len(channel) == 0:
                    channel = masic_output_trans_all[1]
                for expx in channel:
                    exp_name = expx + "_" + exp_l
                    #print(exp_name)
                    if exp_name not in matrix:
                        matrix[exp_name] = {}

                for peptide in msgf_output_trans:
                    s_num = tuple(msgf_output_trans[peptide]["scan"])
                    for expy in channel:
                        if peptide not in matrix[expy + "_" + exp_l]:
                            matrix[expy + "_" + exp_l][peptide] = merge(masic_output_trans, s_num, 1, expy)
                        else:
                            matrix[expy + "_" + exp_l][peptide] += merge(masic_output_trans, s_num, 1, expy)

    else:
        for exp_l in exp_des:
            #print(exp_l)
            if exp_l not in matrix:
                matrix[exp_l] = {}
            for frac in exp_des[exp_l]:
                msgf_file = frac[0]
                masic_file = frac[1]
                msgf_output_trans = msgf_analyse(msgf_file)
                #print(msgf_output_trans)
                masic_output_trans = masic_analyse_label_free(masic_file)
                #print(masic_output_trans)
                for peptide in msgf_output_trans:
                    s_num = tuple(msgf_output_trans[peptide]["scan"])
                    if peptide not in matrix[exp_l]:
                        matrix[exp_l][peptide] = merge(masic_output_trans, s_num, 0)
                    else:
                        matrix[exp_l][peptide] += merge(masic_output_trans, s_num, 0)
    return matrix
