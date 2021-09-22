from map import map
import numpy as np
import re
import time
from readProteinTable import readProteinTable
from header import header

def mean_alternative(input_seq):
    if len(input_seq) == 0:
        return "NA"
    else:
        return (sum(input_seq)/len(input_seq))

def MeanEOS(dictx):
    if type(dictx) == str:
        if dictx == "NA":
            result = 0
        else:
            result = float(dictx)
    else:
        result = np.mean([float(dictx[i]) for i in dictx])
    return result

def Test(value):
    if type(value) != dict:
        if value == "NA":
            result = -1
        else:
            result = float(value)
    else:
        result = np.mean([float(value[i]) for i in value])
    return result

def Quantify(peptides_usedx, quantify_method):
    peptides_used = {}
    if (not peptides_usedx) or sum([peptides_usedx[p]["inten"] for p in peptides_usedx]) == 0:#no characteristic peptides
        return 0
    for pep in peptides_usedx:
        if float(peptides_usedx[pep]["inten"]) != 0:
            peptides_used[pep] = peptides_usedx[pep]
    result = 0
    if quantify_method == 0:
        result = round(np.mean([peptides_used[x]["inten"] for x in peptides_used]), 0)
    elif quantify_method == 1:
        result = np.median([peptides_used[x]["inten"] for x in peptides_used])
    elif quantify_method == 2:
        result = sum([peptides_used[x]["inten"] for x in peptides_used])
    elif quantify_method == 3:
        all_peps = [peptides_used[x]["inten"] for x in peptides_used]
        result = round(np.mean(sorted(all_peps, reverse=True)[0:3]), 0)

    elif quantify_method == 4:
        all_eos = [peptides_used[x]["EOS"] for x in peptides_used]
        all_eos = [MeanEOS(x) for x in all_eos]
        all_inten = [peptides_used[x]["inten"] for x in peptides_used]
        eos_inten = dict(zip(all_eos, all_inten))
        eos_inten_sorted = sorted(eos_inten.items(), key=lambda x: x[0], reverse=True)
        result = round(np.mean([x[1] for x in eos_inten_sorted][0:3]), 0)

    elif quantify_method == 5:
        all_eos = [peptides_used[x]["EOS"] for x in peptides_used]
        all_eos = [MeanEOS(x) for x in all_eos]
        all_inten = [peptides_used[x]["inten"] for x in peptides_used]
        all_eos_inten = [(all_eos[i]*all_inten[i]) for i in range(len(all_eos))]
        eos_inten = dict(zip(all_eos_inten, all_inten))
        eos_inten_sorted = sorted(eos_inten.items(), key=lambda x: x[0], reverse=True)
        result = round(np.mean([x[1] for x in eos_inten_sorted][0:5]), 0)

    else:
        uni_num = 0
        for pep in peptides_used:
            if peptides_used[pep]["uni"] == 1:
                uni_num += 1
                result += peptides_used[pep]["inten"]
        if uni_num == 0:
            result = 0
        else:
            result = round(result/uni_num, 0)
    return result

def protein_matrix(matrix, method, quantify_method, test=1):
    pro_matrix = {}
    for exp in matrix:
        if exp not in pro_matrix:
            pro_matrix[exp] = {}
        for pro in matrix[exp]:
            if pro not in pro_matrix[exp]:
                pro_matrix[exp][pro] = {}
                pro_matrix[exp][pro]["peptides"] = matrix[exp][pro]
            Peptide_Count_DataSet = Unique_Peptide_Count_DataSet = Observable_Peptide_Count_DataSet = Observable_Unique_Peptide_Count_DataSet = 0
            inten = 0
            peptides_used = {}
            flag = False
            if method == 0:
                #flag = False
                for peptide in matrix[exp][pro]:
                    if matrix[exp][pro][peptide]["uni"] == 1 and matrix[exp][pro][peptide]["EOS"] != "NA":
                        flag = True
                        break
                if flag:
                    for peptide in matrix[exp][pro]:
                        if matrix[exp][pro][peptide]["EOS"] != "NA":
                            #inten += matrix[exp][pro][peptide]["inten"]
                            peptides_used[peptide] = matrix[exp][pro][peptide]
                    #for peptide in matrix[exp][pro]:
                        #peptides_used[peptide] = matrix[exp][pro][peptide]

            if method == 1:
                #flag = False
                for peptide in matrix[exp][pro]:
                    if matrix[exp][pro][peptide]["uni"] == 1 and matrix[exp][pro][peptide]["EOS"] != "NA":
                        if (MeanEOS(matrix[exp][pro][peptide]["EOS"]))>= 0.05:
                            flag = True
                            break
                if flag:
                   for peptide in matrix[exp][pro]:
                        if MeanEOS(matrix[exp][pro][peptide]["EOS"]) >= 0.05 and int(matrix[exp][pro][peptide]["N_OBS"]) >= 5 and int(matrix[exp][pro][peptide]["N_EXP"]) >= 4 and matrix[exp][pro][peptide]["ESS"] != "NA":
                            inten += matrix[exp][pro][peptide]["inten"]
                            peptides_used[peptide] = matrix[exp][pro][peptide]

                #else:
                    #del pro_matrix[exp][pro]

            if method == 2:
                flag = True
                for peptide in matrix[exp][pro]:
                    #inten += matrix[exp][pro][peptide]["inten"]
                    peptides_used[peptide] = matrix[exp][pro][peptide]


            if method == 3:
                #flag = False
                for peptide in matrix[exp][pro]:
                    if matrix[exp][pro][peptide]["uni"] == 1:
                        flag = True
                        break
                if flag:
                    for peptide in matrix[exp][pro]:
                        #if matrix[exp][pro][peptide]["EOS"] != "NA":
                            #inten += matrix[exp][pro][peptide]["inten"]
                        peptides_used[peptide] = matrix[exp][pro][peptide]
                #else:
                    #del pro_matrix[exp][pro]
##################################
            ###################################test############################################
###################################

##################################
            ###################################test############################################
###################################
            if not flag:#no characteristic peptides, delete pro
                del pro_matrix[exp][pro]
                continue
            inten = Quantify(peptides_used, quantify_method)
            for peptide in matrix[exp][pro]:
                Peptide_Count_DataSet += 1
                if matrix[exp][pro][peptide]["EOS"] != "NA":
                    Observable_Peptide_Count_DataSet += 1
                    if matrix[exp][pro][peptide]["uni"] == 1:
                        Observable_Unique_Peptide_Count_DataSet += 1
                if matrix[exp][pro][peptide]["uni"] == 1:
                    Unique_Peptide_Count_DataSet == 1
            pro_matrix[exp][pro]["intensity"] = inten
            pro_matrix[exp][pro]["Peptide_Count_DataSet"] = Peptide_Count_DataSet
            pro_matrix[exp][pro]["Unique_Peptide_Count_DataSet"] = Unique_Peptide_Count_DataSet
            pro_matrix[exp][pro]["Observable_Peptide_Count_DataSet"] = Observable_Peptide_Count_DataSet
            pro_matrix[exp][pro]["Observable_Unique_Peptide_Count_DataSet"] = Observable_Unique_Peptide_Count_DataSet

    return pro_matrix

def peptide_matrix(pro_matrix, pro_matrix_trans, method, exps):
    pep_matrix = {}
    for pro in pro_matrix_trans:
        if pro not in pep_matrix:
            pep_matrix[pro] = {}
    for exp in pro_matrix:
        for pro in pro_matrix[exp]:
            if pro_matrix[exp][pro]["intensity"] > 0:
                for pep in pro_matrix[exp][pro]["peptides"]:
                    if method == 0:
                        if pro_matrix[exp][pro]["peptides"][pep]["EOS"] != "NA":
                            if pep not in pep_matrix[pro]:
                                pep_matrix[pro][pep] = exps.copy()
                            pep_matrix[pro][pep][exp] = pro_matrix[exp][pro]["peptides"][pep]["inten"]
                    if method == 2:
                        if pep not in pep_matrix[pro]:
                            pep_matrix[pro][pep] = exps.copy()
                        pep_matrix[pro][pep][exp] = pro_matrix[exp][pro]["peptides"][pep]["inten"]
                    if method == 1:
                        if MeanEOS(pro_matrix[exp][pro]["peptides"][pep]["EOS"]) >= 0.05 and int(
                                pro_matrix[exp][pro]["peptides"][pep]["N_OBS"]) >= 5 and int(
                                pro_matrix[exp][pro]["peptides"][pep]["N_EXP"]) >= 4 and \
                                pro_matrix[exp][pro]["peptides"][pep]["ESS"] != "NA":
                            if pep not in pep_matrix[pro]:
                                pep_matrix[pro][pep] = exps.copy()
                            pep_matrix[pro][pep][exp] = pro_matrix[exp][pro]["peptides"][pep]["inten"]
    return pep_matrix

#def annotation_matrix():
def MakeExpLines(matrixX):
    Line = {}
    for exp in matrixX:
        if exp not in Line:
            Line[exp] = 0
    return Line

def pro_matrix_trans(matrix):
    exps = {}
    matrix_trans = {}
    for exp in matrix:
        if exp not in exps:
            exps[exp] = 0
    for exp in matrix:
        for pro in matrix[exp]:
            if pro not in matrix_trans:
                matrix_trans[pro] = exps.copy()
            matrix_trans[pro][exp] = matrix[exp][pro]["intensity"]
#################
    ##################test
#############
    matrix_trans_copy = matrix_trans.copy()
    for pro in matrix_trans_copy:
        #print(matrix_trans_copy[pro])
        if sum([matrix_trans_copy[pro][exp] for exp in matrix_trans_copy[pro]])== 0:
            del matrix_trans[pro]
    return matrix_trans, exps
#################
    ##################test
#############

def peptide_dir(pep_info, peptide):
    pep_info_single = pep_info[peptide]
    if len(pep_info_single) == 5:
        pros = pep_info_single[2]
        uni = pep_info_single[3]
        PSS = pep_info_single[1]
        ESS = "NA"
        EOS = "NA"
        N_OBS = "NA"
        N_EXP = "NA"
        idx = pep_info_single[-1]
        peptideAtlas_id = "NA"
    else:
        pros = pep_info_single[6]
        uni = pep_info_single[7]
        ESS = pep_info_single[2]
        EOS = pep_info_single[4]
        N_OBS = pep_info_single[3]
        N_EXP = pep_info_single[5]
        PSS = "NA"
        idx = pep_info_single[-1]
        peptideAtlas_id = pep_info_single[0]
    pep_dir = {"ESS": ESS,"EOS": EOS, "EOS": EOS, "N_OBS":N_OBS,
               "PSS": PSS, "uni": uni, "pros": pros, "N_experiments":N_EXP, "ID":idx, "peptideAtlas_id":peptideAtlas_id}
    return pep_dir

def write_list_dict(lidi):
    if type(lidi) == list:
        result = ""
        for i in lidi:
            result = result + str(i) + "; "
        result = result[:-2]
    elif type(lidi) == dict:
        result = ""
        for i in lidi:
            result = result + str(i) + ": " + str(lidi[i]) + "; "
        result = result[:-2]
    else:
        result = str(lidi)
    return result


def write_table(pro_info, pep_info, pep_file, engine, method, quantify_method, is_labeled, main_path, path = "./", tag = ""):
    date = time.strftime("%Y-%m-%d", time.localtime())
    protein_table = pro_info
    map_matrix = map(pep_file, engine, pep_info, main_path, is_labeled)
    print(".....Writing tables......\n")
    pro_inten_matrix = protein_matrix(map_matrix, method, quantify_method)
    #print(pro_inten_matrix)
    ##############################protein matrix writing
    fo_pro_matrix = open(path + "protein_matrix_" + tag + ".txt", "w")
    fo_pro_matrix.write("Protein_ID\tUniprot_Identifier")
    pro_matrix_to_write, experiments = pro_matrix_trans(pro_inten_matrix)
    #print(pro_matrix_to_write)
    for exp in experiments:
        fo_pro_matrix.write("\t" + exp)
    fo_pro_matrix.write("\n")
    for pro in pro_matrix_to_write:
        fo_pro_matrix.write(protein_table[pro]["ID"] + "\t" + pro)
        for exp in pro_matrix_to_write[pro]:
            fo_pro_matrix.write('\t' + str(pro_matrix_to_write[pro][exp]))
        fo_pro_matrix.write("\n")
    fo_pro_matrix.close()
    ####################################
    pep_matrix = peptide_matrix(pro_inten_matrix, pro_matrix_to_write, method, experiments)
    ##############################################
    fo_pep_matrix = open(path + "peptide_matrix_" + tag + ".txt", "w")
    fo_pep_matrix.write("Peptide_ID\tSequence\tProtein_Uniprot_Identifier(s)")
    for exp in experiments:
        fo_pep_matrix.write("\t" + exp)
    fo_pep_matrix.write("\tIs_Unique\tEOS\n")
    for pro in pep_matrix:
        fo_pep_matrix.write("\n")
        for pep in pep_matrix[pro]:
            pep_dir = peptide_dir(pep_info, pep)
            fo_pep_matrix.write(pep_dir["ID"] + "\t" + pep + "\t" + write_list_dict(pep_dir["pros"]))
            for exp in pep_matrix[pro][pep]:
                fo_pep_matrix.write("\t" + str(pep_matrix[pro][pep][exp]))
            fo_pep_matrix.write('\t' + str(pep_dir["uni"]) + '\t' + write_list_dict(pep_dir["EOS"]) + "\n")
    fo_pep_matrix.close()
    ##################################################
    fo_annatation = open(path + "annotation_" + tag + ".txt", "w")
    fo_annatation.write(header())
    for pro in pep_matrix:
        Peptide_Count_DataSet = Unique_Peptide_Count_DataSet = 0
        for pep in pep_matrix[pro]:
            pep_dir = peptide_dir(pep_info, pep)
            Peptide_Count_DataSet += 1
            if str(pep_dir["uni"]) == "1":
                Unique_Peptide_Count_DataSet += 1
        pro_annotation = pro_info[pro]
        fo_annatation.write("\n<Protein start>\n#")
        fo_annatation.write(pro_annotation["ID"] + "\tHomo Sapiens\t" + pro + "\t" + date + "\t" +
        pro_annotation["Protein_Name"] + "\t" + pro_annotation["Gene_Symbol"] + "\t" +
        write_list_dict(pro_annotation["Gene_Groups"]) + '\t' + write_list_dict(pro_annotation["pro_family"]) + "\t" + str(Peptide_Count_DataSet) + "\t" +
        str(Unique_Peptide_Count_DataSet) + "\t" + str(pro_annotation["Peptide_Count"]) + "\t" + str(pro_annotation["Unique_Peptide_count"]) + "\t" +
        pro_annotation["Description"] + "\t" + pro_annotation["Gene_Ontology_Biological_Process"] + pro_annotation["Gene_Ontology_Cellular_Component"] +
        pro_annotation["Gene_Ontology_Molecular_Function"] + "\n")
        fo_annatation.write("\t\t<peptide start>\n\t\t")
        for pep in pep_matrix[pro]:
            fo_annatation.write("!")
            pep_dir = peptide_dir(pep_info, pep)
            Is_Predict = "0" if pep_dir["PSS"] == "NA" else "1"
            Proteins_id = []
            for protein in eval(str(pep_dir["pros"])):
                Proteins_id.append(pro_info[protein]["ID"])
            fo_annatation.write(pep_dir["ID"] + "\t" + pep_dir["peptideAtlas_id"] + "\t" + date + "\t" + write_list_dict(pep_dir["pros"]) +
            "\t" + write_list_dict(Proteins_id) + "\t" + pep_dir["uni"] + "\t" + pep + "\t" + write_list_dict(pep_dir["ESS"]) + "\t" + write_list_dict(pep_dir["N_OBS"]) +
            "\t" + write_list_dict(pep_dir["EOS"]) + "\t" + pep_dir["N_experiments"] + "\t" + write_list_dict(pep_dir["PSS"]) + "\t" + Is_Predict + "\n\t\t")
        fo_annatation.write("<peptide end>\n<protein End>\n\n")
    fo_annatation.close()

    print(".....writing tables done!......\n")




