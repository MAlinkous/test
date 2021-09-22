from readEvidenceFromMaxQuant import read_file
from readPeptideTable import readPeptideTable
from readTandemResults import tandem_matrix_analyse
from readMSGFResults import msgf_matrix_analyse


def exp_map(exp):
    exp_matrix = {}
    #fox = open("EOS_UPS1.txt", 'w')
    for peptide in exp:
        if peptide not in pep_table:
            continue
        info = pep_table[peptide]
        if len(info) == 9:
            pros = info[6]
            uni = info[7]
            ESS = info[2]
            EOS = info[4]
            N_OBS = info[3]
            N_EXP = info[5]
            PSS ="NA"

        elif len(info) == 5:
            pros = info[2]
            uni = info[3]
            PSS = info[1]
            N_EXP = "NA"
            ESS = "NA"
            EOS = "NA"
            N_OBS = "NA"
        if uni == '1':
            if len(pros) > 1:
                print("?????")
            if pros[0] not in exp_matrix:
                exp_matrix[pros[0]] = {}

            exp_matrix[pros[0]][peptide] = {"inten": exp[peptide], "ESS": ESS,
                                             "EOS": EOS, "N_OBS":N_OBS,
                                             "PSS": PSS, "uni": 1, "pros": pros, "N_EXP":N_EXP}
        elif uni == '0':
            for pro in pros:
                if pro not in exp_matrix:
                    exp_matrix[pro] = {}
                exp_matrix[pro][peptide] = {"inten": exp[peptide], "ESS": ESS,
                                             "EOS": EOS, "N_OBS": N_OBS,
                                             "PSS": PSS, "uni": 0, "pros": pros, "N_EXP":N_EXP}

    return exp_matrix

def map(pep_file, engine, pep_info, main_path, is_labeled):
    global pep_table
    print(".....Reading files from search engines......\n")
    if engine == 0:#MaxQaunt out
        peptide_matrix = read_file(pep_file)
    elif engine == 1:#X!tandem out
        peptide_matrix = tandem_matrix_analyse(pep_file, is_labeled, main_path)
    elif engine == 2:#MSGF+
        peptide_matrix = msgf_matrix_analyse(pep_file, is_labeled)
    else:
        print("Please input correct search engine code: 0 = MaxQuant, 1 = X!Tandem, 2 = MSGF+")
    print("File reading done\n")

    print(".....Reading peptide information......\n")
    pep_table = pep_info
    print("Peptide_table reading done\n")

    print(".....Mapping peptides to proteins......\n")
    map_matrix = {}
    ###########################################
    for e in peptide_matrix:
        map_matrix[e] = exp_map(peptide_matrix[e])
    ###############################################################
    print('Mapping finished\n')
    #print(map_matrix)


    return map_matrix



