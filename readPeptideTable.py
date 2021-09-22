#quantifition analysis
import time
import re
import gzip

def readPeptideTable(main_path):
    file_name = main_path + "/docs/peptide_table.gz"
    #print(time.asctime( time.localtime(time.time()) ))
    fi_peptable = gzip.open(file_name, "rt")
    i = 0
    pep_dir = {}
    for line in fi_peptable:
        pep_info = eval(line.strip().split('\t')[0])
        pep_info.append(line.strip().split('\t')[1])
        if re.search('@', line):
            pep_dir[pep_info[1]] = pep_info
        elif re.search('\$', line):
            pep_dir[pep_info[0][1:]] = pep_info

    fi_peptable.close()
    return pep_dir

#readPeptideTable("C:/Users/12891/OneDrive - sjtu.edu.cn/桌面/毕业设计/rebuild/芯片系统rebuild/EPtrans_v2/输出格式/peptide_table.gz")