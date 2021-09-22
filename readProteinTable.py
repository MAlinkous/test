#quantifition analysis
import time
import re
import gzip

def readProteinTable(main_path):
    #print(time.asctime( time.localtime(time.time()) ))
    file_name = main_path + "/docs/protein_table.gz"
    fi_protable = gzip.open(file_name, "rt")
    pro_dir = {}
    for line in fi_protable:
        if line.startswith("ID"):
            continue
        info = line.strip().split('\t')
        ID = info[0]
        Organism = info[1]
        Uniprot_Id = info[2]
        Date = time.strftime("%Y/%m/%d", time.localtime())
        Protein_Name = info[4]
        Gene_Symbol = info[5]
        Peptide_Count = int(info[6])
        Unique_Peptide_count = int(info[7])
        Observable_Unique_Peptide_Count = int(info[9])
        Gene_Groups = info[10]
        pro_family = info[11]
        if not re.search("\[", Gene_Groups):
            Gene_Groups = [Gene_Groups]
        else:
            Gene_Groups = eval(Gene_Groups)
        if not re.search("\[", pro_family):
            pro_family = [pro_family]
        else:
            pro_family = eval(pro_family)
        Description = info[12]
        Gene_Ontology_Biological_Process = info[13]
        Gene_Ontology_Cellular_Component = info[14]
        Gene_Ontology_Molecular_Function = info[15]
        pro_dir[Uniprot_Id] = {"ID":ID, "Organism":Organism, "Uniprot_Id":Uniprot_Id, "Date":Date, "Protein_Name":Protein_Name,
                   "Gene_Symbol":Gene_Symbol, "Peptide_Count":Peptide_Count, "Unique_Peptide_count":Unique_Peptide_count,
                   "Observable_Unique_Peptide_Count":Observable_Unique_Peptide_Count, "Gene_Groups":Gene_Groups,
                   "pro_family":pro_family, "Description":Description, "Gene_Ontology_Biological_Process":Gene_Ontology_Biological_Process,
                   "Gene_Ontology_Cellular_Component":Gene_Ontology_Cellular_Component, "Gene_Ontology_Molecular_Function":Gene_Ontology_Molecular_Function}
    fi_protable.close()
    return pro_dir

#readProteinTable("C:/Users/12891/OneDrive - sjtu.edu.cn/桌面/毕业设计/rebuild/芯片系统rebuild/EPtrans_v2/输出格式/protein_table.gz")
