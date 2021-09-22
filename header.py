#header of annotation file

def header():
    h = (">DATA INFORMATION BEGIN\n\
%Database =\n\
%Data type =\n\
%Data id =\n\
%Data tittle =\n\
%Submission date =\n\
%Last modified date =\n\
%Data contributor =\n\
%Data contact name =\n\
%Data contact email =\n\
%Data contact phone =\n\
%Data contact institute =\n\
%Data contact address =\n\
%Data contact city =\n\
%Data contact country =\n\
%Data describtion =\n\
%Data URL =\n\
%Source raw data type =\n\
%Source raw data URL =\n\
%Source data publication =\n\
%Source data contributor =\n\
%Source data contact name =\n\
%Source data contact email =\n\
%Search enegine =\n\
%Organism =\n\
%Organism part =\n\
%Diseases =\n\
%Instrument =\n\
%Quantification Method=\n\
>DATA INFORMATION END\n\n" + ">ANNOTION BEGIN\n" +
"#ID\tOrganism\tUniprot_Id\tDate\tProtein_Name\tGene_Symbol\tGene_Groups\tProtein_Families\t" +
"Peptide_Count_DataSet\tUnique_Peptide_Count_DataSet\t" +
"Peptide_Count_Theoretically\tUnique_Peptide_Count_Theoretically\t" +
"Description\tGene_Ontology_Biological_Process\tGene_Ontology_Cellular_Component\tGene_Ontology_Molecular_Function\n\t\t" +
"!ID\tPeptideAtlas_id\tDate	Proteins\tProteins_Id\tIs_Unique\tSequence\tESS\tOBS\tEOS\tN_Expriments\tPSS\tIs_predict"
    )

    return h


