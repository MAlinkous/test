from write_tables import write_table
from readProteinTable import readProteinTable
from readPeptideTable import readPeptideTable
import re
import time
import sys
import getopt
import os

def usage():
    print("==========")
    print("Usage:")
    print("-i input file")
    print("-p output filePath")
    print("-t output file tag")
    print("-h help message")
    #print("-m quantificaton method")
    print("-e search engine")
    print("-l is labeled")
    #print("-w characteristic peptide selection criterion")
    print("==========")

def getargv():
    argx = {}
    #argx["main_path"] = re.search("(.+)/[a-zA-Z\.]+$", os.path.abspath('.')).group(1)
    #print(os.path.abspath('.'))
    for root, dirs, files in os.walk(os.path.abspath('.')):
        if "localization.localization" in files:
            localization = root
    argx["main_path"] = localization
    #print(localization)

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:p:t:m:e:w:l:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':

            #print('python ePoTrans.py -i <inputfile> -p <outputfilepath> -t <output file tag> -m <quantifi\
                  #caton method> -e <search engine> -w <characteristic peptide selection criterion> -l <is labeled>')

            print('python ePoTrans.py -i <inputfile> -p <outputfilepath> -t <output file tag> -e <search engine> -l <is labeled>')
            usage()
            sys.exit()
        elif opt in '-i':
            argx['fi'] = arg
        elif opt in '-p':
            argx['fopath'] = arg + '/'
        elif opt in '-t':
            argx['tag'] = arg
        elif opt in '-l':
            argx['is_labeled'] = arg
        elif opt in '-m':
            argx['quantification_method'] = int(arg)
        elif opt in '-e':
            argx['engine'] = int(arg)
        elif opt in '-w':
            argx['selection_method'] = int(arg)
    if 'fopath' not in argx:
        path = argx["main_path"] + "/output/"
        if not os.path.exists(path):
            os.mkdir(path)
        argx['fopath'] = path
    if 'quantification_method' not in argx:
        argx['quantification_method'] = 0
    if 'tag' not in argx:
        argx['tag'] = ""
    if 'selection_method' not in argx:
        argx['selection_method'] = 0
    if 'fi' not in argx:
        print("input file missed! use '-i input_file'")
    if 'fopath' in argx:
        if not os.path.exists(argx['fopath']):
            print('output path not exit!')
            sys.exit(2)
    if argx['is_labeled'] == "1":
        argx['quantification_method'] = 0
    else:
        argx['quantification_method'] = 3
    return argx

def main():
    args = getargv()
    pep_file = args["fi"]
    engine = args["engine"]
    selection_method = args["selection_method"]
    quantification_method = args["quantification_method"]
    path = args["fopath"]
    tag = args["tag"]
    is_labeled = int(args['is_labeled'])
    #test = int(args['test'])
    main_path = args["main_path"]
    pro_info = readProteinTable(main_path)
    pep_info = readPeptideTable(main_path)
    #print(args)

    write_table(pro_info, pep_info, pep_file, engine, selection_method, quantification_method, is_labeled, main_path, path, tag)

if __name__ == "__main__":
    main()

