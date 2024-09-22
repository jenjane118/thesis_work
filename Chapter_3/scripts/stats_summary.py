#!/usr/bin python3
# coding=utf-8

# stats_summary.py
# script to parse and print out summary stats from .tn_stats files generated by
# tpp pre-processing


#**********************************************************************************

def find_files(path):
    """Opens tn_stats text files found in specified directory and appends relevant info
    to dictionary. Uses os.scandir (https://docs.python.org/3/library/os.html).

    Input
    path                path to directory
    """

    import os
    import csv
    import re


    sample_dict = []


    csv_cols = ['read_count', 'density', 'max_count', 'NZ_mean']

    for entry in os.scandir(path):
        if entry.name.endswith((".tn_stats", ".py")):
            base = os.path.basename(entry)

            #path_list.append(entry.name)
            with open(entry, 'r') as file:
                data_dict = {}
                data_dict['sample'] = os.path.splitext(base)[0]
                for line in file:
                    m = re.search(r'(?<=# )\w+', line)
                    if m.group() in csv_cols:
                        #print(m.group())
                        n = re.search(r'(?<=: )\S+', line)
                        #print(n.group())
                        data_dict[m.group()] = n.group()
            sample_dict.append(data_dict)
            continue
        else:
            continue

    return sample_dict

#**********************************************************************************

def parse_tn_stats(sample_dict):
    """Parse files and make dictionary of attributes.
    """
    import csv

    csv_cols = ['sample', 'read_count', 'density', 'max_count', 'NZ_mean']

    with open ("combi_stats.csv", 'w', newline = '') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_cols)
        writer.writeheader()
        for data in sample_dict:
            writer.writerow(data)

    return csvfile


#**********************************************************************************

def parse_dat(path):
    """Parses tnseq_stats file from csv format
    """
    import os
    import csv

    data_dict = {}
    csv_cols = ['sample', 'read_count', 'density', 'max_count', 'NZ_mean']

    with open ("combi_stats.csv", 'a', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=csv_cols)
        for entry in os.scandir(path):
            if entry.name.endswith((".dat", ".py")):
                with open(entry) as csvfile:
                    data_dict = {}
                    reader = csv.reader(csvfile, delimiter='\t', )
                    next(reader, None)
                    for row in reader:
                        try:
                            sample      = row[0]
                            read_count  = row[6]
                            density     = row[1]
                            max_count   = row[5]
                            NZ_mean     = row[3]
                        except NameError as e:
                            print("Error", e)
                        data_dict = {'sample':sample, 'read_count':read_count, 'density':density, 'max_count':max_count,'NZ_mean':NZ_mean}
                        writer.writerow(data_dict)

    return outfile


#**********************************************************************************

def test():
    import csv
    files = find_files("/Users/jenniferstiens/git/tn_seq/lung_tnseq/tpp/")
    res = parse_tn_stats(files)
    #for x in res:
    #    print(x)

    parse = parse_dat("/Users/jenniferstiens/git/tn_seq/lung_tnseq/tpp/")
    #
    # with open("combi_stats.csv", 'r') as f:
    #     reader = csv.reader(f)
    #     for line in reader:
    #         print(line)
    #



########## main ###################################################################

if __name__ == "__main__":

   test()