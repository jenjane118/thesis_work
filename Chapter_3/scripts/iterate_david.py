#!/usr/bin python3
# coding=utf-8

#iterate_modules.py

# wrapper script to get gene lists and submit to David web service


taxon_id = 83332


# create string of gene ids for each module
gene_list = []

with open('module_list.txt', 'r') as f:
    mod_list = f.read().split('\n')
    print(mod_list)

        #delete unwanted output files
        #print("Deleting unwanted output files...")
        #run_fasta.delete_files()
        #append successful job to list
        #id_list.append(fasta_file)
        #write name file (mode='a' appends to file rather than overwrite)
        #with open('bovis_list.txt', mode='a', encoding='utf-8') as name_file:
        #    name_file.write('\n'.join(id_list) + '\n')

print("DONE!")
