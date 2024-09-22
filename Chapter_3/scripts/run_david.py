#!/usr/bin python3
# coding=utf-8

#function to run access DAVID webservice (https://david.ncifcrf.gov/content.jsp?file=WS.html)

def run_david(fas_file, database):
    import os
    #print("~/git/mtb_modules/Scripts/fasta.py --sequence " + fasta_file + " --email j.j.stiens@gmail.com --program ssearch --stype dna --database " + database)
    os.system("~/git/mtb_modules/Scripts/fasta.py --sequence " + fas_file \
    + " --email j.j.stiens@gmail.com --program glsearch --stype dna --database " \
    + database + " --alignments 10 ")

def delete_files():
    # removes unneeded files from fasta.py output
    import os
    import glob
    patterns = ('*.pdf', '*.svg', '*.png', '*.jpg', '*.m9.txt', '*.m10.txt', '*.ids.txt', '*.error.txt', '*.accs.txt', '*.params', '*.xml')
    for p in patterns:
        for f in glob.glob(p):
            try:
                os.remove(f)
            except OSError as e:
                print("Error: %s : %s" % (f, e.strerror))


# *************************************************************************************
########## main ############

if __name__ == "__main__":

    file = "ncRv3461c.fasta"
    db = "233413-g"

    run_fasta(file, db)
    delete_files()
