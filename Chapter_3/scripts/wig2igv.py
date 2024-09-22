#
#               Creative Commons License
#           Attribution-NonCommercial 3.0 Unported
#        Copyright 2014. Michael A. DeJesus & Thomar R. Ioerger.
#

import sys

def hash_genes(path):
    hash = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        if tmp[2] != "gene": continue
        start, end = int(tmp[3]), int(tmp[4])
        features = dict([tuple(f.split("=")) for f in tmp[8].strip().split(";")])
        orf = features.get("ID", "missing")
        if orf == "missing": orf = features.get("Name", "missing")
        for i in range(start,end+1):
            if i not in hash:
                hash[i] = orf
    return hash




if "-f" in sys.argv and  "-gff" in sys.argv:
    path = sys.argv[sys.argv.index("-f")+1]
    gff_path = sys.argv[sys.argv.index("-gff")+1]
else:
  print "usage: python wig2igv.py -f <.wig file> -gff <.gff3 annotation file>"
  sys.exit(0)



# #type=COPY_NUMBER
# Chromosome	Start	End	Feature	reads	TAs
# NC_000962.2	59	61	Rv0001	0	1

print "#type=COPY_NUMBER"
print "Chromosome	Start	End	Feature	reads	TAs"

# note: this will catch overlaps


hash = hash_genes(str(gff_path))

for line in open(path):
    if line.startswith("#"): continue
    if line.startswith("variableStep"):
        chrom = line.strip().split("=")[1]
        continue

    pos, read = int(line.split()[0]), int(line.split()[1])
    orf = hash.get(pos, "non-coding")
    print "%s\t%s\t%s\t%s\t%s\t%s" % (chrom,pos,pos+2,orf,read,1)
