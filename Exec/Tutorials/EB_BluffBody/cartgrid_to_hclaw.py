import glob
from sys import argv
import os

def readandwrite_header(foldername):

    infile=open(foldername+"/Header",'r')

    outfile=open(foldername+"/Header_new","w")

    outfile.write("Hyperclaw-V1.1\n")
    line=infile.readline()
    for line in infile:
        outfile.write(line)

    outfile.close()
    infile.close()

    os.system("cp "+foldername+"/Header_new "+foldername+"/Header")
    os.system("rm "+foldername+"/Header_new")

#main
fpattern=argv[1]
folderlist = glob.glob(fpattern)
folderlist.sort()

print folderlist

for i in range(len(folderlist)):
    readandwrite_header(folderlist[i])
