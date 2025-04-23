import os, sys, glob

def main():

    infile=open('samlist','r')
    sams=[]
    labs=[]
    for line in infile:
        f=line.split()
        sams.append(f[1])
        labs.append(f[0])
    infile.close()

    for i in range(len(labs)):
        samfile=sams[i]
        print('python2 filter_dedup_sam.py {0}'.format(samfile))
        os.system('python2 filter_dedup_sam.py {0}'.format(samfile))
        i=i+1

main()
