import os, sys, glob

def main(argv=sys.argv):
    insamlist=argv[1:]
    for insam in insamlist:
        if insam[-13:]=='_filtered.sam':
            print('input file ends in filtered')
            samfile=insam
            bamfile=insam[:-3]+'bam'
            dedupbamfile=insam[:-4]+'_dedup'
        else:
            print('input file does not end in filtered')
            samfile=insam
            bamfile=insam[:-3]+'bam'
            dedupbamfile=insam[:-4]+'_dedup'
        files=glob.glob('q*')
        i=1
        while 'qd'+str(i) in files:
            i=i+1
        fn='qd'+str(i)
        ofile=open(fn,'w')
        ofile.write('#!/bin/bash\n')
        ofile.write('#SBATCH --mem=30GB\n')
        ofile.write('export TMPDIR=/home/mbeer3/tmp\n')
        ofile.write('samtools view -S -b {0} > {1} \n'.format(insam,bamfile))
        ofile.write('~/work/scripts/dnase_filter.sh {0} 10 1 {1}\n'.format(bamfile,dedupbamfile))
        ofile.close()
#        os.system('sbatch q'+str(i))
        os.system('sbatch -p troctolite qd'+str(i))

main()
