import os, sys, glob

def main():
    samps='ESC DE-12h DE-24h DE-36h DE-48h DE-60h DE-72h'.split()
    i=1
    efile=open('exptlist','w')
    for s in samps:
        bams=glob.glob('ATAC_{0}-?_hg38_dedup.bam'.format(s))
        core='ATAC_{0}'.format(s)
        print(s,len(bams),bams)
        sbam=' '.join(bams)
        ofile=open('qb'+str(i),'w')
        ofile.write('#!/bin/bash\n')
        ofile.write('#SBATCH --time=48:0:0\n')
        ofile.write('#SBATCH --nodes=1\n')
        ofile.write('CORE={0}_hg38\n'.format(core))
        ofile.write('macs2 callpeak -t '+sbam+' --nomodel -g hs -n ${CORE} -B --tempdir /home/mbeer3/tmp\n')
        ofile.write('sort -k1,1 -k2,2n -T /home/mbeer3/tmp ${CORE}_treat_pileup.bdg > ${CORE}_sorted.bdg\n')
        ofile.write('bedGraphToBigWig ${CORE}_sorted.bdg ~/work/scripts/hg38.chrom.sizes ${CORE}.bigwig\n')
        ofile.close()
        os.system('sbatch qb'+str(i))
        efile.write('{0}_hg38\n'.format(core))
        i=i+1
    efile.close()

main()
