import os, sys, glob

def main():
    samps=glob.glob('Sample*.zip')
    print(len(samps), 'samples')
    i=1
    sfile=open('samlist','w')
    for s in samps:
        samp=s[:-4]
        r1files=glob.glob(samp+'/*_R1_001.fastq.gz')
        r2files=glob.glob(samp+'/*_R2_001.fastq.gz')
        r1files.sort()
        r2files.sort()
        files=glob.glob(samp+'/*.fastq.gz')
        print(samp,len(files),len(r1files),len(r2files))
        f=samp.split('_')
        lab=f[1]
        sfile.write(lab+'\tATAC_'+lab+'_hg38.sam\n')
        print(lab)
        if len(r1files)==1 and len(r2files)==1:
            ofile=open('qa'+str(i),'w')
            ofile.write('#!/bin/bash\n')
            ofile.write('#SBATCH --time=48:0:0\n')
            ofile.write('#SBATCH --nodes=1\n')
            ofile.write('CORE=ATAC_{0}_hg38\n'.format(lab))
            ofile.write('~/work/bowtie/bin/bowtie2 -x ~/work/bowtie/genomes/hg38 -1 {0} -2 {1} -S ${2}.sam\n'.format(','.join(r1files),','.join(r2files),'{CORE}'))
            ofile.write('sed -i \'/chrM/d;/random/d;/chrUn/d;/_/d;/EBV/d;/HCV/d;/HIV/d;/SV40/d;/HTLV/d;/CMV/d;/KSH/d;/HPV/d;/HBV/d;/MCV/d\' ${CORE}.sam\n')
            ofile.close()
            os.system('sbatch qa'+str(i))
            i=i+1
    sfile.close()

main()
