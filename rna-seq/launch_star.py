import os, sys, os.path, glob

def main():
    infile=open('ziplist','r')
    for line in infile:
        zfile=line.strip()
        dir1=zfile[0:-4]
        if os.path.isdir(dir1):
#            print(dir1,'already unzipped')
            f=dir1.split('_IGO_')
            core='JY_'+f[0][7:]
            bfile='{0}/{0}_star_anno.bam'.format(core)
            print(bfile)
            if os.path.isfile(bfile):
                print(bfile,'already exists')
            else:
                i=1
                while os.path.isfile('qs'+str(i)):
                    i=i+1
#                print(dir1+'/*R1_001.fastq.gz')
#                print(dir1+'/*R2_001.fastq.gz')
                fq1=glob.glob(dir1+'/*R1_001.fastq.gz')
                fqR1='../'+fq1[0]
                fq2=glob.glob(dir1+'/*R2_001.fastq.gz')
                fqR2='../'+fq2[0]
                qfile=open('qs'+str(i),'w')
                qfile.write('#!/bin/sh\n')
                qfile.write('#SBATCH --mem=40G\n')
                qfile.write('fqR1={0}\n'.format(fqR1))
                qfile.write('fqR2={0}\n'.format(fqR2))
                qfile.write('name={0}\n'.format(core))
                qfile.write('mkdir ${name}\n')
                qfile.write('echo ${name}\n')
                qfile.write('cp lrna_align_star_pe.sh ${name}\n')
                qfile.write('cp STAR ${name}\n')
                qfile.write('cd ${name}\n')
                qfile.write('ln -s ../out out\n')
                qfile.write('lrna_align_star_pe.sh unused ${fqR1} ${fqR2} ${name} 4 30 ${name}\n')
                qfile.close()
                os.system('sbatch -p all2 qs'+str(i))
        else:
            print('unzip {0}'.format(zfile))
#            os.system('unzip {0}'.format(zfile))
    infile.close()
    

main()
