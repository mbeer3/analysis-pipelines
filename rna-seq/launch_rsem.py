import os, sys, os.path, glob

def main():
    t=['NT','NANOG-pro2','OCT4-pro1']
    r=[1,2,3]
    for t1 in t:
        for r1 in r:
            core='JY_{0}-{1}'.format(t1,r1)
            bfile='{0}/{0}_star_anno.bam'.format(core)
            if os.path.isfile(bfile):
                print(bfile,'already exists')
                i=1
                while os.path.isfile('qrs'+str(i)):
                    i=i+1
                qfile=open('qrs'+str(i),'w')
                qfile.write('#!/bin/sh\n')
                qfile.write('#SBATCH --mem=35G\n')
                qfile.write('name={0}\n'.format(core))
                qfile.write('mkdir ${name}\n')


                qfile.write('echo ${name}\n')
                qfile.write('cp lrna_rsem_quantification.sh ${name}\n')
                qfile.write('cp STAR ${name}\n')
                qfile.write('cd ${name}\n')
                qfile.write('ln -s ../long-rna-seq-pipeline/dnanexus/quant-rsem/resources/usr/bin/rsem-calculate-expression rsem-calculate-expression\n')
                qfile.write('lrna_rsem_quantification.sh unused ${name}_star_anno.bam true unstranded 12345 1 ')
                qfile.close()
                os.system('sbatch -p all2 qrs'+str(i))
    

main()
