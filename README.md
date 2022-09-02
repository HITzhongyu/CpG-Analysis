Here is the code for our paper

1. find ONT and PacBio data overlap 

python overlap.py  ONT-methyfile  PacBio-methyfile

2. GTF extract 

note: methy3-file just  extract  chr, stat,end three columns from methy-file 
	if you want to extract the number of UTR you need run this code

python GTF.py  methy3-file  UTR_annotation.gtf   UTR