import sys
import pandas as pd


def merge():
    ont_file = open(sys.argv[1],'r',encoding='UTF-8')
    pacbio_file = open(sys.argv[2],'r',encoding='UTF-8')
    ont_data = pd.read_table(ont_file,sep='\t',names=['chr','stat','end','num','all','true_me','pro','site'])
    pacbio_data = pd.read_table(pacbio_file,sep='\t',
                    names=['chr','stat','end','pro','total','all','true_me','false_me','a','b'])
    ont_data = ont_data[['chr','stat','end']]
    pacbio_data = pacbio_data[['chr','stat','end']]
    overlap = pd.merge(ont_data,pacbio_data,on=['chr','stat','end'])  
    print('overlap number is ', len(overlap))


if __name__ == '__main__':
    if sys.argv[3] == 'overlap':
        merge()