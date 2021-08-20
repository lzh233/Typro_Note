# Rosalind

## Day1 

计算ATGC数量、转录、互补、反向互补，GC含量

```python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 09:06:06 2021
@author: liuzihao

transcription
T>A A>U G>C C>G
"""
class Transcription:
    def __init__(self):
        self.__tran_dic = {"T":"a","A":"U","G":"c","C":"g"}
        self.__comp_dic = {"T":"a", "A":"t","G":"c","C":"g"}
        self.__dna_base = ["A","T","G","C"]
    def get_tarn(self,seq):
        seq = str.upper(seq)
        if "U" in seq:
            print("Please input correct DNA sequence!")
        else:
            for base in self.__dna_base:
                seq = seq.replace(base,self.__tran_dic[base])
            return str.upper(seq)

    def get_complementation(self,seq,reverse = False):
        seq = str.upper(seq)
        if "U" in seq:
            print("Please input correct DNA sequence!")
        else:
            for base in self.__dna_base:
                seq = seq.replace(base,self.__comp_dic[base])
            if reverse == True :
                seq = seq[::-1]
            return str.upper(seq) 
    def dna_coount(self,seq):
        count = {}
        for base in self.__dna_base:
            count[base] = seq.count(base) 
        return count
    def cotent_GC(self,seq):
        seq = str.upper(seq)
        return ((seq.count("G") + seq.count("C"))/len(seq))
            
seq = "ATTGGGCCCC"
def main():
    sequence = Transcription()
    tran = sequence.get_tarn(seq = seq)
    comp = sequence.get_complementation(seq = seq,reverse=False)
    re_comp = sequence.get_complementation(seq = seq,reverse=True)
    cotent_GC = sequence.cotent_GC(seq = seq)
    counts = sequence.dna_coount(seq = seq)
    print(f"\nTranscription \n{tran}")
    print(f"Complementation \n{comp}")
    print(f"Reverse complementation \n{re_comp}\n")
    print(f"DNA Nucleotides \n{counts}\n")
    print(f"GC cotene \n{round(cotent_GC*100,3)}%\n")

if __name__ == '__main__':
    main()
##-------------------output---------------------##
Transcription
UAACCCGGGG
Complementation
TAACCCGGGG
Reverse complementation
GGGGCCCAAT

DNA Nucleotides
{'A': 1, 'T': 2, 'G': 3, 'C': 4}

GC cotene
70.0%

```

