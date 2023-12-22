#!/usr/bin/env python3

import sys

if len(sys.argv) != 2: #teste para veririficar o input
    print('Defina o nome de um e somente um arquivo fasta')
    exit(0)

multi = sys.argv[1] #nome do arquivo a ser lido

dit = {}
with open(multi, 'r') as file:
    for l in file:
        if l.startswith('>'):
            l = l.upper().strip()
            dit[l] = ''
            k = l
        else:
            dit[k] = l.upper().strip().replace('U','T')
            if len(dit[k]) % 3 != 0:
                print('Falha, os ORFs devem ser multiplos de 3','\n','Sequencia problematica:',dit[k],'Resto da divisao por 3:',len(dit[k])%3)
                exit(0)
            elif dit[k][0:3] not in ['ATG']:
                print('Falha, os ORFs devem comecar com um codon de inicio','\n','Sequencia problematica:',dit[k])
                exit(0)
            elif dit[k][-3:] not in ['TAG','TGA','TAA']:
                print('Falha, os ORFs devem terminar com um codon de termino','\n','Sequencia problematica:',dit[k])
                exit(0)

codontab = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}

chaves_dit = list(dit.keys())
valores_dit = list(dit.values())
valor_max = max(valores_dit, key=len)
indice_valor_max = valores_dit.index(valor_max)
seq = chaves_dit[indice_valor_max]

with open('ORF.fna','w') as fna:
    fna.write(seq+'_frame1_'+valor_max[0:3]+'_'+valor_max[-3:]+'\n')
    fna.write(valor_max)

codons = []
for i in range(0, len(valor_max), 3):
    codon = valor_max[i:i+3]
    codons.append(codon)

pep = []
for c in codons:
    pep += codontab[c]

amino = ''.join(pep)

with open('ORF.faa','w') as faa:
    faa.write(seq+'_frame1_'+amino[0]+'_'+amino[-1:])
    faa.write('\n')
    faa.write(amino)

import hashlib

with open("script_getORF.py", "r") as f:
  script = f.read()
md5sum = hashlib.md5(script.encode("utf-8")).hexdigest()
print('O md5sum do script eh:',md5sum)
