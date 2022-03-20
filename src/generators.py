from functools import partial
from src.amino_acid import AminoAcid
from src.constants import *
import random

def _gen_antigen(virus, antigen_len):
    start_ind = random.randint(0,len(virus)-antigen_len)
    antigen = virus[start_ind : start_ind + antigen_len]
    return antigen

def gen_antigen(virus):
    return partial(_gen_antigen, virus, ANTIGEN_LEN)

def gen_antibody():
    num_amino_acids = random.randint(6,8)
    amino_acids = [AminoAcid().gen_random() for _ in range(num_amino_acids)]
    bits = []
    for amino_acid in amino_acids:
        bits.extend(amino_acid.bits)
    return bits
