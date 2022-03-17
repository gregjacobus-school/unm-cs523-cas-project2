import random
class AminoAcid():
    '''
    There are 20 different amino acids.

    codon = combo of 3 letters from {A,G,C,U}
    amino acids: 1 - 6 codons
    '''
    def __init__(self):
        self.codon = None

    def gen_random(self):
        self.codon = [self._random_letter() for _ in range(3)]
        self._convert_to_bits()
        return self

    def _random_letter(self):
        return random.choice(['A', 'C', 'G', 'T'])

    def _convert_to_bits(self):
        '''
        A = 0
        C = 1
        G = 2
        T = 3
        '''
        self.bits = []
        for letter in self.codon:
            if letter == 'A':
                self.bits.append(0)
                self.bits.append(0)
            elif letter == 'C':
                self.bits.append(0)
                self.bits.append(1)
            elif letter == 'G':
                self.bits.append(1)
                self.bits.append(0)
            elif letter == 'T':
                self.bits.append(1)
                self.bits.append(1)
            else:
                raise ValueError(f"Invalid letter '{letter}' in codon")
        assert len(self.bits) == 6
