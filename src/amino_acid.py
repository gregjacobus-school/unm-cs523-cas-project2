import random
class AminoAcid():
    def __init__(self):
        self.codon = None

    def gen_random(self):
        self.codon = [self._random_letter() for _ in range(3)]
        return self

    def _random_letter(self):
        return random.choice(['A', 'C', 'G', 'U'])
