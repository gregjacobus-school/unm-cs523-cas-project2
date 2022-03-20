import random

from dataclasses import dataclass
from typing import Optional

@dataclass
class Antibody:
    bit_string: Optional[str] = ""
    base_sting: Optional[str] = ""
    germinal_center_id: Optional[str] = ""
    generation: Optional[str] = ""
    evolved_against: Optional[str] = ""

    def __iter__(self):
        for c in self.bit_string:
            yield c

    def __len__(self):
        return len(self.bit_string)

    def __post_init__(self):
        if self.bit_string == "":
            self.bit_string = self._get_random_bit_string()
            self.base_string = self._create_base_string_from_bit_string()
            self.evolved_against = self._pick_random_variant()

    def _get_random_bit_string(self):
        num_amino_acids = random.randint(6, 8)
        bitstring = ""
        for i in range(num_amino_acids*3):
            bitstring += str(random.randint(0, 1))
        return bitstring

    def _create_base_string_from_bit_string(self):
        """
        A = 0
        C = 1
        G = 2
        T = 3
        """
        base_string = ""
        for x, y in zip(*[iter(self.bit_string)] * 2):
            if x == "0" and y == "0":
                base_string += "A"
            elif x == "0" and y == "1":
                base_string += "C"
            elif x == "1" and y == "0":
                base_string += "G"
            else:
                base_string += "T"
        return base_string

    def _pick_random_variant(self):
        variants = ["Alpha", "Delta", "Omicron"]
        return random.choice(variants)
