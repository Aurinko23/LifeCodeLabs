import re
from abc import ABC, abstractmethod

class BiologicalSequence(ABC):
    def __init__(self, bioseq):
        self.bioseq = bioseq
    
    def __len__(self):
        return len(self.bioseq)
    
    def __getitem__(self, index):
        return self.bioseq[index]
    

    @abstractmethod
    def verification(self):
        pass

class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, bioseq):
        super().__init__(bioseq)


    def verification(self):
        if not re.fullmatch(r"[ACUTGacutg]+", self.bioseq):
            raise TypeError("not valide NA")

        if "U" in self.bioseq or "u" in self.bioseq:
            if "T" in self.bioseq or "t" in self.bioseq:
                raise TypeError("not valide NA")


    def complement(self):
        table = str.maketrans("AaCcTtGgUu", "TtGgAaCcAa")
        return NucleicAcidSequence(self.bioseq.translate(table))
    

    def __reversed__(self):
        return NucleicAcidSequence(self.bioseq[::-1])
    

    def reverse_complement(self):
        return reversed(self.complement())
    
class RNASequence(NucleicAcidSequence):

    def __init__(self, bioseq):
        super().__init__(bioseq)

class DNASequence(NucleicAcidSequence):
    def __init__(self, bioseq):
        super().__init__(bioseq)

    def transcribe(self):
        return RNASequence(self.bioseq.replace("T", "U").replace("t", "u"))

class AminoAcidSequence(BiologicalSequence):
    def __init__(self, bioseq):
        super().__init__(bioseq)

    def verification(self):
        pass

    def is_basic(self):
        if not re.search(r"[KRHkrh]+", self.bioseq):
            print("There are no basic amino acids here")
        else:
            print("There are basic amino acids here")