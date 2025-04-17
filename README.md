# LifeCodeLabs

## Project Description

LifeCodeLabs is a comprehensive software package for biological data analysis, specializing in working with nucleic acids: DNA and RNA. The project is intended for scientists and researchers who require a reliable tool for genetic data analysis.

The main goal of the project is to simplify complex sequence processing by providing users with a flexible set of tools. These tools allow for quality data filtering and perform basic operations with nucleic acids, such as transcription and translation of DNA, obtaining complementary and reverse complementary sequences for further research and experiments.

LifeCodeLabs is implemented in the Python programming language. The project is geared toward use in educational and scientific research where reliability and reproducibility of results are required.

## Project Structure

- `Abtract_sequences.py` —  tools for working with nucleic and amino acids.
- `fasstq_filter.py` —tools for filtering FASTQ data.


## Functionality

### FASTQ Filtering

The function filters sequences based on quality, GC content, and length.
The program accepts 5 arguments: `input_fastq`, `output_fastq`, `gc_bounds`, `length_bounds`, `quality_threshold`:
- `input_fastq` - a file containing FASTQ sequences with the following structure: sequence identifier, sequence, sequence identifier, quality.
- `gc_bounds` - GC content range (in percentages) for filtering (default is (0, 100)). You can pass a single number as an argument; in that case, the number will be used as the upper bound, and the program will filter out reads with quality less than the given number.
- `length_bounds` - length range for filtering (default is (0, 2**32)). You can pass a single number as an argument; in that case, the number will be used as the upper bound, and the program will filter out reads with length less than the given number.
- `quality_threshold` - threshold value of average read quality for filtering (default is 0, phred33 scale). Reads with an average quality across all nucleotides below the threshold are discarded.

As a result, the program returns a similar file consisting only of those sequences that meet all the conditions.

### Operations with nucleic and amino acids

The library is organized around an abstract base class and several concrete implementations:

- **BiologicalSequence (Abstract Base Class)**  
  Provides common functionality for all biological sequences such as getting the length and accessing sequence elements. It requires subclasses to implement a `verification()` method for validating the sequence.

- **NucleicAcidSequence (Derived Class)**  
  Represents a generic nucleic acid sequence (DNA or RNA). It includes:
  - **Verification:** Checks that the sequence contains only valid nucleotide characters (`A`, `C`, `U`/`T`, `G`), and ensures that RNA and DNA characters are not mixed.
  - **Complement Calculation:** Computes the complement of the sequence using standard nucleotide pairing.
  - **Reversal:** Implements methods for reversing the sequence and obtaining its reverse complement.

- **RNASequence (Derived from NucleicAcidSequence)**  
  Represents an RNA sequence, inheriting the properties and methods from `NucleicAcidSequence`.

- **DNASequence (Derived from NucleicAcidSequence)**  
  Represents a DNA sequence and includes:
  - **Transcription:** Converts a DNA sequence to its corresponding RNA sequence by replacing thymine (`T/t`) with uracil (`U/u`).

- **AminoAcidSequence (Derived Class)**  
  Represents an amino acid sequence and includes a method to check for the presence of basic amino acids (lysine, arginine, and histidine).

---

## Features

- **Sequence Verification:**  
  Ensures that nucleic acid sequences contain only valid characters and that DNA and RNA characters are not mixed.
  
- **Complement and Reverse Complement:**  
  Allows users to compute the complement of a nucleic acid sequence and then reverse it to obtain the reverse complement.
  
- **Transcription (DNA to RNA):**  
  Provides an easy way to transcribe a DNA sequence into its RNA counterpart.
  
- **Basic Amino Acid Detection:**  
  Checks if an amino acid sequence contains any basic amino acids (K, R, H).

## Usage Examples
```python
# Creating a DNA Sequence and Transcribing to RNA
from Abtract_sequences import DNASequence

dna_seq = DNASequence("ATGCGTAC")
dna_seq.verification()
rna_seq = dna_seq.transcribe()
print("Transcribed RNA sequence:", rna_seq.bioseq)
```

```python
# Computing Complement and Reverse Complement
from Abtract_sequences import DNASequence

dna_seq = DNASequence("ATGCGTAC")
complement_seq = dna_seq.complement()
print("Complement:", complement_seq.bioseq)
reverse_complement_seq = dna_seq.reverse_complement()
print("Reverse Complement:", reverse_complement_seq.bioseq)
```

```python
# Working with Amino Acid Sequences
from Abtract_sequences import AminoAcidSequence

aa_seq = AminoAcidSequence("ACDEFGHIKLMNPQRSTVWY")
aa_seq.is_basic()
```

## Authors
Evgeniya Tsymbalova, Saint Petersburg, 2025.

## Software Requirements
The project requires Python 3.9 or higher.
