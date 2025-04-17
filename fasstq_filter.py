from Bio import SeqIO
try:
    from Bio.SeqUtils import gc_fraction
    def GC(sequence):
        return 100 * gc_fraction(sequence, ambiguous="ignore")
except ImportError:
    # Older versions have this:
    from Bio.SeqUtils import GC

def filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0):

    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)

    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)

    s = SeqIO.parse(input_fastq, "fastq")

    filtered_records = []
    for record in s:
        if length_bounds[0] >= len(record.seq) or len(record.seq) >= length_bounds[1]:
            continue

        gc = GC(record.seq)
        if gc_bounds[0] >= gc or gc >= gc_bounds[1]:
            continue

        quality = sum(record.letter_annotations["phred_quality"])

        if quality/len(record.seq) < quality_threshold:
            continue

        filtered_records.append(record)
    
    SeqIO.write(filtered_records, output_fastq, "fastq")