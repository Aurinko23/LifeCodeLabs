from Bio import SeqIO
try:
    from Bio.SeqUtils import gc_fraction
    def GC(sequence):
        return 100 * gc_fraction(sequence, ambiguous="ignore")
except ImportError:
    # Older versions have this:
    from Bio.SeqUtils import GC

import argparse
import logging

logger = logging.getLogger('my_logger')
logger.setLevel(logging.DEBUG)  # Capture all levels

# Create formatters
formatter = logging.Formatter('%(levelname)s - %(message)s')

# Handler for DEBUG, INFO, and WARNING
info_handler = logging.FileHandler('info.log')
info_handler.setLevel(logging.DEBUG)
info_handler.addFilter(lambda record: record.levelno <= logging.WARNING)
info_handler.setFormatter(formatter)

# Handler for ERROR and CRITICAL
error_handler = logging.FileHandler('error.log')
error_handler.setLevel(logging.ERROR)
error_handler.setFormatter(formatter)

# Add handlers to the logger
logger.addHandler(info_handler)
logger.addHandler(error_handler)



def filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0):
    if quality_threshold < 0:
        logger.error(f"quality_threshold={quality_threshold}, can't be lower then 0, change to 0")
        quality_threshold = 0

    if not isinstance(gc_bounds, tuple):
        logger.info(f'You entered one number {gc_bounds}, it will be interpreted as the upper bound')
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


def int_or_range(value):
    try:
        if ',' in value:
            parts = value.split(',')
            if len(parts) != 2:
                raise ValueError
            return tuple(int(x.strip()) for x in parts)
        return int(value)
    except Exception:
        raise argparse.ArgumentTypeError(f"Expected int or int range (start,end), got '{value}'")


def get_args():
    parser = argparse.ArgumentParser(
                        prog='My wonderful parser',
                        description='This tool is needed to parse command line arguments',
                        epilog='Text at the bottom of help')

    parser.add_argument('-i', "--input", required=True, metavar='path', type=str, help='This is fastq format input file (takes value)')
    parser.add_argument('-g', "--gc", type=int_or_range, default=(0, 100), help="Either a single integer (e.g., 5) or a range in the form start,end (e.g., 0,100)")
    parser.add_argument('-l', "--lenb", type=int_or_range, default=(0,  2 ** 32), help="Either a single integer (e.g., 5) or a range in the form start,end (e.g., 0,  2 ** 32)")
    parser.add_argument('-t', "--thresh", type=int, default=0, help='Threshold value of average read quality for filtering (default is 0, phred33 scale)')
    parser.add_argument('-o', "--output", required=True, metavar='path', type=str, help='This is fastq format filtered output file')
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    # args.gggg
    # print(args)
    input_fastq=args.input
    output_fastq=args.output
    gc_bounds=args.gc
    length_bounds=args.lenb
    quality_threshold=args.thresh
    filter_fastq(input_fastq, output_fastq, gc_bounds, length_bounds, quality_threshold)
