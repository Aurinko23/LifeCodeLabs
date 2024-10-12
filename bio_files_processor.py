def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=""):

    if output_fasta == "":
        output_fasta = input_fasta
    one_line = {}
    with open(input_fasta, "r") as file:
        for line in file:
            if line.startswith(">"):
                key = line.strip()
                one_line[key] = ""
            else:
                one_line[key] += line.strip()

    result = []
    for name, seq in one_line.items():
        result.append(name)
        result.append(seq)
    result_string = "\n".join(result)
    with open(output_fasta, "w") as file:
        file.write(result_string)


def parse_blast_output(input_file, output_file):

    with open(input_file) as file:
        f = file.readlines()

    result_string = []
    for index in range(0, len(f)):
        if f[index].startswith("Sequences producing significant alignments"):
            result_string.append(f[index + 3])

    name_protein = []
    for line in result_string:
        split_1 = line.split("]")
        if len(split_1) > 1:
            split_1[0] = split_1[0] + "]"
        split_2 = split_1[0].split(".")
        name_protein.append(split_2[0])

    name_protein.sort()
    result_string = "\n".join(name_protein)
    with open(output_file, "w") as file:
        file.write(result_string)
