def verification(na):
    for nucls in na:
        if "U" in nucls or "u" in nucls:
            if "T" in nucls or "t" in nucls:
                print("Error input")
                exit(0)


def transcribe(*args):
    new_na = []
    for nucls in args:
        for nucl in nucls:
            na_trans = nucl.replace("T", "U").replace("t", "u")
            new_na.append(na_trans)
    return new_na


def reverse(*args):
    new_na = []
    for nucls in args:
        for nucl in nucls:
            na_rever = nucl[::-1]
            new_na.append(na_rever)
    return new_na


def complement(*args):
    new_na = []
    for nucls in args:
        for nucl in nucls:
            change1 = (
                nucl.replace("A", "O")
                .replace("a", "o")
                .replace("C", "H")
                .replace("c", "h")
                .replace("T", "D")
                .replace("t", "d")
                .replace("G", "K")
                .replace("g", "k")
                .replace("U", "J")
                .replace("u", "j")
            )
            change2 = (
                change1.replace("O", "T")
                .replace("o", "t")
                .replace("H", "G")
                .replace("h", "g")
                .replace("D", "A")
                .replace("d", "a")
                .replace("K", "C")
                .replace("k", "c")
                .replace("J", "A")
                .replace("j", "a")
            )
            new_na.append(change2)
    return new_na


def reverse_complement(*args):
    new_na = complement(*args)
    reverse_na = reverse(new_na)
    return reverse_na
