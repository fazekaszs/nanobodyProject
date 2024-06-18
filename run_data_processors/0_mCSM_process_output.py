import os

filenames = os.listdir(".")
filenames = list(filter(lambda x: x.endswith(".csv"), filenames))

print(f"Found: {filenames}")

tlc_to_olc = {
    "GLY": "G",
    "ALA": "A",
    "VAL": "V",
    "ILE": "I",
    "LEU": "L",
    "PHE": "F",
    "SER": "S",
    "THR": "T",
    "TYR": "Y",
    "ASP": "D",
    "GLU": "E",
    "ASN": "N",
    "GLN": "Q",
    "HIS": "H",
    "TRP": "W",
    "LYS": "K",
    "ARG": "R",
    "CYS": "C",
    "MET": "M",
    "PRO": "P"
}

for filename in filenames:

    print(f"Starting {filename}")

    fname_stem = ".".join(filename.split(".")[:-1])

    with open(filename, "r") as f:
        content = f.read()
        
    content = content.split("\n")[1:]
    content = list(filter(lambda x: len(x) != 0 and x[0] != "#", content))
    
    for idx in range(len(content)):
        
        line = content[idx].split(",")
        line = tlc_to_olc[line[1]] + line[0] + line[2] + tlc_to_olc[line[3]] + "\t" + line[5]
        content[idx] = line
        
    content = "\n".join(content)
        
    with open("mCSM." + fname_stem + "_processed.tsv", "w+") as f:
        f.write(content)
