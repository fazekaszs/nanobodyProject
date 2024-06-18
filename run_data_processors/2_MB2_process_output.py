import os

filenames = os.listdir(".")
filenames = list(filter(lambda x: x.endswith(".tsv") and not x.endswith("_processed.tsv"), filenames))

for filename in filenames:

    with open(filename, "r") as f:
        content = f.read()
        
    content = content.split("\n")
    content = list(filter(lambda x: len(x) != 0 and x[0] != "#", content))[1:]
    
    for idx in range(len(content)):
        
        new_line = content[idx].split("\t")[:3]
        new_line = [new_line[1][0] + new_line[0] + new_line[1][1:], new_line[2], ]
        new_line = "\t".join(new_line)
        
        content[idx] = new_line
        
    content = "\n".join(content)
    
    fname_stem = ".".join(filename.split(".")[:-1])
    with open("MB2." + fname_stem + "_processed.tsv", "w+") as f:
        f.write(content)