import os

file_names = os.listdir(".")
file_names = list(filter(lambda x: x.endswith("_interface_mutations.txt"), file_names))

for file_name in file_names:

    with open(file_name, "r") as f:
        content = f.read()
        
    content.replace("\r", "")
    content = content.split("\n")
    content = list(filter(lambda x: len(x) != 0, content))
    
    for idx in range(len(content)):
        
        line = content[idx]
        content[idx] = line[1] + " " + line[2:-1] + " " + line[0] + " " + line[-1]
        
    content = "\n".join(content)
    
    fname_stem = ".".join(file_name.split(".")[:-1])
    with open(f"{fname_stem}_processed.txt", "w+") as f:
        f.write(content)