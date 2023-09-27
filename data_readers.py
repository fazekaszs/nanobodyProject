import os

from collections import namedtuple
from typing import Dict
from pathlib import Path

import pandas as pd


MutationKey = namedtuple("MutationKey", ["method", "protein", "mutant"])
MutationCollection = Dict[MutationKey, float]

MethodFiles = namedtuple("MethodFiles", ["path", "extension"])

READ_PATHS = {
    "mCSM-PPI2": MethodFiles(Path("data/mCSMPPI2"), ".csv"),
    "MutaBind2": MethodFiles(Path("data/MutaBind2/outputs"), ".tsv"),
    "SAAMBE-3D": MethodFiles(Path("data/SAAMBE/outputs"), ".txt")
}

TLC_TO_OLC = {
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


def read_mcsm_ppi2() -> pd.DataFrame:

    prot_and_mut, ddgs = list(), list()

    method_name = "mCSM-PPI2"
    input_path = READ_PATHS[method_name].path
    input_extension = READ_PATHS[method_name].extension

    file_names = os.listdir(input_path)
    file_names = list(filter(lambda x: x.endswith(input_extension), file_names))

    for file_name in file_names:

        protein = file_name.replace(input_extension, "")

        with open(input_path / file_name) as f:
            data = f.read()

        data = list(filter(lambda x: len(x) != 0, data.split("\n")))[1:]

        for line in data:

            line = line.split(",")

            if line[0] != "L":
                continue

            mutant = TLC_TO_OLC[line[1]] + line[0] + line[2] + TLC_TO_OLC[line[3]]
            ddg = -1 * float(line[5])

            prot_and_mut.append(protein + ":" + mutant)
            ddgs.append(ddg)

    return pd.DataFrame(ddgs, index=prot_and_mut, columns=[method_name, ])


def read_mutabind2() -> pd.DataFrame:

    prot_and_mut, ddgs = list(), list()

    method_name = "MutaBind2"
    input_path = READ_PATHS[method_name].path
    input_extension = READ_PATHS[method_name].extension

    file_names = os.listdir(input_path)
    file_names = list(filter(lambda x: x.endswith(input_extension), file_names))

    for file_name in file_names:

        protein = file_name.replace(input_extension, "")
        protein = "_".join(protein.split("_")[:-1])

        with open(input_path / file_name) as f:
            data = f.read()

        data = list(filter(lambda x: len(x) != 0 and not x.startswith("#"), data.split("\n")))[1:]

        for line in data:

            line = line.split("\t")

            if line[2] == "None":
                continue

            mutant = line[1][0] + line[0] + line[1][1:]
            ddg = float(line[2])

            prot_and_mut.append(protein + ":" + mutant)
            ddgs.append(ddg)

    return pd.DataFrame(ddgs, index=prot_and_mut, columns=[method_name, ])


def read_saambe_3d() -> pd.DataFrame:

    prot_and_mut, ddgs = list(), list()

    method_name = "SAAMBE-3D"
    input_path = READ_PATHS[method_name].path
    input_extension = READ_PATHS[method_name].extension

    file_names = os.listdir(input_path)
    file_names = list(filter(lambda x: x.endswith(input_extension), file_names))

    for file_name in file_names:

        protein = file_name.replace(input_extension, "")

        with open(input_path / file_name) as f:
            data = f.read()

        data = list(filter(lambda x: len(x) != 0, data.split("\n")))[1:]

        for line in data:

            line = line.split(" ")

            mutant = line[3] + line[1] + line[2] + line[4]
            ddg = float(line[5])

            prot_and_mut.append(protein + ":" + mutant)
            ddgs.append(ddg)

    return pd.DataFrame(ddgs, index=prot_and_mut, columns=[method_name, ])
