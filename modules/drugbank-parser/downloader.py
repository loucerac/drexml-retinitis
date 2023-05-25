from zenodo_client import Zenodo
import shutil
from dotenv import find_dotenv
from pathlib import Path
import tarfile

DATA_REPOSITORY = Path(find_dotenv()).parent.joinpath("data", "raw")
DATA_REPOSITORY.mkdir(exist_ok=True, parents=True)

files = [
    "amendments_drugActions_drugbank-v050108.tsv",
    "ATC.csv.gz",
    "drugbank-v050108.xml.gz",
    "genes_drugbank-v050108_mygene-20230120.tsv",
    "genes_gtex-V8_mygene-20230120.tsv",
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
    "hp.obo-v1.2-20190906",
    "localPDB.tar.gz",
    "phenotype_annotation20190906.tab",
    "phenotype_to_genes20191010.txt",
    "physiological_paths.tsv",
    "physPathsAnnot.tsv",
    "WHO ATC-DDD 2021-12-03.csv",
    "RP_map_functions_MPC-annot.xlsx"
]

zenodo = Zenodo()

paths = [zenodo.download_latest("7957439", this_file) for this_file in files]

for path in paths: 
    if "localPDB" in path.name:
        this_file = tarfile.open(path.as_posix())
        this_file.extractall(DATA_REPOSITORY)
    else:
        shutil.copyfile(path.as_posix(), DATA_REPOSITORY.joinpath(path.name))
