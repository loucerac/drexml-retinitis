from zenodo_client import Zenodo
import shutil
from dotenv import find_dotenv
from pathlib import Path
import tarfile

DATA_REPOSITORY = Path(find_dotenv()).parent.joinpath("data")
RAW_FOLDER = DATA_REPOSITORY.joinpath("raw")
RAW_FOLDER.mkdir(exist_ok=True, parents=True)
FINAL_FOLDER = DATA_REPOSITORY.joinpath("final")
FINAL_FOLDER.mkdir(exist_ok=True, parents=True)

rp_files = [
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

drexml_files = [
    "expreset_Hinorm_gtexV8.rds.feather",
    "expreset_pathvals_gtexV8.rds.feather"    
]

def download_files(record, files, folder):
    record = str(record)

    zenodo = Zenodo()

    paths = [zenodo.download_latest(record, this_file) for this_file in files]

    for path in paths: 
        if "localPDB" in path.name:
            this_file = tarfile.open(path.as_posix())
            this_file.extractall(folder)
        else:
            shutil.copyfile(path.as_posix(), folder.joinpath(path.name))


if __name__ == "__main__":
    download_files(7957439, rp_files, RAW_FOLDER)
    download_files(7737166, drexml_files, FINAL_FOLDER)
