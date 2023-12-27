from zenodo_client import Zenodo
import shutil
from dotenv import find_dotenv
from pathlib import Path
import tarfile

import pystow
import requests
from urllib.parse import quote


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
    "RP_map_functions_MPC-annot.xlsx",
    "drug_actions_withSimplAction.csv"
]

drexml_files = [
    "expreset_Hinorm_gtexV8.rds.feather",
    "expreset_pathvals_gtexV8.rds.feather"    
]


def get_latest_record(record_id):
    """Get latest zenodo record ID from a given deposition identifier

    Parameters
    ----------
    record_id : str
        deposition identifier

    Returns
    -------
    str
        latest record ID
    
    """

    url = requests.get(f"https://zenodo.org/records/{record_id}", timeout=10).url
    return url.split("/")[-1]


def ensure_zenodo(name, record_id="6020480"):
    """Ensure file availability and download it from zenodo

    Parameters
    ----------
    name : str
        file name
    record_id : str
        deposition identifier

    Returns
    -------
    path : path-like
        PosixPath to downloaded file

    """

    record_id = get_latest_record(record_id)
    print(name, quote(name))
    url = f"https://zenodo.org/records/{record_id}/files/{quote(name)}?download=1"
    #url = f"{quote(url)}?download=1"
    print(url)

    path = pystow.ensure("drexml", "datasets", record_id, url=url)

    return path

def download_files(record, files, folder):
    record = str(record)

    paths = [ensure_zenodo(this_file, record) for this_file in files]

    for path in paths:
        print(path)
        if "localPDB" in path.name:
            this_file = tarfile.open(path.as_posix())
            this_file.extractall(folder)
        else:
            shutil.copyfile(path.as_posix(), folder.joinpath(path.name))


if __name__ == "__main__":
    download_files(7957438, rp_files, RAW_FOLDER)
    download_files(7737166, drexml_files, FINAL_FOLDER)
