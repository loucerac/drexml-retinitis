#!/usr/bin/env python

import collections
import csv
import gzip
import io
import json
import os
import re
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import click
import pandas
import requests
from biothings_client import get_client


def collapse_list_values(row):
    for key, value in row.items():
        if isinstance(value, list):
            row[key] = "|".join(value)
    return row


def convert_uniprot_ids(uniprot_lst):
    mg = get_client("gene")
    # mg.getgenes(["uniprot:P00734"], fields='name,symbol,entrezgene,taxid', as_dataframe=True)
    genes_converted = mg.querymany(
        uniprot_lst,
        scopes="uniprot",
        fields="entrezgene,symbol",
        species="human",
        as_dataframe=True,
    )

    genes_converted = genes_converted.reset_index(names=["uniprot_id"])
    genes_converted = genes_converted.query("notfound!=True")[
        ["uniprot_id", "entrezgene", "symbol"]
    ]
    genes_converted = genes_converted.rename(columns={"entrezgene": "entrez_id"})

    return genes_converted


def build_drug_dataset(root):

    ns = "{http://www.drugbank.ca}"
    inchikey_template = (
        "{ns}calculated-properties/{ns}property[{ns}kind='InChIKey']/{ns}value"
    )
    inchi_template = (
        "{ns}calculated-properties/{ns}property[{ns}kind='InChI']/{ns}value"
    )

    rows = list()
    for i, drug in enumerate(root):
        row = collections.OrderedDict()
        assert drug.tag == ns + "drug"
        row["type"] = drug.get("type")
        row["drugbank_id"] = drug.findtext(ns + "drugbank-id[@primary='true']")
        row["name"] = drug.findtext(ns + "name")
        row["description"] = drug.findtext(ns + "description")
        row["groups"] = [
            group.text for group in drug.findall("{ns}groups/{ns}group".format(ns=ns))
        ]
        row["atc_codes"] = [
            code.get("code")
            for code in drug.findall("{ns}atc-codes/{ns}atc-code".format(ns=ns))
        ]
        row["categories"] = [
            x.findtext(ns + "category")
            for x in drug.findall("{ns}categories/{ns}category".format(ns=ns))
        ]
        row["inchi"] = drug.findtext(inchi_template.format(ns=ns))
        row["inchikey"] = drug.findtext(inchikey_template.format(ns=ns))

        # Add drug aliases
        aliases = {
            elem.text
            for elem in drug.findall(
                "{ns}international-brands/{ns}international-brand".format(ns=ns)
            )
            + drug.findall(
                "{ns}synonyms/{ns}synonym[@language='English']".format(ns=ns)
            )
            + drug.findall(
                "{ns}international-brands/{ns}international-brand".format(ns=ns)
            )
            + drug.findall("{ns}products/{ns}product/{ns}name".format(ns=ns))
        }
        aliases.add(row["name"])
        row["aliases"] = sorted(aliases)

        rows.append(row)

    rows = list(map(collapse_list_values, rows))

    columns = [
        "drugbank_id",
        "name",
        "type",
        "groups",
        "atc_codes",
        "categories",
        "inchikey",
        "inchi",
        "description",
    ]
    drugbank_df = pandas.DataFrame.from_dict(rows)[columns]

    return drugbank_df


def build_protein_df(root, use_groups=False):
    ns = "{http://www.drugbank.ca}"

    protein_rows = list()
    for _, drug in enumerate(root):
        drugbank_id = drug.findtext(ns + "drugbank-id[@primary='true']")
        for category in ["target", "enzyme", "carrier", "transporter"]:
            proteins = drug.findall("{ns}{cat}s/{ns}{cat}".format(ns=ns, cat=category))
            for protein in proteins:
                row = {"drugbank_id": drugbank_id, "category": category}
                row["organism"] = protein.findtext("{}organism".format(ns))
                row["known_action"] = protein.findtext("{}known-action".format(ns))
                actions = protein.findall("{ns}actions/{ns}action".format(ns=ns))
                row["actions"] = "|".join(action.text for action in actions)
                uniprot_ids = [
                    polypep.text
                    for polypep in protein.findall(
                        "{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier".format(
                            ns=ns
                        )
                    )
                ]

                if len(uniprot_ids) != 1:
                    if use_groups:
                        row["is_protein_group_target"] = True
                        row["uniprot_id"] = "|".join(
                            uniprot_id for uniprot_id in uniprot_ids
                        )
                    else:
                        continue
                else:
                    row["is_protein_group_target"] = False
                    row["uniprot_id"] = uniprot_ids[0]

                ref_text = protein.findtext(
                    "{ns}references[@format='textile']".format(ns=ns)
                )
                pmids = re.findall(r"pubmed/([0-9]+)", str(ref_text))
                row["pubmed_ids"] = "|".join(pmids)
                protein_rows.append(row)

    protein_df = pandas.DataFrame.from_dict(protein_rows)
    protein_df.uniprot_id = protein_df.uniprot_id.str.split("|")
    protein_df = protein_df.explode("uniprot_id")

    return protein_df


@click.command()
@click.argument("xml_path", type=click.Path(exists=True))
@click.argument("interim_folder", type=click.Path(exists=False))
@click.argument("final_folder", type=click.Path(exists=False))
@click.option("--use_groups", is_flag=True, default=False, help="number of greetings")
def main(xml_path, interim_folder, final_folder, use_groups):
    xml_path = Path(xml_path)
    interim_folder = Path(interim_folder)
    interim_folder.mkdir(exist_ok=True)

    basename = xml_path.with_name(xml_path.name.partition('.')[0])

    with gzip.open(xml_path) as xml_file:
        tree = ET.parse(xml_file)
    root = tree.getroot()

    drugbank_df = build_drug_dataset(root)
    # write drugbank tsv
    drugbank_df_path = interim_folder.joinpath(f"{basename}_drugs.tsv")
    drugbank_df.to_csv(drugbank_df_path, sep="\t", index=False)
    print(f"Wrote {drugbank_df_path}")

    protein_df = build_protein_df(root, use_groups=use_groups)
    protein_df_path = interim_folder.joinpath(f"{basename}_proteins.tsv")
    protein_df.to_csv(protein_df_path, sep="\t", index=False)
    print(f"Wrote {protein_df_path}")

    db_df = drugbank_df.merge(protein_df, how="inner")

    genes_df = convert_uniprot_ids(db_df.uniprot_id.unique())
    today_str = datetime.today().strftime("%Y%m%d")
    genes_df_fpath = interim_folder.joinpath(f"{basename}_mygene-{today_str}.tsv")
    genes_df.to_csv(genes_df_fpath, sep="\t", index=False)

    db_df = db_df.merge(genes_df, how="left")


if __name__ == "__main__":
    main()
