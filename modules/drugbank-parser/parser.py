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
import pandas as pd
import requests
from biothings_client import get_client


def collapse_list_values(row):
    for key, value in row.items():
        if isinstance(value, list):
            row[key] = "|".join(value)
    return row


def convert_gene_ids(gene_ids, source="entrezgene", target="uniprot,symbol"):
    renamer = {
        "uniprot": "uniprot_id",
        "entrezgene": "entrez_id",
        "symbol": "symbol_id",
    }
    client = get_client("gene")
    genes_converted = client.querymany(
        gene_ids,
        scopes=source,
        fields=target,
        species="human",
        as_dataframe=True,
    )

    genes_converted = genes_converted.reset_index(names=[source])
    cols_query = genes_converted.columns.isin(["uniprot", "entrezgene", "symbol"])
    genes_converted = genes_converted.loc[:, cols_query]
    genes_converted = genes_converted.rename(columns=renamer)

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
    for _, drug in enumerate(root):
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
                f"{ns}international-brands/{ns}international-brand"
            )
            + drug.findall(f"{ns}synonyms/{ns}synonym[@language='English']")
            + drug.findall(f"{ns}international-brands/{ns}international-brand")
            + drug.findall(f"{ns}products/{ns}product/{ns}name")
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
    drugbank_df = pd.DataFrame.from_dict(rows)[columns]

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
                        f"{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier"
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

    protein_df = pd.DataFrame.from_dict(protein_rows)
    protein_df.uniprot_id = protein_df.uniprot_id.str.split("|")
    protein_df = protein_df.explode("uniprot_id")

    return protein_df


@click.group()
@click.option("--update/--no-update", default=False)
@click.option("--drugbank-version")
@click.option("--gtex-version")
@click.argument("interim-folder", type=click.Path(exists=False))
@click.argument("final-folder", type=click.Path(exists=False))
@click.pass_context
def main(ctx, drugbank_version, gtex_version, update, interim_folder, final_folder):
    ctx.ensure_object(dict)
    ctx.obj["DRUGBANK_VERSION"] = drugbank_version
    ctx.obj["GTEX_VERSION"] = gtex_version
    ctx.obj["UPDATE"] = update
    ctx.obj["INTERIM_FOLDER"] = Path(interim_folder)
    ctx.obj["FINAL_FOLDER"] = Path(final_folder)


@main.command()
@click.argument("path", type=click.Path(exists=True))
@click.option("--kind", type=click.Choice(["drugbank", "gtex"], case_sensitive=False))
@click.pass_context
def translate(ctx, path, kind):
    path = Path(path)
    basename = path.name.partition(".")[0]

    kind = kind.lower()
    if kind == "drugbank":
        data = pd.read_csv(path, sep="\t")
        ids = data["uniprot_id"].unique()
        this_source = "uniprot"
        this_target = "entrezgene"
    elif kind == "gtex":
        data = pd.read_feather(path)
        ids = data.columns[data.columns.str.startswith("X")].str.replace("X", "")
        this_source = "entrezgene"
        this_target = "symbol"

    genes_df = convert_gene_ids(ids, source=this_source, target=this_target)
    today_str = datetime.today().strftime("%Y%m%d")
    genes_df_fpath = ctx.obj["INTERIM_FOLDER"].joinpath(
        f"{basename}_mygene-{today_str}.tsv"
    )
    genes_df.to_csv(genes_df_fpath, sep="\t", index=False)
    print(f"Wrote {genes_df_fpath}")


@main.command()
@click.argument("drugbank-path", type=click.Path(exists=True))
@click.argument("drugbank-genes-path", type=click.Path(exists=True))
@click.argument("gtex-genes-path", type=click.Path(exists=True))
@click.argument("gtex-path", type=click.Path(exists=True))
@click.pass_context
def filter(ctx, drugbank_path, drugbank_genes_path, gtex_genes_path):
    data = pd.read_csv(drugbank_path, sep="\t")
    genes_drugbank = pd.read_csv(drugbank_genes_path, sep="\t")
    genes_gtex = pd.read_csv(gtex_genes_path, sep="\t")
    data = (
        data.merge(genes_drugbank, how="inner")
        .merge(genes_gtex, how="inner")
        .query(
            (
                'category=="target"'
                ' & groups.str.contains("^approved")'
                ' & (~groups.str.contains("withdrawn"))'
                ' & organism=="Humans"'
                ' & (known_action=="yes")'
            )
        )
        .sort_values("drugbank_id")
    )

    version = ctx.obj["DRUGBANK_VERSION"]
    genes_gtex[f"approved_{version}"] = genes_gtex.entrez_id.isin(data.entrez_id)


@main.command()
@click.argument("xml-path", type=click.Path(exists=True))
@click.option("--use-groups", is_flag=True, default=False, help="number of greetings")
@click.pass_context
def parse(ctx, xml_path, use_groups):
    xml_path = Path(xml_path)
    interim_folder = Path(ctx.obj["INTERIM_FOLDER"])
    interim_folder.mkdir(exist_ok=True)

    basename = xml_path.name.partition(".")[0]

    with gzip.open(xml_path) as xml_file:
        tree = ET.parse(xml_file)
    root = tree.getroot()

    drugbank_df = build_drug_dataset(root)

    protein_df = build_protein_df(root, use_groups=use_groups)

    db_df = drugbank_df.merge(protein_df, how="inner")
    db_df_path = ctx.obj["INTERIM_FOLDER"].joinpath(f"{basename}.tsv")
    db_df.to_csv(db_df_path, sep="\t", index=False)
    print(f"Wrote {db_df_path}")


if __name__ == "__main__":
    main()
