# %%
import pathlib
import shutil
import urllib.request

import dotenv
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def fdr(pvalues):
    """Benjamini-Hochberg FDR p-value correction for multiple hypothesis testing."""
    return multipletests(pvalues, alpha=0.05, method="fdr_bh")[1]



def main(db_path, atc_path):
    db_path = pathlib.Path(db_path)
    atc_path = pathlib.Path(atc_path)
    np.random.seed(42)

    project_root = pathlib.Path(dotenv.find_dotenv()).absolute().parent
    data_folder = project_root.joinpath("data")
    raw_folder = data_folder.joinpath("raw")
    final_folder = data_folder.joinpath("final")
    results_folder = project_root.joinpath("results")
    tables_folder = results_folder.joinpath("tables")
    tables_folder.mkdir(parents=True, exist_ok=True)


    # %%
    atc_url = "https://raw.githubusercontent.com/fabkury/atcd/master/WHO%20ATC-DDD%202021-12-03.csv"
    
    if not atc_path.exists():
        with urllib.request.urlopen(atc_url) as response, open(atc_path, "wb") as out_file:
            shutil.copyfileobj(response, out_file)

    atc_code_name = pd.read_csv(atc_path, usecols=["atc_code", "atc_name"])
    atc_code_name.head()


    # %%
    shap_selection_df = pd.read_csv(
        results_folder.joinpath("ml", "shap_selection_symbol.tsv"), sep="\t", index_col=0
    )
    drugbank_df = pd.read_csv(
        final_folder.joinpath(db_path), sep="\t"
    ).assign(
        is_selected=lambda x: x.symbol_id.isin(
            shap_selection_df.columns[shap_selection_df.any()]
        )
    )


    # %%
    atc_level_to_len = {1: 1, 2: 3, 3: 4, 4: 5}


    # %%
    results = []
    ora_min_len = 3

    for atc_level, level_len in atc_level_to_len.items():
        tmp_df = (
            drugbank_df.loc[:, ["drugbank_id", "atc_codes", "is_selected"]]
            .dropna()
            .drop_duplicates()
            .groupby("drugbank_id")
            .agg({"atc_codes": list, "is_selected": "any"})
            .reset_index()
            .explode("atc_codes")
            .assign(atc_codes=lambda x: x.atc_codes.str.split("|"))
            .explode("atc_codes")
            .assign(atc_codes=lambda x: x.atc_codes.str[:level_len])
            .drop_duplicates()
        )
        atc_dict = (
            tmp_df[["drugbank_id", "atc_codes"]]
            .groupby("atc_codes")
            .agg({"drugbank_id": list})
            .to_dict()["drugbank_id"]
        )
        background = tmp_df.drugbank_id.unique()
        drugs_selected = tmp_df.drugbank_id[tmp_df.is_selected].unique()

        ora_dict = {}

        for atc_code in atc_dict.keys():
            drug_list_in_atc = np.unique(atc_dict[atc_code])
            n_drug_list_in_atc = drug_list_in_atc.size
            if n_drug_list_in_atc < ora_min_len:
                print(f"Ignore ATC codes with less than {ora_min_len} drugs")
            else:
                selected_drugs_in_atc = np.intersect1d(drugs_selected, drug_list_in_atc)
                selected_drugs_notin_atc = np.setdiff1d(drugs_selected, drug_list_in_atc)

                drugs_in_atc_not_selected = np.intersect1d(
                    drug_list_in_atc, np.setdiff1d(background, drugs_selected)
                )
                drugs_notin_atc_not_selected = np.setdiff1d(
                    np.setdiff1d(background, drugs_selected), drug_list_in_atc
                )

                contingency_table = np.array(
                    [
                        [selected_drugs_in_atc.size, drugs_in_atc_not_selected.size],
                        [selected_drugs_notin_atc.size, drugs_notin_atc_not_selected.size],
                    ]
                )

                odds_ratio, pvalue = stats.fisher_exact(
                    contingency_table, alternative="greater"
                )
                ora_dict[atc_code] = {
                    "ora_pval": pvalue,
                    "ora_unconditional_odds_ratio": odds_ratio,
                }

        this_results = (
            pd.DataFrame(ora_dict)
            .T.reset_index(names=["atc_code"])
            .assign(ora_bylevel_pval_adj=lambda x: fdr(x.ora_pval))
            .assign(atc_level=atc_level)
            .merge(atc_code_name, how="left")
        )

        results.append(this_results)


    # %%
    col_order = [
        "atc_code",
        "atc_level",
        "atc_name",
        "ora_unconditional_odds_ratio",
        "ora_pval",
        "ora_bylevel_pval_adj",
        "ora_pval_adj",
    ]

    results = (
        pd.concat(results, axis=0, ignore_index=True)
        .assign(ora_pval_adj=lambda x: fdr(x.ora_pval))
        .loc[:, col_order]
    )

    results.to_csv(
        tables_folder.joinpath("selected_drugs_atc_ora.tsv"),
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    import sys
    _, db_path, atc_path = sys.argv
    print(sys.argv)
    main(db_path, atc_path)
