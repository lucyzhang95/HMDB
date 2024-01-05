# from hmdb_metabolites_parser import load_hmdb_data
import pandas as pd
import pickle


# all_data = []
# data = load_hmdb_data()
# for obj in data:
#     print(obj)
#     if "associated_microbes" in obj:
#         all_data.append(obj)
# with open("hmdb_metabolites_output.pkl", "wb") as handle:
#     pickle.dump(all_data, handle, protocol=pickle.HIGHEST_PROTOCOL)


def convert_to_int(val):
    return pd.to_numeric(val, errors='coerce').astype('Int64')


metabolite_microbe_pairs = []
with open("hmdb_metabolites_output.pkl", "rb") as handle:
    data = pickle.load(handle)
    df = pd.json_normalize(data)
    microbes = pd.json_normalize(df["associated_microbes"].explode())
    col_to_convert = ["taxid", "parent_taxid"]
    microbes[col_to_convert] = microbes[col_to_convert].apply(convert_to_int)
    microbes["rank"] = microbes["rank"].replace('', pd.NA)
    microbes["scientific_name"] = microbes["scientific_name"].str.capitalize()
    microbe_export = microbes.to_csv("hmdb_microbes.csv")

    for hmdb_d in data:
        for microbe_d in hmdb_d["associated_microbes"]:
            metabolite_microbe_d = {"scientific_name": microbe_d["scientific_name"].capitalize(),
                                    "name": hmdb_d["name"],
                                    "chemical_formula": hmdb_d["chemical_formula"],
                                    "smiles": hmdb_d["xrefs"]["smiles"],
                                    "inchikey": hmdb_d["xrefs"]["inchikey"]}
            if "taxid" in microbe_d:
                metabolite_microbe_d["taxid"] = microbe_d["taxid"]
            if "rank" in microbe_d:
                metabolite_microbe_d["rank"] = microbe_d["rank"]
            if "pubchem_cid" in hmdb_d["xrefs"]:
                metabolite_microbe_d["pubchem_cid"] = hmdb_d["xrefs"]["pubchem_cid"]
            if "monoisotopic_mw" in hmdb_d:
                metabolite_microbe_d["monoisotopic_mw"] = hmdb_d["monoisotopic_mw"]
            metabolite_microbe_pairs.append(metabolite_microbe_d)

microbe_metabolites_df = pd.json_normalize(metabolite_microbe_pairs)
microbe_metabolites_df[["taxid"]] = microbe_metabolites_df[["taxid"]].apply(convert_to_int)
microbe_metabolites_df.to_csv("hmdb_microbe_metabolite.csv")


