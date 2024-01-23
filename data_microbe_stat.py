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


def get_node_pair(node_str: str, key: str, filter_keys: list | None) -> list:
    output_pair = []
    with open("hmdb_metabolites_output.pkl", "rb") as handle:
        data = pickle.load(handle)
        for hmdb_d in data:
            if node_str in hmdb_d:
                for obj in hmdb_d[node_str]:
                    if isinstance(obj, dict):
                        pair_d = {
                            "metabolite_name": hmdb_d["name"],
                            "chemical_formula": hmdb_d["chemical_formula"],
                            "smiles": hmdb_d["xrefs"]["smiles"],
                            "inchikey": hmdb_d["xrefs"]["inchikey"],
                            key: obj[key]
                        }
                        for filter_key in filter_keys:
                            if filter_key in obj:
                                pair_d[filter_key] = obj[filter_key]
                        output_pair.append(pair_d)
                    else:
                        for item in hmdb_d[node_str]:
                            pair_d = {
                                "metabolite_name": hmdb_d["name"],
                                "chemical_formula": hmdb_d["chemical_formula"],
                                "smiles": hmdb_d["xrefs"]["smiles"],
                                "inchikey": hmdb_d["xrefs"]["inchikey"],
                                key: item
                            }
                            output_pair.append(pair_d)
    return output_pair


def export_metabolite_pair_to_csv(input_pair_list, output_file, int_row):
    df = pd.json_normalize(input_pair_list)
    if int_row is not None:
        df[[int_row]] = df[[int_row]].apply(convert_to_int)
    df.to_csv(output_file, index=False)
    return df


microbe_metabolite = get_node_pair(node_str="associated_microbes", key="scientific_name", filter_keys=["taxid", "rank"])
disease_metabolite = get_node_pair(node_str="associated_diseases", key="name", filter_keys=["omim"])
pathway_metabolite = get_node_pair(node_str="associated_pathways", key="name", filter_keys=["smpdb_id", "kegg_map_id"])
protein_metabolite = get_node_pair(node_str="associated_proteins", key="name", filter_keys=["uniprotkb"])
biospecimen_metabolite = get_node_pair(node_str="biospecimen_samples", key="biospecimen_sample", filter_keys=None)
# print(microbe_metabolite)
print(f"microbe-metabolite: {len(microbe_metabolite)}")
# print(disease_metabolite)
print(f"disease-metabolite: {len(disease_metabolite)}")
# print(pathway_metabolite)
print(f"pathway-metabolite: {len(pathway_metabolite)}")
# print(protein_metabolite)
print(f"protein-metabolite: {len(protein_metabolite)}")
# print(biospecimen_metabolite)
print(f"biospecimen-metabolite: {len(biospecimen_metabolite)}")

df_microbe_metabolite = export_metabolite_pair_to_csv(input_pair_list=microbe_metabolite,
                                                      output_file="hmdb_microbe_metabolite.csv",
                                                      int_row="taxid", )
df_disease_metabolite = export_metabolite_pair_to_csv(input_pair_list=disease_metabolite,
                                                      output_file="hmdb_disease_metabolite.csv",
                                                      int_row=None)
df_pathway_metabolite = export_metabolite_pair_to_csv(input_pair_list=pathway_metabolite,
                                                      output_file="hmdb_pathway_metabolite.csv",
                                                      int_row=None)
df_protein_metabolite = export_metabolite_pair_to_csv(input_pair_list=protein_metabolite,
                                                      output_file="hmdb_protein_metabolite.csv",
                                                      int_row=None)
df_biospecimen_metabolite = export_metabolite_pair_to_csv(input_pair_list=biospecimen_metabolite,
                                                          output_file="hmdb_biospecimen_metabolite.csv",
                                                          int_row=None)



