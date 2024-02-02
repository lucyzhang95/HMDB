import os
import pickle

import pandas as pd

from hmdb_metabolites_parser import load_hmdb_data


class SaveData:
    def __init__(self):
        self.data = load_hmdb_data()

    def save_all_hmdb_data_to_pkl(self, output_path, node_str=None):
        if node_str is None:
            data_to_save = [obj for obj in self.data]
        else:
            data_to_save = [obj for obj in self.data if node_str in obj]
        with open(output_path, "wb") as handle:
            pickle.dump(data_to_save, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return data_to_save


class DataManipulation:
    def __init__(self):
        with open("data/hmdb_microbe_metabolites.pkl", "rb") as handle:
            self.hmdb_data = pickle.load(handle)

    def get_node_pair(self, node_str: str, key: str, filter_keys: list | None) -> list:
        output_pair = []
        for hmdb_d in self.hmdb_data:
            if node_str in hmdb_d:
                for obj in hmdb_d[node_str]:
                    pair_d = {
                        "metabolite_name": hmdb_d["name"],
                        "chemical_formula": hmdb_d["chemical_formula"],
                        "smiles": hmdb_d["xrefs"]["smiles"],
                        "inchikey": hmdb_d["xrefs"]["inchikey"],
                    }
                    if isinstance(obj, dict):
                        pair_d[key] = obj[key]
                        for filter_key in filter_keys:
                            if filter_key in obj:
                                pair_d[filter_key] = obj[filter_key]
                        output_pair.append(pair_d)
                    else:
                        for item in hmdb_d[node_str]:
                            pair_d[key] = item
                            output_pair.append(pair_d)
        print(f"Count of unique metabolite-{node_str}: {len(output_pair)}")
        return output_pair


# def export_metabolite_pair_to_csv(input_pair_list, output_file, int_row):
#     df = pd.json_normalize(input_pair_list)
#     if int_row is not None:
#         df[[int_row]] = df[[int_row]].apply(convert_to_int)
#     df.to_csv(output_file, index=False)
#     return df


if __name__ == "__main__":
    data = SaveData()
    if not os.path.exists("data/hmdb_metabolites_output.pkl"):
        hmdb_data_pkl = data.save_all_hmdb_data_to_pkl(
            "data/hmdb_metabolites_output.pkl", node_str=None
        )
    if not os.path.exists("data/hmdb_microbe_metabolites.pkl"):
        hmdb_microbe_pkl = data.save_all_hmdb_data_to_pkl(
            "data/hmdb_microbe_metabolites.pkl", node_str="associated_microbes"
        )

    manipulation = DataManipulation()
    microbe_metabolite = manipulation.get_node_pair(
        node_str="associated_microbes", key="scientific_name", filter_keys=["taxid", "rank"]
    )
    disease_metabolite = manipulation.get_node_pair(
        node_str="associated_diseases", key="name", filter_keys=["omim"]
    )
    pathway_metabolite = manipulation.get_node_pair(
        node_str="associated_pathways", key="name", filter_keys=["smpdb_id", "kegg_map_id"]
    )
    protein_metabolite = manipulation.get_node_pair(
        node_str="associated_proteins", key="name", filter_keys=["uniprotkb"]
    )
    biospecimen_metabolite = manipulation.get_node_pair(
        node_str="biospecimen_samples", key="biospecimen_sample", filter_keys=None
    )

    # df_microbe_metabolite = export_metabolite_pair_to_csv(input_pair_list=microbe_metabolite,
    #                                                       output_file="data/hmdb_microbe_metabolite.csv",
    #                                                       int_row="taxid", )
    # df_disease_metabolite = export_metabolite_pair_to_csv(input_pair_list=disease_metabolite,
    #                                                       output_file="data/hmdb_disease_metabolite.csv",
    #                                                       int_row=None)
    # df_pathway_metabolite = export_metabolite_pair_to_csv(input_pair_list=pathway_metabolite,
    #                                                       output_file="data/hmdb_pathway_metabolite.csv",
    #                                                       int_row=None)
    # df_protein_metabolite = export_metabolite_pair_to_csv(input_pair_list=protein_metabolite,
    #                                                       output_file="data/hmdb_protein_metabolite.csv",
    #                                                       int_row=None)
    # df_biospecimen_metabolite = export_metabolite_pair_to_csv(input_pair_list=biospecimen_metabolite,
    #                                                           output_file="data/hmdb_biospecimen_metabolite.csv",
