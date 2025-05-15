import os
import pickle
import json

import biothings_client
import pandas as pd
from ete3 import NCBITaxa

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

    def save_all_hmdb_data_to_json(self, output_path, node_str=None):
        if node_str is None:
            data_to_save = [obj for obj in self.data]
        else:
            data_to_save = [obj for obj in self.data if node_str in obj]
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(data_to_save, f, ensure_ascii=False, indent=4)
        return data_to_save


class DataManipulation:
    def __init__(self):
        with open("data/hmdb_microbe_metabolites.pkl", "rb") as handle:
            self.hmdb_data = pickle.load(handle)
            print(len(self.hmdb_data))

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
    data_to_save = [obj for obj in load_hmdb_data()]
    with open("data/hmdb_metabolites_full.json", "w", encoding="utf-8") as f:
        json.dump(data_to_save, f, ensure_ascii=False, indent=4)
    with open("data/hmdb_metabolites_full.pkl", "wb") as handle:
        pickle.dump(data_to_save, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # data = SaveData()
    # if not os.path.exists("data/hmdb_metabolites_output.pkl"):
    #     hmdb_data_pkl = data.save_all_hmdb_data_to_pkl(
    #         "data/hmdb_metabolites_output.pkl", node_str=None
    #     )
    # if not os.path.exists("data/hmdb_microbe_metabolites.pkl"):
    #     hmdb_microbe_pkl = data.save_all_hmdb_data_to_pkl(
    #         "data/hmdb_microbe_metabolites.pkl", node_str="associated_microbes"
    #     )
    # hmdb_microbe_json = data.save_all_hmdb_data_to_pkl(
    #     "data/hmdb_microbe_metabolites.json", node_str="associated_microbes"
    #         )

    manipulation = DataManipulation()
    # microbe_metabolite = manipulation.get_node_pair(
    #     node_str="associated_microbes",
    #     key="scientific_name",
    #     filter_keys=["taxid", "rank", "lineage"],
    # )
    # disease_metabolite = manipulation.get_node_pair(
    #     node_str="associated_diseases", key="name", filter_keys=["omim"]
    # )
    # pathway_metabolite = manipulation.get_node_pair(
    #     node_str="associated_pathways", key="name", filter_keys=["smpdb_id", "kegg_map_id"]
    # )
    # protein_metabolite = manipulation.get_node_pair(
    #     node_str="associated_proteins", key="name", filter_keys=["uniprotkb"]
    # )
    # biospecimen_metabolite = manipulation.get_node_pair(
    #     node_str="biospecimen_samples", key="biospecimen_sample", filter_keys=None
    # )

    # df_microbe_metabolite = export_metabolite_pair_to_csv(input_pair_list=microbe_metabolite,
    #                                                       output_file="data/hmdb_microbe_metabolite.csv",
    #                                                       int_row="taxid", "lineage")
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

    # microbe_names = []
    # for d in microbe_metabolite:
    #     microbe_names.append(d["scientific_name"])
    #
    # taxids_query = []
    #
    # ncbi = NCBITaxa()
    # ncbi.update_taxonomy_database()
    #
    # # get the taxid for each scientific name
    # get_taxid = ncbi.get_name_translator(set(microbe_names))
    # for _name, taxid in get_taxid.items():
    #     for _id in taxid:
    #         taxids_query.append(_id)
    #
    # # query the taxid to get the lineage taxids
    # t = biothings_client.get_client("taxon")
    # query = t.gettaxa(taxids_query, fields=["scientific_name", "lineage", "rank"])
    #
    # # lineage_d example: {1526: [1526, 1506553, 186803, 186802, 186801, 1239, 1783272, 2, 131567, 1]}
    # lineage_d = {int(t["query"]): t["lineage"] for t in query if "notfound" not in t.keys()}
    # lineage_l = [taxid for k, v in lineage_d.items() for taxid in v]
    #
    # # query all lineage taxids to get their scientific name and rank
    # lineage_query = t.gettaxa(set(lineage_l), fields=["scientific_name", "rank"])
    # # lineage_name_d example: {2: {'query': '2', '_id': '2', '_version': 1, 'rank': 'superkingdom', 'scientific_name': 'bacteria'}}
    # lineage_name_d = {int(t["query"]): t for t in lineage_query if "notfound" not in t.keys()}
    #
    # all_mapped_lineage = {}
    # # check if each taxid in the lineage_d value list in lineage_name_d
    # for taxid, lineage in lineage_d.items():
    # # mapped_lineage example: {[clostridium] aminophilum: [{'species': '[clostridium] aminophilum'},
    # # {'genus': 'lachnoclostridium'}, {'family': 'lachnospiraceae'},
    # # {'order': 'eubacteriales'}, {'class': 'clostridia'}, {'phylum': 'bacillota'},
    # # {'clade': 'terrabacteria group'}, {'superkingdom': 'bacteria'}]}
    #     mapped_lineage = [
    #         {lineage_name_d[_id]["rank"]: lineage_name_d[_id]["scientific_name"]}
    #         for _id in lineage
    #         if _id in lineage_name_d
    #         if lineage_name_d[_id]["rank"] != "no rank"
    #     ]
    #     all_mapped_lineage.update({lineage_name_d[taxid]["scientific_name"]: mapped_lineage})
    #
    # # print(microbe_metabolite)
    #
    # for metabolite_d in microbe_metabolite:
    #     if metabolite_d["scientific_name"] in all_mapped_lineage:
    #         metabolite_d["lineage_rank_name"] = all_mapped_lineage[metabolite_d["scientific_name"]]
    #
    # df = pd.json_normalize(microbe_metabolite)
    # print(df)
