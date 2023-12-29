from pathlib import Path
import os
import zipfile
import xml.etree.ElementTree as ET
import pickle
import biothings_client
from collections import defaultdict


def strip_tag_namespace(tag) -> str:
    idx = tag.rfind("}")
    # rfind() method "not found" == -1
    if idx != -1:  # if idx is not "not found"
        tag = tag[idx + 1:]
    return tag


def get_all_microbe_names(input_xml) -> str:
    if not os.path.exists("hmdb_mapped_taxon.pkl"):
        for event, elem in ET.iterparse(input_xml, events=("start", "end")):
            if event == 'end' and elem.tag.endswith('metabolite'):
                for metabolite in elem:
                    tagname = strip_tag_namespace(metabolite.tag)
                    if tagname == "ontology":
                        for descendant in metabolite.iter("{http://www.hmdb.ca}descendant"):
                            term = descendant.findall("{http://www.hmdb.ca}term")
                            if term and term[0].text == "Microbe":
                                microbe_descendants = descendant.findall(".//{http://www.hmdb.ca}term")
                                for microbe_name in microbe_descendants[1:]:
                                    # microbes_l.append(microbe_name.text)
                                    yield microbe_name.text
                                    # print(microbe_name.text) yield set(microbes_l)


def get_taxon_info(microbial_names):
    t = biothings_client.get_client("taxon")
    if not os.path.exists("hmdb_mapped_taxon.pkl"):
        taxon_info = t.querymany(microbial_names,
                                 scopes="scientific_name",
                                 fields=["_id", "scientific_name", "lineage", "parent_taxid", "rank"])
        # print(taxon_info)

        unique_taxon_d = {}
        taxon_d = defaultdict(list)
        for d in taxon_info:
            if "notfound" not in d:
                # taxid 2 is bacteria super kingdom on ncbi taxon browser
                if 2 in d["lineage"] and d["rank"] != "subgenus":
                    taxon_d[d["query"]].append((d["_score"]))

        # Take the highest score associated with the query microbial name
        max_score = dict([(name, max(score)) for name, score in taxon_d.items()])
        # print(max_score)

        for d in taxon_info:
            if d["query"] in max_score and d["_score"] == max_score[d["query"]]:
                unique_taxon_d[d["query"]] = {
                    "taxid": d["_id"],
                    "scientific_name": d["scientific_name"],
                    "lineage": d["lineage"],
                    "parent_taxid": d["parent_taxid"],
                    "rank": d["rank"]
                }
        # print(unique_taxon_d)
        yield unique_taxon_d


def save_mapped_taxon_to_pkl(input_xml, output_pkl):
    microbes_l = []
    microbes = get_all_microbe_names(input_xml)
    for microbe in microbes:
        microbes_l.append(microbe)

    unique_microbes = set(microbes_l)
    taxon = get_taxon_info(unique_microbes)

    check_pkl_file = os.path.exists(output_pkl)
    if not check_pkl_file:
        for taxon_d in taxon:
            with open(output_pkl, "wb") as handle:
                pickle.dump(taxon_d, handle, protocol=pickle.HIGHEST_PROTOCOL)
                return output_pkl
    else:
        pass


def remove_empty_values(new_association_list) -> list:
    filtered_association = []
    if new_association_list:
        for d in new_association_list:
            if any(d.values()):
                filtered_association.append({k: v for k, v in d.items() if v})
        new_association_list = filtered_association
        return new_association_list


def remove_duplicate_microbe(microbe_list) -> list:
    unique_microbes_l = []
    for microbe in microbe_list:
        is_unique = True
        for other_microbe in microbe_list:
            if microbe != other_microbe and microbe in other_microbe:
                is_unique = False
                break
        if is_unique:
            unique_microbes_l.append(microbe)
    return unique_microbes_l


path = Path.cwd()
file_path = os.path.join(path, "hmdb_metabolites.zip")
file_name = "hmdb_metabolites.xml"


with zipfile.ZipFile(file_path, "r") as zip_f:
    with zip_f.open(file_name) as xml_f:
        for event, elem in ET.iterparse(xml_f, events=("start", "end")):
            if event == 'end' and elem.tag.endswith('metabolite'):
                output = {"_id": None, "name": None}
                pathway_dict = {"name": [], "kegg_map_id": [], "smpdb_id": []}
                disease_dict = {"name": [], "omim": [], "pmid": []}
                protein_dict = {"name": [], "uniprotkb": []}
                for metabolite in elem:
                    tname = metabolite.tag.split("}")[-1]
                    tname_list = ["accession", "name", "chebi_id", "pubchem_compound_id",
                                  "chemical_formula", "chemspider_id", "drugbank_id",
                                  "foodb_id", "kegg_id"]
                    if tname == "accession":
                        output["_id"] = metabolite.text
                    elif tname in tname_list and not None:
                        output[tname] = metabolite.text if metabolite.text else None
                    elif tname == "ontology":
                        for descendant in metabolite.iter("{http://www.hmdb.ca}descendant"):
                            term = descendant.findall("{http://www.hmdb.ca}term")
                            if term and term[0].text == "Microbe":
                                microbe_descendants = descendant.findall(".//{http://www.hmdb.ca}term")
                                # The output of microbe_descendants = ['Microbe', 'Alcaligenes', 'Alcaligenes eutrophus']
                                microbe_names = [microbe.text for microbe in microbe_descendants[1:]]
                                output["associated_microbe"] = microbe_names
                    elif tname == "biological_properties":
                        for child in metabolite.iter("{http://www.hmdb.ca}biological_properties"):
                            biospec_loc = child.find("{http://www.hmdb.ca}biospecimen_locations")
                            if biospec_loc:
                                specimens = [specimen.text for specimen in biospec_loc]
                                output["sample_name"] = specimens

                            pathways = child.find("{http://www.hmdb.ca}pathways")
                            if pathways:
                                for pathway in pathways.iter("{http://www.hmdb.ca}pathway"):
                                    pathway_dict["name"].append(pathway.findtext("{http://www.hmdb.ca}name"))
                                    pathway_dict["kegg_map_id"].append(
                                        pathway.findtext("{http://www.hmdb.ca}kegg_map_id"))
                                    pathway_dict["smpdb_id"].append(pathway.findtext("{http://www.hmdb.ca}smpdb_id"))
                    elif tname == "diseases":
                        for diseases in metabolite.iter("{http://www.hmdb.ca}disease"):
                            if diseases:
                                disease_dict["name"].append(diseases.findtext("{http://www.hmdb.ca}name"))
                                disease_dict["omim"].append(diseases.findtext("{http://www.hmdb.ca}omim_id"))
                                for ref in diseases.findall(".//{http://www.hmdb.ca}pubmed_id"):
                                    disease_dict["pmid"].append(ref.text)
                    elif tname == "protein_associations":
                        for proteins in metabolite.iter("{http://www.hmdb.ca}protein"):
                            if proteins:
                                protein_dict["name"].append(proteins.findtext("{http://www.hmdb.ca}name"))
                                protein_dict["uniprotkb"].append(proteins.findtext("{http://www.hmdb.ca}uniprot_id"))

                pathway_dict = {key: [value for value in values if value] for key, values in pathway_dict.items() if
                                values != []}
                disease_dict = {key: [value for value in values if value] for key, values in disease_dict.items() if
                                values != []}
                protein_dict = {key: [value for value in values if value] for key, values in protein_dict.items() if
                                values != []}

                output["pathways"] = pathway_dict
                output["diseases"] = disease_dict
                output["associated_proteins"] = protein_dict
                print(output)







