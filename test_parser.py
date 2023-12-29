import xml.etree.ElementTree as ET
import os
import time
import biothings_client
import pickle
from collections.abc import Iterator
from collections import defaultdict

xml_path = "/Users/bailinzhang/Documents/Wu_Lab/Projects/HMDB_data"
xml_name = "metabolites_19014.xml"
# xml_name = "hmdb_metabolites.xml"
xml_file = os.path.join(xml_path, xml_name)

"""
# Too slow, not ideal for big files
tree = ET.parse("metabolites_microbe_example.xml")
root = tree.getroot()

output = {
    "_id": None,
    "hmdb": {
        "chebi": None,
        "pubchem_cid": None,
        "kegg_id": None,
        "food_id": None,
        "name": None,
        "description": None,
        "chemical_formula": None,
        "state": None,
        "associated_microbe": None

    }
}
for child in root:
    for element in child:
        if element.tag == "{http://www.hmdb.ca}name":
            name = element.text
            output["hmdb"]["name"] = name
        elif element.tag == "{http://www.hmdb.ca}description":
            description = element.text
            output["hmdb"]["description"] =description
        elif element.tag == "{http://www.hmdb.ca}chemical_formula":
            chem_form = element.text
            output["hmdb"]["chemical_formula"] = chem_form
        elif element.tag == "{http://www.hmdb.ca}state":
            state = element.text
            output["hmdb"]["state"] = state
        elif element.tag == "{http://www.hmdb.ca}kegg_id":
            kegg_id = element.text
            output["hmdb"]["kegg_id"] = kegg_id
        elif element.tag == "{http://www.hmdb.ca}foodb_id":
            food_id = element.text
            output["hmdb"]["food_id"] = food_id
        elif element.tag == "{http://www.hmdb.ca}chebi_id":
            chebi = element.text
            output["hmdb"]["chebi"] = chebi
            output["_id"] = chebi
        elif element.tag == "{http://www.hmdb.ca}pubchem_compound_id":
            pubchem_cid = element.text
            output["hmdb"]["pubchem_cid"] = pubchem_cid
        elif element.tag == "{http://www.hmdb.ca}ontology":
            for child_elem in element:
                if child_elem.tag == "{http://www.hmdb.ca}root":
                    for descendant in child_elem.iter("{http://www.hmdb.ca}descendant"):
                        term = descendant.find("{http://www.hmdb.ca}term")
                        if term.text == "Microbe":
                            descen_terms = descendant.findall(".//{http://www.hmdb.ca}term")
                            if descen_terms:
                                microbe = descen_terms[-1].text
                                print(microbe)
                                output["hmdb"]["associated_microbe"] = microbe

print(output)

'''
output:
{
    "_id":"17066",
    "hmdb":{
        "chebi":"17066",
        "pubchem_cid":"92135",
        "kegg_id":"C01089",
        "food_id":"FDB021869",
        "name":"3-Hydroxybutyric acid",
        "description":"3-Hydroxybutyric acid (CAS: 300-85-6), also known as beta-hydroxybutanoic acid, is a typical partial-degradation product of branched-chain amino acids (primarily valine) released from muscle for hepatic and renal gluconeogenesis. This acid is metabolized by 3-hydroxybutyrate dehydrogenase (catalyzes the oxidation of 3-hydroxybutyrate to form acetoacetate, using NAD+ as an electron acceptor). The enzyme functions in nervous tissues and muscles, enabling the use of circulating hydroxybutyrate as a fuel. In the liver mitochondrial matrix, the enzyme can also catalyze the reverse reaction, a step in ketogenesis. 3-Hydroxybutyric acid is a chiral compound having two enantiomers, D-3-hydroxybutyric acid and L-3-hydroxybutyric acid, and is a ketone body. Like the other ketone bodies (acetoacetate and acetone), levels of 3-hydroxybutyrate in blood and urine are raised in ketosis. In humans, 3-hydroxybutyrate is synthesized in the liver from acetyl-CoA and can be used as an energy source by the brain when blood glucose is low. Blood levels of 3-hydroxybutyric acid levels may be monitored in diabetic patients to look for diabetic ketoacidosis. Persistent mild hyperketonemia is a common finding in newborns. Ketone bodies serve as an indispensable source of energy for extrahepatic tissues, especially the brain and lung of developing mammals. Another important function of ketone bodies is to provide acetoacetyl-CoA and acetyl-CoA for the synthesis of cholesterol, fatty acids, and complex lipids. During the early postnatal period, acetoacetate (AcAc) and beta-hydroxybutyrate are preferred over glucose as substrates for the synthesis of phospholipids and sphingolipids in accord with requirements for brain growth and myelination. Thus, during the first two weeks of postnatal development, when the accumulation of cholesterol and phospholipids accelerates, the proportion of ketone bodies incorporated into these lipids increases. On the other hand, an increased proportion of ketone bodies is utilized for cerebroside synthesis during the period of active myelination. In the lung, AcAc serves better than glucose as a precursor for the synthesis of lung phospholipids. The synthesized lipids, particularly dipalmitoylphosphatidylcholine, are incorporated into surfactant, and thus have a potential role in supplying adequate surfactant lipids to maintain lung function during the early days of life (PMID: 3884391). 3-Hydroxybutyric acid is found to be associated with fumarase deficiency and medium-chain acyl-CoA dehydrogenase deficiency, which are inborn errors of metabolism. 3-Hydroxybutyric acid is a metabolite of Alcaligenes and can be produced from plastic metabolization or incorporated into polymers, depending on the species (PMID: 7646009, 18615882).",
        "chemical_formula":"C4H8O3",
        "state":"Solid",
        "associated_microbe":"Alcaligenes eutrophus" # Need to read paper to figure out the edge info
    }
}
'''
"""

# The difficulty of processing this data is that there are a lot of same element tags for different biological entities
# For example, descendant and descendants tags under different levels have different element text contents


start_time = time.time()


def strip_tag_namespace(tag) -> str:
    idx = tag.rfind("}")
    # rfind() method "not found" == -1
    if idx != -1:  # if idx is not "not found"
        tag = tag[idx + 1:]
    return tag


'''
microbe_list = ['Alcaligenes', 'Clostridium neopropionicum']
taxon = get_taxon_info(microbe_list)
'''


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
                                    # print(microbe_name.text)yield set(microbes_l)


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


def save_mapped_taxon_to_pkl(output_pkl):
    microbes_l = []
    microbes = get_all_microbe_names(xml_file)
    for microbe in microbes:
        microbes_l.append(microbe)

    unique_microbes = set(microbes_l)
    taxon = get_taxon_info(unique_microbes)

    check_pkl_file = os.path.isfile(output_pkl)
    if check_pkl_file is False:
        for taxon_d in taxon:
            with open(output_pkl, "wb") as handle:
                pickle.dump(taxon_d, handle, protocol=pickle.HIGHEST_PROTOCOL)
                return output_pkl
    else:
        pass


"""
mapped_taxon = save_mapped_taxon_to_pkl("hmdb_mapped_taxon.pkl")


with open("hmdb_mapped_taxon.pkl", "rb") as handle:
    taxon = pickle.load(handle)
    print(taxon)


taxon_l = ["Enterococcus", "Escherichia coli"]
test = get_taxon_info(taxon_l)
for microbe in test:
    print(microbe)



microbes = get_all_microbe_names(xml_file)
for microbe in microbes:
    taxon = get_taxon_info(microbe)
    for taxon_d in taxon:
        print(taxon_d)
"""

for event, elem in ET.iterparse(xml_file, events=("start", "end")):
    if event == 'end' and elem.tag.endswith('metabolite'):
        output = {"_id": None,
                  "name": None,
                  "xrefs": {},
                  "associated_microbes": [],
                  "associated_pathways": [],
                  "associated_diseases": [],
                  "associated_proteins": []}
        for metabolite in elem:
            tname = strip_tag_namespace(metabolite.tag)
            # TODO: add name, description, status, state, chemical_formula
            tname_list = ["chebi_id", "pubchem_compound_id", "chemspider_id", "drugbank_id",
                          "foodb_id", "kegg_id", "bigg_id", "smiles", "inchikey", "pdb_id"]
            if tname == "accession":
                if metabolite.text:
                    output["_id"] = metabolite.text
            elif tname in tname_list:
                if metabolite.text:
                    xref_ids = {tname: metabolite.text}
                    output["xrefs"].update(xref_ids)
            elif tname == "ontology":
                for descendant in metabolite.iter("{http://www.hmdb.ca}descendant"):
                    term = descendant.findall("{http://www.hmdb.ca}term")
                    if term and term[0].text == "Microbe":
                        microbe_descendants = descendant.findall(".//{http://www.hmdb.ca}term")
                        # The output of microbe_descendants = ['Microbe', 'Alcaligenes', 'Alcaligenes eutrophus']
                        microbe_names = [microbe.text for microbe in microbe_descendants[1:]]
                        unique_microbes = remove_duplicate_microbe(microbe_names)
                        with open("hmdb_mapped_taxon.pkl", "rb") as handle:
                            taxon = pickle.load(handle)
                            for microbe in unique_microbes:
                                if microbe in taxon:
                                    output[microbe] = output["associated_microbes"].append(taxon[microbe])
            elif tname == "biological_properties":
                for child in metabolite.iter("{http://www.hmdb.ca}biological_properties"):
                    biospec_loc = child.find("{http://www.hmdb.ca}biospecimen_locations")
                    if biospec_loc:
                        specimens = [specimen.text for specimen in biospec_loc]
                        output["sample_name"] = specimens

                    pathways = child.find("{http://www.hmdb.ca}pathways")
                    if pathways:
                        for pathway in pathways.iter("{http://www.hmdb.ca}pathway"):
                            pathway_dict = {"name": pathway.findtext("{http://www.hmdb.ca}name"),
                                            "kegg_map_id": pathway.findtext("{http://www.hmdb.ca}kegg_map_id"),
                                            "smpdb_id": pathway.findtext("{http://www.hmdb.ca}smpdb_id")}
                            output["associated_pathways"].append(pathway_dict)
            elif tname == "diseases":
                for diseases in metabolite.iter("{http://www.hmdb.ca}disease"):
                    if diseases:
                        disease_dict = {"name": diseases.findtext("{http://www.hmdb.ca}name"),
                                        "omim": diseases.findtext("{http://www.hmdb.ca}omim_id"),
                                        "pmid": []}
                        for ref in diseases.findall(".//{http://www.hmdb.ca}pubmed_id"):
                            if ref.text:
                                disease_dict["pmid"].append(ref.text)
                        output["associated_diseases"].append(disease_dict)
            elif tname == "protein_associations":
                for proteins in metabolite.iter("{http://www.hmdb.ca}protein"):
                    if proteins:
                        protein_dict = {"name": proteins.findtext("{http://www.hmdb.ca}name"),
                                        "uniprotkb": proteins.findtext("{http://www.hmdb.ca}uniprot_id")}
                        output["associated_proteins"].append(protein_dict)

        output["associated_pathways"] = remove_empty_values(output["associated_pathways"])
        output["associated_diseases"] = remove_empty_values(output["associated_diseases"])
        output["associated_proteins"] = remove_empty_values(output["associated_proteins"])
        output = {k: v for k, v in output.items() if v}
        print(output)

end_time = time.time()
total_time = end_time - start_time
print(f"The process takes {total_time} sec.")
