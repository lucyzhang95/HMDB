from pathlib import Path
import os
import zipfile
import xml.etree.ElementTree as ET

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







