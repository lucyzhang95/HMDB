from pathlib import Path
import os
import zipfile
import xml.etree.ElementTree as ET

path = Path.cwd()
file_path = os.path.join(path, "hmdb_metabolites.zip")
file_name = "hmdb_metabolites.xml"

"""
with zipfile.ZipFile(metabolites_file_path, "r") as zip_f:
    with zip_f.open(metabolites_file_name) as xml_f:
        partial_xml = ""
        for line in xml_f:
            line = line.decode("utf-8")
            partial_xml += line.strip()
            if partial_xml.startswith("<metabolite>") and partial_xml.endswith("</metabolite>"):
                root = ET.fromstring(line)
                metabolite_output = {}
                for metabolite in root.findall(".//{http://www.hmdb.ca}metabolite"):
                    metabolite_output["_id"] = metabolite.findtext("{http://www.hmdb.ca}accession")
                    name = metabolite.findtext("{http://www.hmdb.ca}name")
                    description = metabolite.findtext("{http://www.hmdb.ca}description")
                    chebi = metabolite.findtext("{http://www.hmdb.ca}chebi_id")
                    hmdb_data = {
                        "chebi": chebi,
                        "name": name,
                        "description": description
                    }
                    metabolite_output["hmdb"] = hmdb_data
                    print(metabolite_output)
"""

with zipfile.ZipFile(file_path, "r") as zip_f:
    with zip_f.open(file_name) as xml_f:
        for event, elem in ET.iterparse(xml_f, events=("start", "end")):
            if event == 'end' and elem.tag.endswith('metabolite'):
                output = {"_id": None}
                for metabolite in elem:
                    tname = metabolite.tag.split("}")[-1]
                    if tname == "accession":
                        output["_id"] = metabolite.text
                    elif tname in ["accession", "name", "chebi_id"]:
                        output[tname] = metabolite.text if metabolite.text else None
                    elif tname == "ontology":
                        for descendant in metabolite.iter("{http://www.hmdb.ca}descendant"):
                            term = descendant.find("{http://www.hmdb.ca}term")
                            for child in descendant.iter("{http://www.hmdb.ca}Microbe"):
                                print(child.tag, child.text)
                            if term is not None and term.text == "Microbe":
                                microbe = descendant.findall(".//{http://www.hmdb.ca}term")
                                species = microbe[-1].text
                                if species != "Microbe":
                                    output["microbe"] = species
                print(output)







