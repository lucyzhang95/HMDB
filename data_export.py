from hmdb_metabolites_parser import load_hmdb_data

data = load_hmdb_data()
for obj in data:
    if "associated_microbes" in obj:
        print(obj)