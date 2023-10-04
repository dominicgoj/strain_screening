import pandas as pd

def importFile(file_name):
    pd_df = pd.read_csv("raw_data/"+file_name+".csv", encoding = 'unicode_escape', sep=";")
    return pd_df

def exportFile(file_name, pd_df):
    try:
        pd_df.to_csv("exports/"+file_name+".csv", sep=";")
        print("File %s.csv was exported to folder /exports/" %file_name)
    except Exception as e:
        print("Error was raised during export: ", e)

#importing file which contains information about which strain in database has expresses which relevant protein, clean and then export it
enzyme_list = importFile("ncc2").drop_duplicates(["EC"], keep="last").drop(["Gene", "Species", "Strain", "Comment", "Enzyme desired", "id"], axis=1).reset_index(drop=True)
export_file_name = input("Please enter name of exported enzyme file: ")
exportFile(export_file_name, enzyme_list)
