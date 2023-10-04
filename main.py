import pandas as pd
import numpy as np

db = pd.read_csv("raw_data/ncc.csv", sep=";").drop("Species and Strain", axis=1) #drop species and strain column

def exportFile(file_name, pd_df):
    pd_df.to_csv("exports/"+file_name+".csv", sep=";")

negatives = db[db["desired"]=="no"] #select series of all strains which are not desired according to column

excluded_strains = db[db['Species'].isin(negatives['Species'])&db['Strain'].isin(negatives['Strain'])] #create df excluding all negatives

positive_strains = db[~db['id'].isin(excluded_strains.id)] #create df excluding all negatives

exportFile("negative_traits", negatives)
exportFile("excluded_strains", excluded_strains)

exclude_multiplicate_strains = positive_strains.drop_duplicates(subset=['Strain','EC','Compound'], keep='last').drop(["desired", "Gene"], axis=1) #duplicates with multiple genes having same strain ec and compound are dropped
hotlist = exclude_multiplicate_strains['Strain'] # count which strain occurs most often
hotlist = hotlist.to_frame()
hotlist[['label', 'nr']] = hotlist['Strain'].str.split("NCC", 1, expand=True)
hotlist_count = hotlist["nr"].value_counts(dropna=False)



hotlist_count_df = hotlist_count.to_frame()

hotlist_count_df.columns = ["Nr. of enzymes"]
hotlist_count_df['Strain'] = hotlist_count_df.index
hotlist_count_df[["Genus", "Species"]] = np.nan



ncc_tax_df = pd.read_csv("raw_data/ncc_tax.csv", sep=";")
ncc_tax_df.columns = ["Strain", "Genus", "Species"]
ncc_tax_strains = ncc_tax_df.Strain

for i in hotlist_count_df.Strain:
    if(i):
        selection = ncc_tax_df.loc[ncc_tax_df["Strain"] == int(i)]

        hotlist_count_df.at[i, "Genus"] = selection.Genus
        hotlist_count_df.at[i, "Species"] = selection.Species
        #print(hotlist_count_df.loc[[i]])
exportFile("positive_strains", exclude_multiplicate_strains)
exportFile("final_ranking", hotlist_count_df)
