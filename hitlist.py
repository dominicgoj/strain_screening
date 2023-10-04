import pandas as pd
import numpy as np
import traceback
import warnings

def importFile(file_name):
    pd_df = pd.read_csv("raw_data/"+file_name+".csv", encoding = 'unicode_escape', sep=";")
    return pd_df

def exportFile(file_name, pd_df):
    try:
        pd_df.to_csv("exports/"+file_name+".csv", sep=";")
        print("File %s.csv was exported to folder /exports/" %file_name)
    except Exception as e:
        print("Error was raised during export: ", e)

export_file_name = input("Please enter name of exported hitlist file: ")

warnings.filterwarnings("ignore", category=FutureWarning)
db = importFile("ncc2")
ncc_tax = importFile("ncc_tax")
list_of_individual_strains = db.drop_duplicates("Strain", keep="last")
strainlist_no_gene_duplicates = db.drop_duplicates(["Strain", "EC", "Substrate", "Product"], keep="last") #does not contain duplicate strains and compounds because genes are irrelevant
final_df = pd.DataFrame()


def calcHits():
    rows_list = []
    for i in list_of_individual_strains.Strain:

        selection = strainlist_no_gene_duplicates.loc[strainlist_no_gene_duplicates["Strain"] == i] #find each strain from unique strain list in general list, gene duplicates not allowed

        positives = selection.loc[selection["Enzyme desired"] == "yes"]
        negatives = selection.loc[selection["Enzyme desired"] == "no"]
        # select rows with positive and negative negative_traits
        positive_score = 0
        negative_score = 0
        for p in positives["Weight"]:
            positive_score += p
        for n in negatives["Weight"]:
            negative_score += n
        sum_score = positive_score - negative_score

        # sum up all positive and negative scores and subtract

        positive_hits = len(positives.index)
        negative_hits = len(negatives.index)
        #count how many positive negative hits for this strain
        if((positive_hits+negative_hits)!=0):
            ratio = positive_hits/(positive_hits+negative_hits)
        else:
            ratio = 0
        #calculate ratio of positives to negatives

        species = list_of_individual_strains.loc[list_of_individual_strains["Strain"] == i].Species
        dataset = {"Strain": i, "Positives": positive_hits, "Negatives": negative_hits, "P/N ratio": ratio, "Score": sum_score}
        rows_list.append(dataset)
        #final_df = pd.concat([final_df, cache_df], ignore_index=True)

    final_df = pd.DataFrame(rows_list)
    final_df = final_df.dropna(subset=['Strain'])
    
    return final_df

def searchKeys(final_df):
    for strain in final_df["Strain"]:

        #select genus and species from key list
        if(strain):
            try:
                enzyme_df = db.loc[db["Strain"] == strain].drop_duplicates(["EC"], keep="last")
                enzymes_values = enzyme_df["EC"] #get all EC numbers that strain has
                desired_enzymes = enzyme_df.loc[enzyme_df["Enzyme desired"] == "yes"]["EC"]
                undesired_enzymes = enzyme_df.loc[enzyme_df["Enzyme desired"] == "no"]["EC"]
                de = ""
                ude = ""
                ## iterating through all (un)desired enzymes and appending them to string for later insertion
                for i in desired_enzymes:
                    de = de + i + ", "
                for i in undesired_enzymes:
                    ude = ude + i + ", "


                label, number = strain.split("NCC")
                selected_tax = ncc_tax.loc[ncc_tax["Key"] == int(number)]
                genus, species = selected_tax["GENUS"].item(), selected_tax["SPECIES"].item() #use item to get value of column genus and species
                origin = selected_tax["ORIGIN"].item()
                found_id = selected_tax["FOUND"].item()
                found_dict = {0: "not meat", 1: "meat", 2: "fish"}
                #select row in final_df where strain is equal, insert genus and species
                insert = final_df.loc[final_df["Strain"] == strain]
                insert_index = insert.index
                final_df.loc[insert_index, "Genus"] = genus
                final_df.loc[insert_index, "Species"] = species
                final_df.loc[insert_index, "Origin"] = origin
                final_df.loc[insert_index, "Desired enzymes"] = de
                final_df.loc[insert_index, "Undesired enzymes"] = ude
                final_df.loc[insert_index, "Found on"] = found_dict[found_id]
                #this comes last because nan values are not found in dict
                #set value NaN to genus and species
                #sorting


            except Exception as e:
                #traceback.print_exc()
                label, number = ["", ""]
                continue

print("Calculating hits...")
final_df = calcHits()

print("Assigning hits to NCC database...")
searchKeys(final_df)

def renameColumns(df):
    df.rename(columns = {"Positives":"Positive Traits", "Negatives":"Negative Traits"}, inplace=True)
    return df
export_df = renameColumns(final_df)

exportFile(export_file_name , export_df)
