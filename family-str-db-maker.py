import random
from typing import Dict, List, Any
from scipy import stats
import numpy as np
import pandas as pd


profiles_df = pd.read_csv("profiles_3500.tsv", sep="\t")
print(profiles_df.head())

example_families_df = pd.read_csv("familyDBTest2.csv")


# Prep the frequency table for use
allele_freq_pop_stats_csv = "Revised-Allele-Freqs-PopStats-filled-in.csv"
allele_freq_table = pd.read_csv(allele_freq_pop_stats_csv)
allele_freq_table.fillna(value=0, inplace=True)

example_profile = {
    'Sample Name': '100F',
    'AMEL': ['X'],
    'CSF1PO': [10.0, 13.0],
    'D10S1248': [14.0, 14.0],
    'D12S391': [17.0, 18.0],
    'D13S317': [11.0, 10.0],
    'D16S539': [10.0, 12.0],
    'D18S51': [17.0, 13.0],
    'D19S433': [14.0, 14.0],
    'D1S1656': [16.0, 12.0],
    'D21S11': [32.2, 28.0],
    'D22S1045': [17.0, 12.0],
    'D2S1338': [17.0, 24.0],
    'D2S441': [11.0, 10.0],
    'D3S1358': [15.0, 15.0],
    'D5S818': [12.0, 12.0],
    'D7S820': [12.0, 9.0],
    'D8S1179': [16.0, 12.0],
    'FGA': [22.0, 23.0],
    'Penta_D': [9.0, 14.0],
    'Penta_E': [12.0, 18.0],
    'TH01': [9.3, 6.0],
    'TPOX': [8.0, 8.0],
    'vWA': [19.0, 17.0]
}

allele_lists = {
    'Sample Name': 'blank',
    'AMEL': [],
    'CSF1PO': [],
    'D10S1248': [],
    'D12S391': [],
    'D13S317': [],
    'D16S539': [],
    'D18S51': [],
    'D19S433': [],
    'D1S1656': [],
    'D21S11': [],
    'D22S1045': [],
    'D2S1338': [],
    'D2S441': [],
    'D3S1358': [],
    'D5S818': [],
    'D7S820': [],
    'D8S1179': [],
    'FGA': [],
    'Penta_D': [],
    'Penta_E': [],
    'TH01': [],
    'TPOX': [],
    'vWA': []
}

# print("Example Profile:\n", example_profile)


def make_full_profile(partner_name):
    partner_sex = partner_name[-1]
    mock_profile = {"Sample Name": partner_name + "_SP"}
    if partner_sex == 'F':
        mock_profile["AMEL"] = ["X", "Y"]
    else:
        mock_profile["AMEL"] = ["X"]
    for col in allele_freq_table:
        # print(col)
        if col == 'Allele':
            pass
        # elif col == "Sample Name":
        elif col == "AMEL":
            pass
        else:
            mock_profile[col] = random.choices(allele_freq_table['Allele'],
                                               allele_freq_table[col], k=2)
    return mock_profile


def fill_in_half_profile(half_profile):
    full_profile: Dict[Any, List[Any]] = {}
    for col in allele_freq_table:
        if col == 'Allele':
            pass
        elif col == 'Sample Name':
            pass
        elif col == 'AMEL':
            pass
        else:
            full_profile[col] = [half_profile[col]] + \
                                random.choices(allele_freq_table['Allele'],
                                               allele_freq_table[col], k=1)
    return full_profile


def split_profile(profile):
    """take a profile and split it in half
    paternal and maternal"""
    # print("in split profile")
    parent1 = {}
    parent2 = {}
    for locus in profile:
        if locus == 'Sample Name':
            # print("in name")
            parent1_name = profile[locus] + "_P1"
            parent2_name = profile[locus] + "_P2"
            # print(parent1_name, parent2_name)
            parent1[locus] = parent1_name
            parent2[locus] = parent2_name
        elif locus == "AMEL":
            # print("in AMEL")
            parent1[locus] = ["X"]
            parent2[locus] = ["X", "Y"]
        elif locus == "DYS391":
            pass
        else:
            alleles = profile[locus]
            position = random.randint(0, 1)
            parent1[locus] = alleles[position]
            parent2[locus] = alleles[1 - position]
    return [parent1, parent2]


def make_parent_profiles(child_profile):
    """take a profile and produce parent profiles"""
    child_name = child_profile["Sample Name"]
    # print("Child Name: ", child_name)
    half_parent_profiles = split_profile(child_profile)
    parent1 = fill_in_half_profile(half_parent_profiles[0])
    parent2 = fill_in_half_profile(half_parent_profiles[1])
    parent1["AMEL"] = ["X"]
    parent2["AMEL"] = ["X", "Y"]
    if child_name[-1] in ["F", "M"]:
        parent1["Sample Name"] = child_name + "_P1"
        parent2["Sample Name"] = child_name + "_P2"
    elif child_name.split("_")[1] in ["P1", "P2"]:
        child_number = child_name.split("_")[1][-1]
        parent1["Sample Name"] = \
            child_name.split('_')[0] + f"_GP{child_number}a"
        parent2["Sample Name"] = \
            child_name.split('_')[0] + f"_GP{child_number}b"
    return [parent1, parent2]


def make_child_profiles(parent1, parent2, number):
    """take two parent profiles and
    produce a number of children"""
    children = []
    while number > 0:
        child = {}
        for locus in parent1:
            if locus == "Sample Name":
                if len(parent1[locus].split('_')) == 2:
                    if parent1[locus].split('_')[1] in ["P1", "P2"]:
                        child[locus] = parent1[locus].split('_')[0] + "_S"
                    elif parent1[locus].split('_')[1] in ["GP1a", "GP1b"]:
                        child[locus] = parent1[locus].split('_')[0] + "_1AU"
                    elif parent1[locus].split('_')[1] in ["GP2a", "GP2b"]:
                        child[locus] = parent1[locus].split('_')[0] + "_2AU"
                else:
                    child[locus] = parent1[locus] + "_C"
            elif locus == "AMEL":
                sex = random.choice(['X', 'Y'])
                # print(sex)
                if sex == 'X':
                    child[locus] = ['X']
                else:
                    child[locus] = ['X', sex]
            elif locus == "DYS391":
                pass
            else:
                allele1 = parent1[locus][random.randint(0, 1)]
                allele2 = parent2[locus][random.randint(0, 1)]
                child[locus] = [allele1, allele2]
        children.append(child)
        number -= 1
    return children


def make_family(starting_profile):
    """Makes a family based on one donor profile should make
    parents, grandparents, aunts, uncles, siblings, cousins, children
    start with making parents -> grandparents -> siblings of parents ->
    siblings of donor -> children of donor"""
    family = []
    family.append(starting_profile)
    parents = make_parent_profiles(starting_profile)
    for parent in parents:
        family.append(parent)
    one_grandparents = make_parent_profiles(parents[0])
    for grandparent in one_grandparents:
        family.append(grandparent)
    two_grandparents = make_parent_profiles(parents[1])
    for grandparent in two_grandparents:
        family.append(grandparent)
    one_siblings = make_child_profiles(one_grandparents[0],
                                       one_grandparents[1], 2)
    for siblings in one_siblings:
        family.append(siblings)
    two_siblings = make_child_profiles(two_grandparents[0],
                                       two_grandparents[1], 2)
    for siblings in two_siblings:
        family.append(siblings)
    siblings = make_child_profiles(parents[0], parents[1], 2)
    for sibling in siblings:
        family.append(sibling)
    spouse = make_full_profile(starting_profile["Sample Name"])
    family.append(spouse)
    children = make_child_profiles(starting_profile, spouse, 3)
    for child in children:
        family.append(child)
    return family


def profiles_db_to_family_db(profiles_db):
    """Takes a profilesDB (Pandas DataFrame of profiles)
    and calls make family
    on each profile"""
    print("in profiledb to familydb")
    families_temp = profiles_df.iloc[0:0]
    for profile in profiles_db.to_dict(orient='Records'):
        # print(profile)
        for key in profile:
            #print(key)
            if key not in ["Sample Name", "AMEL", "DYS391"]:
                #print(profile[key])
                temp = sorted(profile[key].split(','))
                # print(temp)
                temp = [float(x) for x in temp]
                # print(temp)
                profile[key] = temp
                # print(profile[key])
                if len(profile[key]) == 1:
                    allele = profile[key][0]
                    if allele.is_integer():
                        allele = int(allele)
                    profile[key] = [allele, allele]
                else:
                    #if profile[key][0] == profile[key][1]:
                    allele1 = float(profile[key][0])
                    if allele1.is_integer():
                        allele1 = int(allele1)
                    allele2 = float(profile[key][1])
                    if allele2.is_integer():
                        allele2 = int(allele2)
                    profile[key] = sorted([allele1, allele2])
            elif key == "AMEL":
                temp = sorted(profile[key].split(','))
                #temp = ','.join(profile[key])
                #print("AMEL: ", temp)
                profile[key] = temp
        # print(profile)
        family_temp = make_family(profile)
        for family_profile in family_temp:
            #print(family_profile)
            families_temp = families_temp.append(family_profile, ignore_index=True)
    return families_temp




def family_db_to_csv(families_db):
    """Takes a familiesDB (Pandas DataFrame of profiles)
    and makes a CSV in the correct format
    this means sorting alleles so they are in numerical order
    and also make output as genotype (12,12)
    or phenotype (12)"""

    print("In families to csv")

    families_df = profiles_df.iloc[0:0]

    for profile in families_db.to_dict(orient='Records'):
        # print(profile)
        for key in profile:
            if key not in ["Sample Name", "AMEL", "DYS391"]:
                #print(key)
                #print(profile[key])
                temp = sorted(profile[key])
                # print(temp)
                profile[key] = temp
                # print(profile[key])
                if profile[key][0] == profile[key][1]:
                    allele = profile[key][0]
                    if isinstance(allele, float):
                        if allele.is_integer():
                            allele = int(allele)
                    profile[key] = str(allele)
                else:
                    allele1 = profile[key][0]
                    allele2 = profile[key][1]
                    if isinstance(allele1, float):
                        if allele1.is_integer():
                            allele1 = int(allele1)
                    if isinstance(allele2, float):
                        if allele2.is_integer():
                            allele2 = int(allele2)
                    # print(allele1, allele2)
                    # print(str(allele1) + ',' + str(allele2))
                    profile[key] = str(allele1) + ',' + str(allele2)
            elif key == "AMEL":
                temp = profile[key]
                temp = ','.join(profile[key])
                profile[key] = temp
        # print(profile)
        families_df = families_df.append(profile, ignore_index=True)

    families_output = families_df.set_index("Sample Name")
    families_output.to_csv("familyDBTest2.csv")



def dataframetofrequency(df):
    # names
    tempallelesdictionary = allele_lists
    for name in df.columns:
        tempalleles = []
        print(name)
        if isinstance(df[name][0], str):
            print("string")
            for alleles in df[name]:
                alleleslist = alleles.split(',')
                print(alleleslist)
                tempalleles.extend(alleleslist)
        elif isinstance(df[name][0], list):
            print("list")
            for alleles in df[name]:
                print(alleles)
                tempalleles.extend(alleles)
        tempallelesdictionary[name] = tempalleles
    for item in tempallelesdictionary:
        temparray = np.array(tempallelesdictionary[item])
        frequencyarray = np.unique(temparray, return_counts=True)
        tempsum = np.sum(frequencyarray[1])
        topercentage = np.vectorize(lambda a, b: a / b)
        percentages = topercentage(frequencyarray[1], tempsum)
        print(tempsum)
        print(frequencyarray[1])
        print(f"{item}\n", frequencyarray, percentages)



        # print(df[name][0])



# example_family = make_family(example_profile)
# print("Making Family\n", example_family)

# families_df = profiles_db_to_family_db(profilesDF)

# print("Families DF\n", families_df)

# family_db_to_csv(families_df)

# dataframetofrequency(profilesDF)

# dataframetofrequency(families_df)

dataframetofrequency(profiles_df)
