#!/usr/bin/env python
# coding: utf-8


#!/usr/bin/env python
import pandas as pd
import matplotlib as mp
import numpy as np
from scipy.stats import spearmanr
from collections import defaultdict
import operator
import seaborn as sp
import os


'''xl = pd.ExcelFile("Lung_Original_File.xlsx")
print(xl.sheet_names)
df1 = xl.parse('Queried Samples')'''


'''df = xl.parse('Queried Samples')
top5 = df.head()
print(top5)'''



'''df_dtypes = df.dtypes
df_shape = df.shape

print(df_dtypes)
print(df_shape)'''



'''df1= df.drop(columns = ['xrn', 'disease_group_term', 'local_met_status', 'chromosomal_sex', 'msi_status','bait_set','tTMB','rna_qc_status','IPAT_akt_pathway_call','AKT1_or_PIK3CA','PTEN_LOF','PTEN_het_del','PTEN_homo_del','PTEN_SV_LOH','study','mutation_load','pathologyPercentTumorNuclei','percent_genome_loh', 'non_human_variants','hlas','hla_loh_status','hla_loh_lost_alleles','hla_loh_status_imputed','sample_qc_status'])
#print(df1.head)'''


'''top5 = df1.head()
df = df1
print(top5)'''



'''print(df.dtypes)'''

#def setup1():
'''
The setup function creates the environment for the program to start. It includes loading the files and creating
the DataFrame in while the work is done
:return:
'''


#Removing the columns that aren't being actively used and we are naming this DataFrame as df1


if os.path.isfile('df1.csv'):

    df1 = pd.read_csv('df1.csv')
  #  df2 = pd.read_csv('df2.csv')
    print("Done importing")

else:
    # Import the excel file, look at the sheets and create a Data Frame with a certain sheet.

    xl = pd.ExcelFile("Lung_Original_File.xlsx")
    print(xl.sheet_names)
    df = xl.parse('Variants')

    # Assign variables to the dataframe's shape and types so that it's easier to print.
    df_dtypes = df.dtypes
    df_shape = df.shape

    # print(df_dtypes)
    print(df_shape)

    dfx = df.drop(columns=['xrn', 'tissue', 'chromosomal_sex', 'disease_ontology_term', 'disease_group_term', 'transcript',
                 'bait_set', 'fraction_reads', 'coding_type', 'zygosity', 'somatic/germline', 'study',
                 'nonhuman_reads_per_million', 'refSeq', 'altSeq', 'original_coding_type'])

    df1 = pd.DataFrame(dfx)
    df1.to_csv('df1.csv', index=False)

#return df1
#def setup2():

if os.path.isfile('df2.csv'):
    df2 = pd.read_csv('df2.csv')
    print("Done importing 2nd file")

else:

    xl = pd.ExcelFile("Lung_Original_File.xlsx")
    print(xl.sheet_names)
    dfy = xl.parse('Queried Samples')

    # Assign variables to the dataframe's shape and types so that it's easier to print.
    dfy_dtypes = dfy.dtypes
    dfy_shape = dfy.shape

    # print(df_dtypes)
    print(dfy_shape)

    df2 = pd.DataFrame(dfy)
    df2.to_csv('df2.csv', index=False)



#return df1,df2
#Creating a function that create dictionaries for CNAs, Short Variants and Rearrangement.
#def mutation_dictionary(df1):

'''
Creating a function that create dictionaries for CNAs, Short Variants and Rearrangement.
'''

def setup3():

    #Create empty dictionaries to store the values for CNA, Short Variant and rearrangement
    my_dict_CN = {}
    my_dict_SV = {}
    my_dict_RE = {}

    #iterrows is a method in Pandas that iterates over the rows in a given DataFrame(df1). Index and Row are returned from
    #iterrows. Index stores the index label of the current row being iterated and row is the value of each column

    for index, row in df1.iterrows():

        # Checks the column 'alteration_type' and looks for if it contains the word CN
        if row['alteration_type'] == 'CN':

        #Every dictionary has key and value. The key is unique and can have multiple values. Here the 'trf' which is essentially
        #the ORD number is being assigned as the key and the value of that is the "long_id" which contains data in this format
        #CDKN2A:deletion:0.0,  KEAP1:V155F:463G>T:0.0847, ALK:EML4 for CNA, SV and rearrangement

            key = row['trf']
            value = row['long_id']

        #Every row is being checked and if the key(ORD number) isn't present in the dictionary then the key is being added
        #in the dictionary and the value for the key "long id" is being assgined to that key.
            if key not in my_dict_CN:
                my_dict_CN[key] = []

        # Essentially my_dict_CN is being checked for the key and using .append the value ("long id") is being assigned.
            my_dict_CN[key].append(value)


        #Checking the row for they keyword SV.
        if row['alteration_type'] == 'SV':


        #KEAP1:V155F:463G>T:0.0847 is an example of the long id in SV. We are interested in the gene name and the protein
        #effect that the mutation has caused. We are ignoring the CDS effect and the %reads in this sample.

            key = row['trf']
        # We use a for loop to go through the long id and split it at every ",". In the above example the long_id gets
        #split into 4 different values. The first 2 elements are joined using ":" . In this case we would obtain KEAP1:V155F
        #We assign this new value to a variable called "gene_prot_effect".
            for value in row['long_id'].split(','):
                gene_prot_effect = ':'.join(value.split(':')[:2])


        #Same as above, if the key isn't present in the dictonary add it and then add the gene_prot_effect as value of key.
            if key not in my_dict_SV:
                my_dict_SV[key] = []

            my_dict_SV[key].append(gene_prot_effect)


        #In this case we are looking at 2 different rows. First we check if the alteration type has "RE" and if the does
        #then we check the column 'Var_Status' to see if the value is "known". If Yes, then the loop starts.

        if row['alteration_type'] == 'RE' and row['var_status'] == 'known':


        #Assigning the key as ORD number present in the column 'trf' and assigning it the value of 'long_id' present.
            key = row['trf']
            value = row['long_id']

            if key not in my_dict_RE:
                my_dict_RE[key] = []

            my_dict_RE[key].append(value)


    '''        
    #print(my_dict)
    first2pairs_RE = {k: my_dict_RE[k] for k in list(my_dict_RE)[:5]}
    print(first2pairs_RE)
    
    first2pairs_SV = {k: my_dict_SV[k] for k in list(my_dict_SV)[:5]}
    print(first2pairs_SV)
    
    first2pairs_CN = {k: my_dict_CN[k] for k in list(my_dict_CN)[:5]}
    print(first2pairs_CN)
    '''

    '''
    We create dictionaries called first2pairs_RE, SV and CN to a result that we obtain from running the codes above. 
    Running the code above we get 3 unique dictionaries with key value pairs. This is the original dictionary that 
    contains all of the info. The codes below creates 3 new dictionaries which is basically a copy of the dictionary above.
    We can then alter this new dictionary being fully confident that our original dictionary won't get changed. 
    Also this new dictionary is helpful if we want to maybe look at only the first 10 elements of each dictionary. 
    In this case we can work with a much smaller dictionary and it's faster. 
    '''

    # We create a new dictionary called first2pairs_RE.
    # list(my_dict_RE) creates a list of the keys present in the dictionary my_dict_RE
    # for k in list(my_dict_RE) creates a for loop that goes through the entire dictionary
    # K: my_dict_RE[k] looks for the 'k' which is the keys in the original dictionary and then adds the key:value pairs

    # Example: {'ORD-0599028-01': ['ALK:EML4']} is the first key:value pair in the original dictionary. When we run this code
    # k = ORD-0599028-01. It searches the my_dict_RE for ORD-0599028-01. The code becomes my_dict_RE[ORD-0599028-01] and
    # then adds this new key and it's associated value, which is 'ALK:EML4', to the new dictionary called first2pairs_RE

    # the len function is then used to calculate the length of the dictionary. The dictionary store a huge number of
    # key:value pairs so printing them might take a while and crashes the computer. Printing the length helps us know that
    #dictionary has been setup.

    first2pairs_RE = {k: my_dict_RE[k] for k in list(my_dict_RE)}
    #print(len(first2pairs_RE))

    print(" ")
    first2pairs_SV = {k: my_dict_SV[k] for k in list(my_dict_SV)}
    #print(len(first2pairs_SV))

    print(" ")
    first2pairs_CN = {k: my_dict_CN[k] for k in list(my_dict_CN)}
    #print(len(first2pairs_CN))

   # return (first2pairs_CN,first2pairs_RE,first2pairs_SV)


    #def ancestry_dictionary(df2):
    '''
     This function creates a new dictionary that has the ORD number as key and genetic ancestry as the value.
    :return:
    '''

    #Creating a new dictionary for genetic ancestry
    my_dict_ancestry = {}


    # Using iterrows to go through the df2 DataFrame and assigning keys as the items present in the column 'trf'. In this
    # case it's the ORD number. The value for the keys is assigned as the item found in "genetic_ancestry" column. The
    # values are either EAS(East Asian), SAS(South Asian), AFR(African), EUR (European), AMR(Mixed American) or the column
    # was blank and those column are given a value of nan.

    # df2 is the DataFrame created from the "Queried Sample" sheet of our excel file.

    for index, row in df2.iterrows():
        key = row['trf']
        value = row['genetic_ancestry']


        my_dict_ancestry[key] = value


    #Creates a new dictionary from the my_dict_ancestry above for same reason as creating secondary dictionaries for
    # rearrangement, short variant and CNA.

    print(" ")
    first2pairs_ancestry = {k: my_dict_ancestry[k] for k in list(my_dict_ancestry)}
    #print(len(first2pairs_ancestry))


    '''
    Since the dictionary only has 6 values of genetic ancestry. It can be helpful if we have a dictionary that has the
    genetic ancestry as the key and the ORD numbers as the values. So creating a new dictionary for that.
    
    defaultdict is a present in the module collections. You will get an error if you don't import collections. defualtlist
    is useful when we don't have a key. If we search for a non existing key like ASIA then it will return an empty list.
    '''

    my_dict_ancestry1 = defaultdict(list)

    # Goes through the original dictionary and assigns key value pairs to the new dictionary called my_dict_ancestry1.
    # The only difference is here the key and values are reversed.
    # E.g first2pairs_ancestry = {[TRF000121: AFR],[ORD093489: EUR], [ORD289832 :EUR]}
    #     first2pairs_ancestry1 = {[AFR: TRF000121], [EUR: 'ORD093489', 'ORD289832']}

    for key, value in my_dict_ancestry.items():
        my_dict_ancestry1[value].append(key)

    first2pairs_ancestry1 = {k: my_dict_ancestry1[k] for k in list(my_dict_ancestry1)}
    #print(first2pairs_ancestry1)


    # Provides a list of the various keys in my_dict_ancestry1 and assigns it to a variable called ancestry_list
    ancestry_list = list(my_dict_ancestry1)

    #print(ancestry_list)


    # We can now search the dictionary with keywords 'AFR', 'EUR' etc. It will then give us the list of samples that
    # contains the provided genetic ancestry. Again len is used so that we get the length. Removing len gives us the
    #complete list, which can be very long.

    '''
    print(len(my_dict_ancestry1['AFR']))
    print(len(my_dict_ancestry1['EUR']))
    print(len(my_dict_ancestry1['EAS']))
    print(len(my_dict_ancestry1['AMR']))
    print(len(my_dict_ancestry1['SAS']))
    '''



    #def reverse_mutation_dictionary(first2pairs_RE, first2pairs_SV,first2pairs_CN):



    '''
    We create a new dictionary where mutation is the key and ords are the values. This way we can find which samples
    have a certain mutation.
    '''

    # mutation_dictionary()

    #Initializing new dictionaries to be used

    new_dict_RE = {}
    new_dict_SV = {}
    new_dict_CN = {}


    #Using a for loop to go through the dictionary first2parirs_RE that we initialized in mutation_dictionary function.
    # Another for loop to go inside the values since a ORD (key) can have multiple mutations( value). This 2nd for loop
    # iterates through this value. If the mutation isn't present assigns it as the key. If it's present then it assigns
    # the adds the original key(ORD number) as the value.

    for key, value in first2pairs_RE.items():
        for val in value:
            if val not in new_dict_RE:
                new_dict_RE[val] = [key]
            else:
                new_dict_RE[val].append(key)

    #Printing the length due to a long list. Output is something like {'ALK:EML4': ['ORD-05028', 'ORD-0667']} when len is removed
    #print(len(new_dict_RE))



    #Very similar to the above code. Iterates through the short variant dictionary and switches the key:value pair so that
    # SV mutation is the key and the samples that have that mutation are the values.

    for key, value in first2pairs_SV.items():
        for val in value:
            if val not in new_dict_SV:
                new_dict_SV[val] = [key]
            else:
                new_dict_SV[val].append(key)

    #print(len(new_dict_SV))


    #Same as above. Switches the CNA mutation with the ORD number as key:value pairs

    for key, value in first2pairs_CN.items():
        for val in value:
            if val not in new_dict_CN:
                new_dict_CN[val] = [key]
            else:
                new_dict_CN[val].append(key)

    #print(len(new_dict_CN))



    #def sorted_mutations_top10(new_dict_RE, new_dict_SV, new_dict_CN):

    # Creating 3 new array to store the top 10 values for highest number of mutations
    top_mutation_SV = []
    top_mutation_CN = []
    top_mutation_RE = []


    # Sorts the new_dict_RE created in reverse_mutation function in reverse order. It then takes the top 10 items in the
    # sorted dictionary. The key=lambda x:len(x[1]) specifies that the sorting should happen based on the value for
    # each key-value pair. In this case the number of key with the most value

    # We then assign the key and length of the values of the top 10 to a new variable called mut_counting_Re and print it.
    # Every key is then added to a list called top_mutation_RE

    for key, value in sorted(new_dict_RE.items(), key=lambda x: len(x[1]), reverse=True)[:10]:
        mut_counting_RE = (f"{key}: {len(value)}")
        print(mut_counting_RE)
        top_mutation_RE.append(key)

    print("")
    print("--------------------------------------------")
    print("")

    # Sorts the new_dict_SV created in reverse_mutation function in reverse order. It then takes the top 10 items in the
    # sorted dictionary. The key=lambda x:len(x[1]) specifies that the sorting should happen based on the value for
    # each key-value pair. In this case the number of key with the most value

    # We then assign the key and length of the values of the top 10 to a new variable called mut_counting and prints it.
    # Every key is then added to a list called top_mutation_SV

    for key, value in sorted(new_dict_SV.items(), key=lambda x: len(x[1]), reverse=True)[:10]:

        mut_counting = (f"{key}: {len(value)}")
        print(mut_counting)
        top_mutation_SV.append(key)

    '''
            For giving only the output of EGFR in the top 10 list.
                if 'EGFR' in key:
            mut_counting = (f"{key}: {len(value)}")
            print(mut_counting)
            top_mutation_SV.append(key)
        
    
    #For printing all of the mutations that has EGFR
    
    app_keys = [key for key in new_dict_SV.keys() if "EGFR" in key]
    
    print(len(app_keys))
    
    '''

    print("")
    print("--------------------------------------------")
    print("")

    # Same as above but kidna doen'st work because there are same mutations that have different amplifcation values
    # We do get a top 10 list of mutations. Losses are all at 0.0 so their data is accurate.

    for key, value in sorted(new_dict_CN.items(), key=lambda x: len(x[1]), reverse=True)[:10]:
        if 'EGFR' in key:
            mut_counting1 = (f"{key}: {len(value)}")
            print(mut_counting1)
            top_mutation_CN.append(key)

    print(top_mutation_CN)



def get_sampleID_for_rearrangement(mutation, new_dict_RE):
    '''
    This function takes a mutation as argument and gives the sample that have the given mutation.
    :param mutation:
    :return:
    '''

    key_to_get = mutation

    #print(new_dict.keys())

    if key_to_get in new_dict_RE:
        print(len(new_dict_RE[key_to_get]))
        print(new_dict_RE[key_to_get])
    else:
        print(key_to_get, "mutation not found in Rearrangement.")
        print("Input mutation should be in the format of 'IGH:BCL2'")

    print("-----------------------------------------")




def get_sampleID_for_CNA(mutation, new_dict_CN):
    '''
    This function takes a mutation as argument and gives the samples that have the given mutation.
    :param mutation: The mutation that you want to check
    :return:
    '''

    key_to_get = mutation

    #Find the ORD for a given CNA amp/loss. If the key to get is in the dictionary
    if key_to_get in new_dict_CN:
        print(len(new_dict_CN[key_to_get]))
       # print(new_dict_CN[key_to_get1])
    else:
        print(key_to_get, "amp/loss not found in CNA")
        print("Input mutation should be in the format of 'CDKN2A:deletion:0.0'")





def get_sampleID_for_CNA_fixed(mutation, new_dict_CN):
    '''
    Does similar stuff to the function but is slightly more powerful. It fixes the "drawback" of the above code. When
     we aren't looking particularly for a specific mutation, this code can be used. E.g when the mutation for the above
     fucntions get_sampleID_for_CNA is just "CDKN2" instead of "CDKN2A:deletion:0.0", the function above gives us empty list
     but this function gives all of the samples that contain the CDKN2A mutation + CDKN2B mutation + CDKN2C mutation etc.

    :param mutation: The mutation to be searched.
    :return:
    '''

    print("-----------------------------------------")
    mut_arr = []
    for key in new_dict_CN:
        if 'CDKN2B' in key:
            mut_arr.extend(new_dict_CN[key])

    print(len(mut_arr))





def mutation_per_ancestry_RE(mutation, first2pairs_RE, ancestry_list, my_dict_ancestry1):
    '''
    This function returns that percentage breakdown of the mutation per ancestry
    :return:
    '''
    #Prints the rearrangement mutation percentage per ancestry

   # mutation_count = 'RET:KIF5B'
    #Assisgn mutaiton_count variable to the mutation provided. Initialized other variables like count and count_total to 0

    mutation_count = mutation
    print("Mutation = ", mutation_count)
    print(" ")
    count = 0
    count_total = 0


    '''
    Goes through the ancestry_list. Each element of that list is assigned the variable ancestry while going through the
    for loop. It then has another for loop where my_dict_ancestry is searched for ancestry from the first line. The ancestry
    is compared with the key of my_dict_ancestry. We then get a ORD number and then assign the samples(values) to from the
    dictionary first2pairs_RE. Then if the mutation_count is present in the values, the variable count is increased by 1
    and the percentage is counted. 
    
    Example: AFR is obtained from ancestry_list, my_dict_ancestry1 is searched for keys that are AFR. If AFR is present
     values is assigned to mutations that have AFR ancestry. If our mutation to be counted is in value then variable
    count is increased by 1 and percentage is calculated. The loop is then repeated for EUR, EAS etc. Finally the results
    are displayed out.
    '''

    for ancestry in ancestry_list:
        for key in my_dict_ancestry1[ancestry]:
            if key in first2pairs_RE:
                values = first2pairs_RE[key]
               # print(values)
                if mutation_count in values:
                    count += values.count(mutation_count)
                    perct = round((count*100)/ len(my_dict_ancestry1[ancestry]),2)

        print(ancestry, "has a total of [",count,"] mutations at a percentage of :", perct)
        count_total += count
        count = 0

    print("Total mutations = ", count_total)
    print(" ")





def mutation_per_ancestry_SV(mutation, ancestry_list, my_dict_ancestry1, first2pairs_SV, ):
    '''
    Prints SV percentage per ancestry
    :param mutation:
    :return:
    '''

    #mutation_count = 'EGFR:E746_A750del'
    # Assisgn mutaiton_count variable to the mutation provided. Initialized other variables like count and count_total to 0

    mutation_count = mutation
    print("Mutation = ", mutation_count)
    print(" ")
    count = 0
    count_total = 0

    '''
       Goes through the ancestry_list. Each element of that list is assigned the variable ancestry while going through the
    for loop. It then has another for loop where my_dict_ancestry is searched for ancestry from the first line. The ancestry
    is compared with the key of my_dict_ancestry. We then get a ORD number and then assign the samples(values) to from the
    dictionary first2pairs_SV. Then if the mutation_count is present in the values, the variable count is increased by 1
    and the percentage is counted. 
    
    Example: AFR is obtained from ancestry_list, my_dict_ancestry1 is searched for keys that are AFR. If AFR is present
     values is assigned to mutations that have AFR ancestry. If our mutation to be counted is in value then variable
    count is increased by 1 and percentage is calculated. The loop is then repeated for EUR, EAS etc. Finally the results
    are displayed out.
    '''

    for ancestry in ancestry_list:
        for key in my_dict_ancestry1[ancestry]:
            if key in first2pairs_SV:
                values = first2pairs_SV[key]
                if mutation_count in values:
                    count += values.count(mutation_count)
                    perct = round((count*100)/ len(my_dict_ancestry1[ancestry]),2)

        print(ancestry, "has a total of [",count,"] mutations at a percentage of :", perct)
        count_total += count
        count = 0
    print(" ")
    print("Total mutations = ", count_total)



def mutation_per_array(top_mutation_SV, ancestry_list, my_dict_ancestry1, first2pairs_SV):
    '''
    #Prints out the SV percentage per the array generated in  get_sampleID_for_CNA_fixed(
    :param mutd_array:
    :return:
    '''

    count = 0
    count_total = 0

    for mut in top_mutation_SV:
        mutation_count = mut
        #print("Mutation = ", mutation_count)
        print(" ")


        for ancestry in ancestry_list:
            for key in my_dict_ancestry1[ancestry]:
                print(key)
                if key in first2pairs_SV:
                    values = first2pairs_SV[key]
                    if mutation_count in values:
                        count += values.count(mutation_count)
                        perct = round((count*100)/ len(my_dict_ancestry1[ancestry]),2)



            print(mut, ancestry, "has a total of [",count,"] mutations at a percentage of :", perct)
            count_total += count
            count = 0
           # print(" ")
        print("Total mutations = ", count_total)
        count_total = 0



mutations_array = ['EGFR:G719A', 'EGFR:R108K', 'EGFR:E746_A750del', 'EGFR:L858R', 'EGFR:V742I', 'EGFR:L692V', 'EGFR:L833V', 'EGFR:T790M', 'EGFR:V774_C775insHV', 'EGFR:V774M', 'EGFR:H773R', 'EGFR:C797S', 'EGFR:L747_S752del', 'EGFR:D770_N771insSVD', 'EGFR:L747_A750>P', 'EGFR:V769_D770insASV', 'EGFR:S768I', 'EGFR:V769L', 'EGFR:H773_V774insTH', 'EGFR:H773_V774insNPH', 'EGFR:E746_S752>V', 'EGFR:H773_V774insH', 'EGFR:G724S', 'EGFR:L861Q', 'EGFR:E746_A750>IP', 'EGFR:L747_T751>P', 'EGFR:L747_P753>S', 'EGFR:E709A', 'EGFR:G719C', 'EGFR:E709G', 'EGFR:E709K', 'EGFR:G719D', 'EGFR:L747_T751del', 'EGFR:E709_T710>D', 'EGFR:N771_P772insN', 'EGFR:E746_T751>I', 'EGFR:V834L', 'EGFR:N771_P772insT', 'EGFR:D761Y', 'EGFR:L833F', 'EGFR:V765M', 'EGFR:E746_T751>V', 'EGFR:K754E', 'EGFR:N771>GF', 'EGFR:K714E', 'EGFR:N771_P772insH', 'EGFR:S768T', 'EGFR:L718V', 'EGFR:L718Q', 'EGFR:D770_N771insG', 'EGFR:T751_I759>S', 'EGFR:D770>GY', 'EGFR:E746_P753>VS', 'EGFR:G721C', 'EGFR:L747P', 'EGFR:L747_P753>Q', 'EGFR:F712L', 'EGFR:V292L', 'EGFR:S720F', 'EGFR:V843I', 'EGFR:E746_P753>IS', 'EGFR:L704F', 'EGFR:S768_V769>IL', 'EGFR:L747_S752>SRD', 'EGFR:K745_E746insIPVAIK', 'EGFR:P589L', 'EGFR:C797N', 'EGFR:V774L', 'EGFR:G719S', 'EGFR:D770_N771insGD', 'EGFR:E114K', 'EGFR:H773_V774insDPHPH', 'EGFR:D770H', 'EGFR:G719R', 'EGFR:E746_P753>VNS', 'EGFR:E746_T751>A', 'EGFR:H773L', 'EGFR:T854A', 'EGFR:Y998fs*4', 'EGFR:R831H', 'EGFR:D379Y', 'EGFR:A289V', 'EGFR:C775_R776insPHVC', 'EGFR:K1179_E1180>N*', 'EGFR:L833_V834>FL', 'EGFR:L747_K754>SPE', 'EGFR:G779F', 'EGFR:L747_S752>Q', 'EGFR:G1185V', 'EGFR:H870Q', 'EGFR:H773>PNPY', 'EGFR:T751_I759>N', 'EGFR:D770_N771insT', 'EGFR:E709V', 'EGFR:N771_P772insHN', 'EGFR:A289T', 'EGFR:R776H', 'EGFR:L861R', 'EGFR:E746_A750>AP', 'EGFR:D855N', 'EGFR:H835L', 'EGFR:L62R', 'EGFR:S752_I759del', 'EGFR:L703V', 'EGFR:E746_S752>I', 'EGFR:N700S', 'EGFR:R748K', 'EGFR:E709Q', 'EGFR:H773_V774>LM', 'EGFR:D770_N771insSTH', 'EGFR:R222C', 'EGFR:Q105H', 'EGFR:S229C', 'EGFR:G729R', 'EGFR:E746_T751>VP', 'EGFR:G873E', 'EGFR:E746_T751>VA', 'EGFR:G735S', 'EGFR:V786M', 'EGFR:N771_P772insVDN', 'EGFR:L747_A755>SQQG', 'EGFR:E412K', 'EGFR:splice site 2947-237_3114+116del521', 'EGFR:L747_S752>QH', 'EGFR:V769M', 'EGFR:A763_Y764insFQEA', 'EGFR:V689L', 'EGFR:E1196*', 'EGFR:H773_V774insPH', 'EGFR:R836H', 'EGFR:A871G', 'EGFR:R108G', 'EGFR:E746_A750>QP', 'EGFR:I759T', 'EGFR:E330K', 'EGFR:I759M', 'EGFR:I732T', 'EGFR:S645C', 'EGFR:N771_P772insC', 'EGFR:N771_P772insV', 'EGFR:I491V', 'EGFR:D770_N771insY', 'EGFR:H870R', 'EGFR:T725M', 'EGFR:V769_D770insGTL', 'EGFR:E746_T751>L', 'EGFR:P373Q', 'EGFR:E758_I759>DF', 'EGFR:D770_N771insGL', 'EGFR:T263P', 'EGFR:A871E', 'EGFR:D379E', 'EGFR:A743T', 'EGFR:H773>YPNPY', 'EGFR:L718M', 'EGFR:P772>HT', 'EGFR:L747_T751>Q', 'EGFR:K745_E746insLPVAIK', 'EGFR:G598V', 'EGFR:S768_V769insLDS', 'EGFR:S492R', 'EGFR:A289D', 'EGFR:D770>NNPH', 'EGFR:H773_V774insAH', 'EGFR:G796S', 'EGFR:S752F', 'EGFR:E545Q', 'EGFR:V769_D770insGSV', 'EGFR:E711K', 'EGFR:S306L', 'EGFR:V769_D770insGT', 'EGFR:E709_T710>A', 'EGFR:P772_H773insQPNP', 'EGFR:T790S', 'EGFR:K467E', 'EGFR:S768R', 'EGFR:K745_E746insVPVAIK', 'EGFR:H773Q', 'EGFR:L844P', 'EGFR:L838V', 'EGFR:E746_T751>IP', 'EGFR:R324L', 'EGFR:R836C', 'EGFR:L792H', 'EGFR:G796C', 'EGFR:H893Y', 'EGFR:E746_T751>AP', 'EGFR:N771>GY', 'EGFR:D770>ASVG', 'EGFR:splice site 2947-1G>C', 'EGFR:L747_P753>A', 'EGFR:R414S', 'EGFR:N771_P772>SVDNA', 'EGFR:I744M', 'EGFR:A763V', 'EGFR:D46N', 'EGFR:H805L', 'EGFR:H773>YNPY', 'EGFR:S752Y', 'EGFR:L747_T751>S', 'EGFR:V651M', 'EGFR:P772_H773insGHP', 'EGFR:A767_S768insVA', 'EGFR:P772_H773insDNP', 'EGFR:Q787R', 'EGFR:R776C', 'EGFR:splice site 3272-14_3279del22', 'EGFR:P753_I759>A', 'EGFR:L747_T751>N', 'EGFR:L747S', 'EGFR:E746_A750>ANKE', 'EGFR:K745E', 'EGFR:P772_H773insTHP', 'EGFR:P772T', 'EGFR:P596L', 'EGFR:K757R', 'EGFR:E709_T710>G', 'EGFR:K714N', 'EGFR:K860I', 'EGFR:V769_D770insSIV', 'EGFR:Q1173*', 'EGFR:E746_S752>L', 'EGFR:I759F', 'EGFR:P386L', 'EGFR:E758D', 'EGFR:L747_A755>SKD', 'EGFR:Q1113*', 'EGFR:D770A', 'EGFR:P772_H773insQ', 'EGFR:E734Q', 'EGFR:R705G', 'EGFR:N771>FH', 'EGFR:L747_E749del', 'EGFR:E872K', 'EGFR:T854S', 'EGFR:L858M', 'EGFR:L844V', 'EGFR:L747_K754>A', 'EGFR:D1168N', 'EGFR:T415M', 'EGFR:L387V', 'EGFR:D770>CVE', 'EGFR:E746_S752>T', 'EGFR:E709_T710>N', 'EGFR:A289P', 'EGFR:P772R', 'EGFR:Q1174*', 'EGFR:splice site 3114+1G>T', 'EGFR:H773_V774insGNPH', 'EGFR:E746_A750>L', 'EGFR:I759L', 'EGFR:T354M', 'EGFR:E1196fs*3', 'EGFR:E746_T751>VQ', 'EGFR:H773_V774insPHPH', 'EGFR:S768L', 'EGFR:E746_S752>A', 'EGFR:E865K', 'EGFR:N771>GT', 'EGFR:D770_N771insS', 'EGFR:E746_L747>APG', 'EGFR:E746_K754>ISE', 'EGFR:S752C', 'EGFR:E746_A750>KP', 'EGFR:N771>KT', 'EGFR:D761N', 'EGFR:T785I', 'EGFR:D314N', 'EGFR:L747_A755>SRD', 'EGFR:M766T', 'EGFR:E602Q', 'EGFR:P753L', 'EGFR:T751S', 'EGFR:S768_V769insMDS', 'EGFR:V769_D770insSSV', 'EGFR:L747_K754>QR', 'EGFR:E746_S752>QQ', 'EGFR:L747_P753>PKA', 'EGFR:D1012H', 'EGFR:H47Y', 'EGFR:G1189fs*10', 'EGFR:E746_A750>NP', 'EGFR:E709D', 'EGFR:G465V', 'EGFR:V769_D770insGVV', 'EGFR:D770_N771insD', 'EGFR:E1180*', 'EGFR:D770_N771>GVV', 'EGFR:L747_P753>QH', 'EGFR:L858W', 'EGFR:Q1143*', 'EGFR:A289I', 'EGFR:A289Y', 'EGFR:S1028fs*8', 'EGFR:L858Q', 'EGFR:V774A', 'EGFR:L747_A755>SKS', 'EGFR:E746D', 'EGFR:A750_I759>SN', 'EGFR:V689M', 'EGFR:splice site 3115-2A>T', 'EGFR:I744V', 'EGFR:C636Y', 'EGFR:D770_N771insGF', 'EGFR:N1053fs*25', 'EGFR:H773_V774insGHPH', 'EGFR:splice site 3115-2A>G', 'EGFR:A750_I759>PD', 'EGFR:P772_H773insA', 'EGFR:S768_V769>PL', 'EGFR:T751_A755del', 'EGFR:L747_P753>QQ', 'EGFR:A743_I744insVKIPVA', 'EGFR:L747V', 'EGFR:E746_T751>Q', 'EGFR:G719V', 'EGFR:S752_A755del', 'EGFR:Y727C', 'EGFR:L747_S752>ANN', 'EGFR:V769>LRD', 'EGFR:A763_Y764insFFEA', 'EGFR:splice site 3272-1G>C', 'EGFR:K1182*', 'EGFR:H850N', 'EGFR:splice site 3272-1G>A', 'EGFR:D770_N771insGN', 'EGFR:E746_T751>IA', 'EGFR:A289S', 'EGFR:E1193*', 'EGFR:D770N', 'EGFR:E204K', 'EGFR:D770>STLS', 'EGFR:A767V', 'EGFR:H773_V774insGTDNPH', 'EGFR:H773>TY', 'EGFR:K1061_E1062>M*', 'EGFR:E746del', 'EGFR:L747_K754>SQ', 'EGFR:A750_I759>QN', 'EGFR:A289L', 'EGFR:G1209fs*36', 'EGFR:R832C', 'EGFR:I740V', 'EGFR:D770>ES', 'EGFR:S768>II', 'EGFR:P772_H773insANP', 'EGFR:N771>KG', 'EGFR:N771>HPH', 'EGFR:N771_P772insG', 'EGFR:I475V', 'EGFR:P772_H773insGCP', 'EGFR:V802I', 'EGFR:P373L', 'EGFR:L747_A755>TRD', 'EGFR:A864V', 'EGFR:D770_N771insGT', 'EGFR:L703F', 'EGFR:P386F', 'EGFR:H773_V774insY', 'EGFR:I491L', 'EGFR:N771_P772insDN', 'EGFR:P772_H773insGNP', 'EGFR:P772_H773insNV', 'EGFR:splice site 425-1G>A', 'EGFR:N771>PH', 'EGFR:P772_H773insNPHP', 'EGFR:P772L', 'EGFR:M766I', 'EGFR:H773_V774insNPHPH', 'EGFR:Q432K', 'EGFR:splice site 3114+1G>A', 'EGFR:G863D', 'EGFR:L747_P753>QA', 'EGFR:V843L', 'EGFR:splice site 559+1G>A', 'EGFR:Q1126*', 'EGFR:splice site 3272-1G>T', 'EGFR:splice site 3115-1G>A', 'EGFR:D770G', 'EGFR:G863S', 'EGFR:P772>HR', 'EGFR:E746_S752>D', 'EGFR:P772_H773insNP', 'EGFR:D770_N771insP', 'EGFR:A750_I759>SS', 'EGFR:E746_T751>AA', 'EGFR:S1030*', 'EGFR:T415A', 'EGFR:D379N', 'EGFR:L1034fs*3', 'EGFR:E749K', 'EGFR:E746_S752>IPVA', 'EGFR:V769_D770insEASV', 'EGFR:P753S', 'EGFR:E749_E758>DP', 'EGFR:D770_N771insGTD', 'EGFR:L747_P753>N', 'EGFR:R414K', 'EGFR:A750_I759>PT', 'EGFR:A767D', 'EGFR:Y764_V765insHH', 'EGFR:E746V', 'EGFR:D837E', 'EGFR:H773_V774insHPH', 'EGFR:R252C', 'EGFR:L747_T751>AR', 'EGFR:D256A', 'EGFR:W731L', 'EGFR:S784F', 'EGFR:L747F', 'EGFR:H773_V774insSH', 'EGFR:G63R', 'EGFR:V774_C775insHNPHV', 'EGFR:A871P', 'EGFR:N771_P772>SEDNS', 'EGFR:N771_P772insGTDN', 'EGFR:E1015*', 'EGFR:A750_E758>P', 'EGFR:V323I', 'EGFR:E931G', 'EGFR:S151fs*17', 'EGFR:Y813F', 'EGFR:D770_N771>GSVDN', 'EGFR:R377S', 'EGFR:H773_V774>NPNPYV', 'EGFR:746_751ELREAT>L', 'EGFR:V769_H773>LDNPNPH', 'EGFR:H773_V774insNPY', 'EGFR:H773_V774insVH', 'EGFR:L747_A755>SKG', 'EGFR:D770>EG', 'EGFR:L747_A755>SNS', 'EGFR:M793_P794>ITQLMP', 'EGFR:H1124fs*44', 'EGFR:D770_N771insH', 'EGFR:R521fs*18', 'EGFR:E746_P753>LS', 'EGFR:N771_P772>SHP', 'EGFR:L747_A755>NFNS', 'EGFR:N771_P772>SPHP', 'EGFR:N1184fs*14', 'EGFR:P1019L', 'EGFR:E1062fs*16', 'EGFR:V769_D770insSLRD', 'EGFR:K745_A750>KIP', 'EGFR:H773_V774insPY', 'EGFR:P772_H773insGDP', 'EGFR:W410C', 'EGFR:L858_A859>RS', 'EGFR:P992fs*10', 'EGFR:D770_N771insAVD', 'EGFR:E746_R748>IP', 'EGFR:N771_H773>HHPH', 'EGFR:D984fs*18', 'EGFR:E746_S752>C', 'EGFR:L747_P753>QIA', 'EGFR:A750V', 'EGFR:D770_N771>GVMDN', 'EGFR:Q1159*', 'EGFR:E746_A750>RP', 'EGFR:L747_A750>S', 'EGFR:A767_S768insTLA', 'EGFR:V769_D770insGTV', 'EGFR:D770_N771>EDN', 'EGFR:N771_P772>HHP', 'EGFR:V769_D770insSASV', 'EGFR:A871S', 'EGFR:N771_P772insNN', 'EGFR:E746_T751>FPS', 'EGFR:splice site 559+2T>G', 'EGFR:E746_A750>VP', 'EGFR:D770_N771insSMD', 'EGFR:N771_P772>SVDNP', 'EGFR:D770_N771>ASVDS', 'EGFR:E746_T751>S', 'EGFR:P772_H773insQSPNP', 'EGFR:E746_T751>LS', 'EGFR:L747_K754>G', 'EGFR:S1071*', 'EGFR:T751_I758>D', 'EGFR:E749Q', 'EGFR:K757M', 'EGFR:I491T', 'EGFR:P772_H773insPNP', 'EGFR:D770>GTH', 'EGFR:E746K', 'EGFR:E1079*', 'EGFR:H805R', 'EGFR:splice site 3238_3271+5del39', 'EGFR:T790L', 'EGFR:F968L', 'EGFR:E746_L747>IP', 'EGFR:A750_L760>GDI', 'EGFR:L858S', 'EGFR:splice_site_3162+1_3162+1delG', 'EGFR:T751_L760>NL', 'EGFR:L747_N756>SPT', 'EGFR:C620W', 'EGFR:Q1164*', 'EGFR:splice site 2947-1G>A', 'EGFR:T415K', 'EGFR:E746A', 'EGFR:K708N', 'EGFR:S752_I759>QG', 'EGFR:L747_P753>QK', 'EGFR:E746_T751>VN', 'EGFR:T790I', 'EGFR:Y1125*', 'EGFR:H773_V774insGH', 'EGFR:K757E', 'EGFR:C797G', 'EGFR:I744_A750>VKIP', 'EGFR:L747_T751>A', 'EGFR:E746_A755>EAN', 'EGFR:D770_N771>RDG', 'EGFR:T751I', 'EGFR:A750_I759>GD', 'EGFR:E749_A755>DS', 'EGFR:R533M', 'EGFR:E746_A750>T', 'EGFR:E746_S752>EQH', 'EGFR:A750_E758>Q', 'EGFR:T790A', 'EGFR:T710_F712>P', 'EGFR:A1102fs*65', 'EGFR:A750_I759>AN', 'EGFR:L833W', 'EGFR:V769_N771>LDNRGH', 'EGFR:Y998*', 'EGFR:splice site 3163-1delG', 'EGFR:E746_A750>LP', 'EGFR:H773>QW', 'EGFR:W410L', 'EGFR:V769_D770insCV', 'EGFR:H433L', 'EGFR:A763_Y764insLQEA', 'EGFR:V774_C775insDNPHV', 'EGFR:L747_A755>QHLRS', 'EGFR:E746_T751>QP', 'EGFR:D800N', 'EGFR:P772S', 'EGFR:G721A', 'EGFR:I744_E749>MKI']
print(len(mutations_array))


def mutation_per_custom_array(top_mutation_SV, ancestry_list, first2pairs_SV, my_dict_ancestry1):

    count = 0
    af_count = eu_count = nan_count = eas_count = amr_count = sas_count = 0
    count_total = 0



    top_mutation_EGFR = ['EGFR:L858R','EGFR:E746_A750del','EGFR:T790M', 'EGFR:G719A', 'EGFR:L747_P753>S']

    for mut in top_mutation_SV:
        mutation_to_count = mut


        for ancestry in ancestry_list:
            for key in my_dict_ancestry1[ancestry]:
                if key in first2pairs_SV:
                    values = first2pairs_SV[key]
                    if mutation_to_count in values:
                        count += values.count(mutation_to_count)
                        perct = round((count*100)/ len(my_dict_ancestry1[ancestry]),2)
                        if ancestry == 'AFR':
                            af_count += 1
                        elif ancestry == 'EUR':
                            eu_count += 1
                        elif ancestry == 'SAS':
                            sas_count += 1
                        elif ancestry == 'EAS':
                            eas_count += 1
                        elif ancestry == 'AMR':
                            amr_count += 1
                        else:
                            nan_count += 1


            print(mut, ancestry, "has a total of [",count,"] mutations at a percentage of :", perct)
            count_total += count
            count = 0
        print(" ")
        print("Total mutations = ", count_total)
        count_total = 0

    print("af", af_count)
    print("eu", eu_count)
    print(" ")
    print("eas", eas_count)
    print("sas", sas_count)
    print("")
    print("amr", amr_count)
    print("nan", nan_count)


    #print("African EGFR mutation percentage:", (len(my_dict_ancestry1['AFR'])))

    #print("African EGFR mutation percentage:", (round((af_count*100)/ len(my_dict_ancestry1['AFR']),2)))

    #print(af_count/, eu_count, nan_count, eas_count, amr_count, sas_count)

    a = (nan_count*100)/3690
    print(a)

    # Looks like sample bias. A sample that has 2 different EGFR mutation gets counted twice.
    """
    
    print(" ")
    print("African EGFR mutation percentage:", (round((af_count*100)/ len(my_dict_ancestry1['AFR']),2)))
    print("European EGFR mutation percentage:", (round((eu_count*100)/ len(my_dict_ancestry1['EUR']),2)))
    print("East Asian EGFR mutation percentage:", (round((eas_count*100)/ len(my_dict_ancestry1['EAS']),2)))
    print("South Asian EGFR mutation percentage:", (round((sas_count*100)/ len(my_dict_ancestry1['SAS']),2)))
    print("AMR EGFR mutation percentage:", (round((amr_count*100)/ len(my_dict_ancestry1['AMR']),2)))
    print("nan EGFR mutation percentage:", (a))
    
    """




def mutation_gene_CNA(gene1, new_dict_CN,first2pairs_ancestry, my_dict_ancestry1):

    '''
    #Takes in a gene name and then searches the total amplification or deletions present. It counts the number of samples
    #that have the said gene amp or deletion and gives a percentage of the samples that have that gene's amp or deletion per ancestry.
    :param gene1:
    :return:
    '''

    #Initializing variables and assigninig it to 0 to start the program
    #Creates empty array and searches new_dict_CN and looks for a gene. If it contains the gene then adds the
    #sample ID in the new array.
    #

    af_count = eu_count = nan_count = eas_count = amr_count = sas_count = 0
    gene1_arr = []
    for key in new_dict_CN:
        if gene1 in key:
            gene1_arr.extend(new_dict_CN[key])


    print("Total number of mutations:", len(gene1_arr), "out of 87,846 samples")


    #Each element of gene1_arr is checked as the key and then the ancestry is searched. Based on the ancestry
    # the count of the ancestry is updated.

    for key in gene1_arr:

        if first2pairs_ancestry[key] == 'AFR':
            af_count += 1
        elif first2pairs_ancestry[key] == 'EUR':
            eu_count += 1
            #print("eu", key)
        elif first2pairs_ancestry[key] == 'SAS':
            sas_count += 1
        elif first2pairs_ancestry[key] == 'EAS':
            eas_count += 1
        elif first2pairs_ancestry[key] == 'AMR':
            amr_count += 1
        else:
            nan_count += 1

    print(af_count, eu_count, nan_count, eas_count , amr_count , sas_count)

    # The percentage of the mutation in a certain ancestry is calculated and assigned to a variable and printed.

    afc = round(af_count*100/8784,2)
    print(afc)
    euc = round(eu_count*100/70297,2)
    print(euc)
    eas = round(eas_count*100/3690,2)
    print(eas)
    sas = round(sas_count*100/565,2)
    print(sas)
    amr = round(amr_count*100/2407,2)
    print(amr)
    nan = round(nan_count*100/3690,2)
    print(nan)


    print("-------------------------------------------------")
    print("African",gene1, " mutation percentage:", afc)
    print("European",gene1, "mutation percentage:", (round((eu_count*100)/ len(my_dict_ancestry1['EUR']),2)))
    print("East Asian",gene1," mutation percentage:", (round((eas_count*100)/ len(my_dict_ancestry1['EAS']),2)))
    print("South Asian",gene1,"mutation percentage:", (round((sas_count*100)/ len(my_dict_ancestry1['SAS']),2)))
    print("AD Mixed American",gene1," mutation percentage:", (round((amr_count*100)/ len(my_dict_ancestry1['AMR']),2)))


def main():


   # setup1()
    #setup2()
    setup3()

    #Runs the mutation_dictionary function that starts the whole process of creating dictionaries


    gene1 = 'IGH:BCL2'
   # mutation_gene_CNA(gene1, new_dict_CN, first2pairs_ancestry, my_dict_ancestry1)
    #get_sampleID_for_rearrangement(gene1,new_dict_RE, first2pairs_RE)



if __name__ == '__main__':
    main()