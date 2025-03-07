#Cytogenetics ONLY tool for MDS v1.2.1

import streamlit as st
import pandas as pd
import numpy as np
import re


# CG STRING PROCESSING

#--idem processing

def process_idem(row):
    cytogenetics_value = str(row[cytovar])

    cytogenetics_value = cytogenetics_value.replace(" ", "") #whitespace elimination
    cytogenetics_value = cytogenetics_value.replace("{", "[" ) #bracket change
    cytogenetics_value = cytogenetics_value.replace("}", "]" ) #bracket change
    cytogenetics_value = cytogenetics_value.replace("]", "]/") #add / after bracet for it to detect segment
    cytogenetics_value = cytogenetics_value.replace("];", "]/") #semi-colon change to fit split below
    cytogenetics_value = cytogenetics_value.replace("],", "]/") #comma change to fit split below
    cytogenetics_value = cytogenetics_value.replace("]//", "]/") #double / elimination fit split below

  
    
    if cytogenetics_value.endswith("/"):
        cytogenetics_value = cytogenetics_value.rstrip("/")    

    cytogenetics_value = "".join([char.lower() if char not in ['X', 'Y'] else char for char in cytogenetics_value])

    segments = cytogenetics_value.split('/')
    
    first_segment_abnl = re.sub(r'\[\d+\]', '', ",".join(segments[0].split(',')[1:]))
    
    #cum abnl, first segment
    cumulative_abnl = first_segment_abnl

    for i in range(1, len(segments)): #skip first segment 
        if 'idem' in segments[i]:
            #get the chromosomes
            chromosome_count = segments[i].split(',')[0]
            #replace idem with cum abnl
            segments[i] = chromosome_count + "," + cumulative_abnl + "," + ",".join(segments[i].split(',')[1:]).replace('idem,', '')

        #update the cumm abnormalities without adding duplicates
        current_segment_abnl = re.sub(r'\[\d+\]', '', ",".join(segments[i].split(',')[1:]))
        for abnl in current_segment_abnl.split(','):
            if abnl and abnl not in cumulative_abnl:
                cumulative_abnl += "," + abnl
    return '/'.join(segments)



def count_abn(row):
    cytogenetics_value = str(row["processed_cg"])
    
    segments = cytogenetics_value.split('/')

    max_abn = 0

    for segment in segments:
        if segment.endswith('[1]'): #ignoring the [1] segments
            continue

        abns = segment.split(',')[2:]

        seg_abn_count = len(abns)

        if seg_abn_count > max_abn:
            max_abn = seg_abn_count
    
    return max_abn


# VERY GOOD RISK ABNORMALITIES# %%
#"""-Y(minusy)"""

def minusy(row):
    cytogenetics_value = row["processed_cg"]

    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue

        if isinstance(segment, str):
            pattern = r'-Y|- Y'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0


# """del11q(delelevenq)"""

def elevenq(row):
    cytogenetics_value = row["processed_cg"]

    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'): 
            continue

        if isinstance(segment, str):
            pattern = r'del\(11\)\(q'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0

# GOOD RISK ABNORMALITIES

# """"del5q"""

def delfiveq(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue

        if isinstance(segment, str):
            pattern = r'del\(5\)\(q'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0


# """"del12p"""

def deltwelvep(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue

        if isinstance(segment, str):
            pattern = r'del\(12\)\(p'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0

#
# """del20q"""

def deltwentyq(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue
        
        if isinstance(segment, str):
            pattern = r'del\(20\)\(q'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0


# INTERMEDIATE RISK ABNORMALITIES

# """del7q"""

def delsevenq(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue
        if isinstance(segment, str):
            pattern = r'del\(7\)\(q'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0


# """"plus8"""

def pluseight(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue
        if isinstance(segment, str):
            pattern = r'\+8'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0


# """plus19"""

def plusnineteen(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue
        if isinstance(segment, str):
            pattern = r'\+19'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0


# """i17q"""

def iseventeenq(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue
        if isinstance(segment, str):
            pattern = r'i\(17\)\(q'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0



# POOR RISK ABNORMALITIES

# """minus7"""

def minusseven(row):
    cytogenetics_value = row["processed_cg"]

    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue

        if isinstance(segment, str):
            pattern = r'\-7'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0



# """inv3/t3q/del3q (inv_del_t_3q)"""

def chr3abn(row):
    cytogenetics_value = row["processed_cg"]

    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue

        if isinstance(cytogenetics_value, str):
            pattern = r'inv\(3|t\(3|del\(3\)\(q'
            match = re.search(pattern, cytogenetics_value)
            if match:
                return 1
    return 0



# """"DEL 17p AND -17""" #not in CG risk but useful for IPSS-M calculation

def delseventeen(row):
    cytogenetics_value = row["processed_cg"]

    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue
        if isinstance(cytogenetics_value, str):
            pattern = r'del\(17\)\(p|-17|- 17'
            match = re.search(pattern, cytogenetics_value)
            if match:
                return 1
    return 0


# """Diploid or 46,XX or 46,XY """

othervars = ["del77q", "del55q", "del1717p", "elevenq", "minusy"]

def diploid(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue
    
        pattern = r'(46,XX|46,XY)\[(\d+)\]'
        match = re.search(pattern, segment)
        if match:
            metaphase_count = int(match.group(2))
            if metaphase_count >= 10:
                return 1  # High confidence
            else:
                return 2 # Low confidence

    return 0  # Not diploid

def metaphase(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue

        pattern = r'\[(?:cp)?(\d+)\]'
        match = re.search(pattern, segment)

        if match:
            metaphase_count = int(match.group(1))
            return metaphase_count
                


# CG-RISK CALCULATION


def segments(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split("/")
    
    count = 0
    
    for segment in segments:
        if not segment.endswith('[1]'):
            count += 1
            
    return count


# """CALCULATE CG RISK"""

good_vars = ["del20q", "del12p", "del5q"]

intermediate_or_higher_vars = ["del7q", "minus7", "plus8", "plus19", "i17q", "inv_del_t_3q"]

higher_vars = ["minus7", "inv_del_t_3q"]

def cg_risk(row): #this is the right order of resolving the if statements to call cg_risk 10/2023
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split("/")

    for segment in segments:

        if segment.endswith("[1]"):
            continue

        #easy: 4 "Very Poor" if highly complex more than 3 abnormalities (see count_abn function)
        if row["abn_total"] >3:
            return 4
            
        #easy: "Very Good" if isolated -Y or isolated del11q
        if ((row["loss_of_y"] == 1) and (row[good_vars].any() ==0) and (row[intermediate_or_higher_vars].any() ==0)) or \
            ((row["del11q"] == 1) and (row[good_vars].any() ==0) and (row[intermediate_or_higher_vars].any() ==0)):
            return 0 
        
        #step 1 complexity, 
        #diploid must be 1  
        #del12p must be isolated + abn_total = 1, 
        #del20p must be isolated + abn_total = 1,  
        #del5q must be isolated, or del5q + other clone, except clones with intermediate or high-risk abnormalities
        if (row["diploid"] == 1) or \
            ((row["del12p"] == 1) and (row["abn_total"] == 1)) or \
            ((row["del20q"]) == 1 and (row["abn_total"] == 1)) or \
            ((row["del5q"] ==1) and (row["clone_total"] <=2) and (row[intermediate_or_higher_vars].any() == 0)):
            return 1
        
        #step 2 complexity, minus7 or del7 + another, -7 alone, abn3 any, or 3 abnormalities (complex=3)
        if (row["minus7"] == 1) or \
            (row["inv_del_t_3q"] == 1) or \
            ((row["minus7"] == 1) and (row["clone_total"] == 2)) or \
            ((row["del7q"] == 1) and (row["abn_total"] >= 2)) or \
            (row["abn_total"] == 3 ):
            return 3

        #basically almost anything else, isolated del7q, isolated +8, isolated +19, isolated i17q, or clone >=2 and abn total =2, or any other clone,
        #with diploid, del5q, intermediate or higher vars have to be 0
        if ((row["del7q"] == 1) and (row["abn_total"] == 1)) or\
            ((row["plus8"] == 1) and (row["abn_total"] == 1)) or \
            ((row["plus19"] == 1) and (row["abn_total"] == 1)) or \
            ((row["i17q"] == 1) and (row["abn_total"] == 1)) or \
            ((row["clone_total"] >= 2) and (row["abn_total"] <= 2) and (row["diploid"] == 0) and (row[higher_vars].any() == 0) and (row["del5q"]== 0)) or \
            ((row["clone_total"] == 2) and (row[higher_vars].any() == 0)):
            return 2
    
        
        if not any(char.isdigit() for char in segment):
            return "manual check required"
        
        else:
            return "manual check required"    
    
    return "manual check required"    


def check_numeric_and_sum(row):
    if row.apply(lambda x: isinstance(x, str)).any():
        return "manual check required"
    else:
        return row.sum()

        
        
def master_function(data):
    data["processed_cg"] = data.apply(process_idem, axis=1)    
    data["abn_total"] = data.apply(count_abn, axis=1)
    data["clone_total"] = data.apply(lambda row: segments(row), axis=1)
    data["loss_of_y"] = data.apply(lambda row: minusy(row) if row["abn_total"] == 1 else 0, axis=1)
    data["del11q"] = data.apply(lambda row: elevenq(row) if row["abn_total"] == 1 else 0, axis=1)
    data['del5q'] = data.apply(lambda row: delfiveq(row), axis=1)
    data['del12p'] = data.apply(lambda row: deltwelvep(row), axis=1)
    data['del20q'] = data.apply(lambda row: deltwentyq(row), axis=1)
    data["del7q"] = data.apply(lambda row: delsevenq(row), axis=1)
    data["plus8"] = data.apply(lambda row: pluseight(row), axis=1)
    data["plus19"] = data.apply(lambda row: plusnineteen(row), axis=1)
    data["i17q"] = data.apply(lambda row: iseventeenq(row), axis=1)
    data["minus7"] = data.apply(lambda row: minusseven(row), axis=1) 
    data["inv_del_t_3q"] = data.apply(lambda row: chr3abn(row), axis=1) #includes those with any amount of abnormalities 
    data["del17or17p"] = data.apply(lambda row: delseventeen(row), axis=1)
    data["diploid"] = data.apply(lambda row: diploid(row) if row['abn_total'] == 0 else 0, axis=1) #this function may need work
    data["metaphase"] = data.apply(lambda row: metaphase(row), axis =1)
    data["cg_risk_score"] = data.apply(lambda row: cg_risk(row), axis=1)

    cols = [cytovar, "metaphase", "loss_of_y", "del11q","del5q", "del12p", "del20q", "del7q", "plus8", "plus19", "i17q", "minus7", "inv_del_t_3q", "del17or17p", "diploid", "cg_risk_score"]
    return data[cols]


st.title("Cytogenetic risk calculator for myelodysplastic syndrome")

st.markdown("""Please enter the column names in your file without spaces below: 
            Please make sure your columns contain **only** numeric values.""")
cytovar = st.text_input("Enter the name of the cytogenetics column (better performance with ISCN 2020)")


u_file = st.file_uploader("Please upload a CSV file. Please convert your Excel file to CSV prior to uploading.", type =['csv'])

if st.button('Disclaimer'):
    st.write("This is a first relaease app. Please check the results. The cytogenetic risk caller performs better (>99%) with ISCN 2020 nomenclature. If the app is unsure it will call for a manual check for the specific row. Please email me with questions to the address below for support.")

st.markdown("""
#### Contact Information:
- **Created by:** Samuel Urrutia
- **Email:** surrutia@wustl.org
- All files are deleted upon closing the browser.
""")

if u_file is not None:
    if u_file.name.endswith('.csv'):
        data = pd.read_csv(u_file)
        st.write("Original Data")
        st.write(data)
     
    try:
            processed_data = master_function(data.copy())
            st.write("Processed Data")
            st.write(processed_data)

            proc_data = processed_data.to_csv(index=False)

            st.download_button(
                label="Download Processed Data",
                data=proc_data,
                file_name="processed_cytogenetic_data.csv",
                mime="text/csv"
            )
    except Exception as e:
        st.write(f"An error occurred: {e}")
else:
    st.write("Please upload a file.")

    
