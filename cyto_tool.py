#Urrutia Cytogenetics tool

import streamlit as st
import pandas as pd
import re

def process_idem(row):
    cytogenetics_value = str(row["Cyto at DX"])

    cytogenetics_value = cytogenetics_value.replace(" ", "")

    segments = cytogenetics_value.split('/')
    
    first_segment_abnl = re.sub(r'\[\d+\]', '', ",".join(segments[0].split(',')[1:]))
    
    #cum abnl, first segment
    cumulative_abnl = first_segment_abnl

    for i in range(1, len(segments)): #skip first segment 
        if 'idem' in segments[i]:
            #get the chromosomes
            chromosome_count = segments[i].split(',')[0]
            #replace i dem with cum abnl
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


#VERY GOOD RISK ABNORMALITIES

#-Y(minusy)

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

#del11q(delelevenq)

def elevenq(row):
    cytogenetics_value = row["processed_cg"]

    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):  #this terms makes them ignore those with [1] abnormality. 
            continue

        if isinstance(segment, str):
            pattern = r'del\(11\)\(q'
            match = re.search(pattern, segment)
            if match:
                return 1
    return 0

#GOOD RISK ABNORMALITIES

#del5q

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

#del12p

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



#del20q

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

#INTERMEDIATE RISK ABNORMALITIES

#del7q 

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


# %%
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


#plus19

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



#i17q

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



#POOR RISK ABNORMALITIES 

#minus7

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



#inv3/t3q/del3q (inv_del_t_3q)

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



#DEL 17p AND -17 #note in CG risk but useful for IPSS-M calculation later

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


#Diploid or 46,XX or 46,XY

othervars = ["del77q", "del55q", "del1717p", "elevenq", "minusy"]

def diploid(row):
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split('/')

    for segment in segments:
        if segment.endswith('[1]'):
            continue
        if isinstance(cytogenetics_value, str):
            pattern = r'46,XX\[|46,XY\[|Diploid CG'
            match = re.search(pattern, cytogenetics_value)
            if match:
                return 1
    return 0


#CG-RISK CALCULATION


def segments(row): #this function is not callad in master function
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split("/")
    return len(segments)


#CALCULATE CG RISK

intermediate_or_higher_vars = ["del7q", "minus7", "plus8", "plus19", "i17q", "inv_del_t_3q"]

def cg_risk(row): #this is the right order of resolving the if statements to call cg_risk 10/2023
    cytogenetics_value = row["processed_cg"]
    segments = cytogenetics_value.split("/")

    for segment in segments:

        if segment.endswith("[1]"):
            continue

        if row["abn_total"] >3:
            return 4
        
        if (row["minusy"] == 1) or (row["delelevenq"] == 1):
            return 0 
        
        if (row["diploid"] == 1) or ((row["del12p"] == 1) and (row[intermediate_or_higher_vars].any() == 0)) or (row["del20q"] == 1 and (row[intermediate_or_higher_vars].any() == 0)) \
        or ((row["del5q"] ==1) and (len(segments) <= 6) and (row["abn_total"]==1) and (row["diploid"] == 0) and (row[intermediate_or_higher_vars].any() == 0)):
            return 1
        
        elif (row["minus7"] == 1) or (row["inv_del_t_3q"] == 1) or ((row["minus7"] == 1) and (len(segments) == 2))\
                or ((row["del7q"] == 1) and (len(segments) > 2)) or (row["abn_total"] == 3 ):
            return 3
        
        elif (row["del7q"] == 1) or (row["plus8"] == 1) or (row["plus19"] == 1) or (row["i17q"] == 1) \
        or ((len(segments) >=2) and (row["abn_total"] <= 2) and (row["diploid"] == 0) and (row["del5q"]== 0)) or (row["abn_total"] == 2) or (row["abn_total"] == 1):
            return 2
        
        else:
            return "manual check required"        
        


def master_function(data):
    data["processed_cg"] = data.apply(process_idem, axis=1)
    data["abn_total"] = data.apply(count_abn, axis=1)
    data["minusy"] = data.apply(lambda row: minusy(row) if row["abn_total"] == 1 else 0, axis=1)
    data["delelevenq"] = data.apply(lambda row: elevenq(row) if row["abn_total"] == 1 else 0, axis=1)
    data['del5q'] = data.apply(lambda row: delfiveq(row) if row["abn_total"] == 1 else 0, axis=1)
    data['del12p'] = data.apply(lambda row: deltwelvep(row) if row["abn_total"] == 1 else 0, axis=1)
    data['del20q'] = data.apply(lambda row: deltwentyq(row) if row["abn_total"] == 1 else 0, axis=1)
    data["del7q"] = data.apply(lambda row: delsevenq(row) if row["abn_total"] == 1 else 0, axis=1)
    data["plus8"] = data.apply(lambda row: pluseight(row) if row['abn_total'] == 1 else 0, axis=1)
    data["plus19"] = data.apply(lambda row: plusnineteen(row) if row['abn_total'] == 1 else 0, axis=1)
    data["i17q"] = data.apply(lambda row: iseventeenq(row) if row['abn_total'] == 1 else 0, axis=1)
    data["minus7"] = data.apply(lambda row: minusseven(row) if row['abn_total'] <=2 else 0, axis=1) 
    data["inv_del_t_3q"] = data.apply(lambda row: chr3abn(row) if row["abn_total"] <=2 else 0, axis=1) #includes those with any amount of abnormalities 
    data["del1717p"] = data.apply(lambda row: delseventeen(row), axis=1)
    data["diploid"] = data.apply(lambda row: diploid(row) if row['abn_total'] == 0 else 0, axis=1) #this function mayu need work needs work
    data["cg_risk"] = data.apply(lambda row: cg_risk(row), axis=1)

    return data

#This is the streamlit webapp try
   
st.title("MDS CG Processor")

u_file = st.file_uploader("Please upload a CSV file (please convert your Excel file to CSV prior to uploading and name your CG column 'cytogenetics')", type =['csv'])

if u_file is not None:
    if u_file.name.endswith('.csv'):
        data = pd.read_csv(u_file)

    st.write("Original Data")
    st.write(data)

    if 'Cyto at DX' in data.columns:
        processed_data = master_function(data.copy())
        st.write("Processed Data")
        st.write(processed_data)
    else:
        st.error ("The file does not have a Cyto at DX column.")

