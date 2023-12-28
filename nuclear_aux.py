# -*- coding: utf-8 -*-
"""
Auxiliary function for the MESA nuclear network retriever.
process_line(line) processes an isotope 'line' and returns the element name and isotopes included.
read_file(data_dict,file_path,net_name) 
"""

import re
import os
import argparse

# define your command line arguments
parser = argparse.ArgumentParser(description='Auxiliary function for the MESA nuclear network retriever.\n \
# process_line(line): \n \
    processes an isotope \'line\' and returns the element name and isotopes included.\n \
# read_file(data_dict,file_path,net_name): \n \
    accepts (and returns) a dict (\'data_dict\') with element names as the keys and isotopes as values.\n \
    \'file path\' is used to indicate network file location, \n \
    \'network\' the requested network name.')


def process_line(line):
    parts = line.strip().split()
        
    # Check for "Y" format (i.e., 'neut', 'prot')
    if len(parts) == 1 and parts[0][-1].isalpha():
        return parts[0], []
    
    # Check for "X#" format (i.e., 'fe54', 'h2'
    if len(parts) == 1:
        return re.findall("[a-zA-Z]+", parts[0])[0], [int(re.findall(r'\d+', parts[0])[0])]
    
    # Check for "X #1 #2" format (i.e., 'fe 52 54', 'he 3 4')
    if len(parts) == 3 and parts[0].isalpha() and all(part.isdigit() for part in parts[1:]):
        x = parts[0]
        start, end = int(parts[1]), int(parts[2])
        return x, list(range(start, end + 1))
    
    # If the format is not recognized, return an empty list for the values
    return "", []

def read_file(data_dict,file_path,net_name):
    # append the .net at the end if missing
    if net_name[-4:]!='.net':
        net_name = net_name + '.net'
    # exclude the approx networks:
    if net_name=="approx19.net": #h1, he3, he4, c12, n14, o16, ne20, mg24, si28, s32, ar36, ca40, ti44, cr48, fe52, fe54, ni56, neut
        data_dict = data_dict | {'neut':[],'h':[1],'he':[3,4],'c':[12],'n':[14],'o':[16],'ne':[20],'mg':[24],'si':[28],'s':[32],'ar':[36],'ca':[40],'ti':[44],'cr':[48],'fe':[52,54],'ni':[56]}
        return data_dict
    elif net_name=="approx20.net":
        data_dict = data_dict | {'neut':[],'h':[1],'he':[3,4],'c':[12],'n':[14],'o':[16],'ne':[20],'mg':[24],'si':[28],'s':[32],'ar':[36],'ca':[40],'ti':[44],'cr':[48],'fe':[52,54,56],'ni':[56]}
        return data_dict
    elif net_name=="approx21.net":
        data_dict = data_dict | {'neut':[],'h':[1],'he':[3,4],'c':[12],'n':[14],'o':[16],'ne':[20],'mg':[24],'si':[28],'s':[32],'ar':[36],'ca':[40],'ti':[44],'cr':[48,56],'fe':[52,54,56],'ni':[56]}
        return data_dict
    # make sure the net exists:
    if not os.path.isfile(os.path.join(file_path,net_name)):
    #.. and if it doesn't, suggest a fix:
        if net_name[0:4]=='mesa':
            net_name_corrected = 'mesa_'+net_name[4:]
            if os.path.isfile(os.path.join(file_path,net_name_corrected)):
                print('Requested network {} is missing, I\'ve allowed myself to check for {} which was present'.format(net_name,net_name_corrected))
                net_name = net_name_corrected
                del net_name_corrected
            else:
                print('Requested network {} is missing'.format(net_name))
                return data_dict
        else:
            # otherwise print error
            print('Requested network {} is missing'.format(net_name))
            return data_dict
    # for the approx21(stuff)+moreStuff networks:
    if re.findall('approx21',net_name):
        for i,term in enumerate(net_name.split('_')):
            if i==0 and term=="approx21":
                data_dict = read_file(data_dict,file_path,term+'.net')
            elif term=="plus":
                continue
            elif len(term)>3:
                if term[-4:]=='.net':
                    term=term[0:-4]
                key, values = process_line(term)
                if key not in data_dict:
                    data_dict[key] = []
                data_dict[key].extend(values)
        return data_dict
    # for all other cases:
    with open(os.path.join(file_path,net_name), "r") as file:
        first_line = True
        for line in file:
            if line=='\n':
                continue
            if line.lstrip(' ')[0]=='!':
                continue
            if re.findall('include',line):
                data_dict = read_file(data_dict,file_path,line.strip().split()[1][1:-1])
                continue
            if re.findall('add_iso\(',line) and re.findall('\)',line):
                key, values = process_line(line[line.find("(")+1:line.find(")")])
                if key not in data_dict:
                    data_dict[key] = []
                data_dict[key].extend(values)
                break
            if re.findall('add_isos\(',line) and re.findall('\)',line):
                for term in line[line.find("(")+1:line.find(")")].split(', '):
                    key, values = process_line(term)
                    if key not in data_dict:
                        data_dict[key] = []
                    data_dict[key].extend(values)
                break
            if first_line and (re.findall("add_isos\(",line) or re.findall("add_isos_and_reactions\(",line)):
                first_line = False
                continue
            if re.findall("\)",line):
                break
            
            key, values = process_line(line.split('!')[0])
            if key not in data_dict:
                data_dict[key] = []
            data_dict[key].extend(values)

    return data_dict

def process_reac(line):
    # reaction lines have each side with ' ' (one space) separating reaction elements.
    # between them, are at least '  ' (double space), sometimes more.
    # It seems to be reasonable to assume no temporary element names, so only XX### (no XXX###).
    parts = line.strip().split('  ').remove('')
    
    lhs = parts[0]
    rhs = parts[1]
    
    def tryIt(item,pat):
        try: 
             return re.findall(pat,item)[0]
        except IndexError: 
             return re.findall(pat,item)
    
    lhs = [tryIt(it,r"[a-zA-Z]+\d+") or it for it in lhs.split(' ')]
    rhs = [tryIt(it,r"[a-zA-Z]+\d+") or it for it in rhs.split(' ')]
    
    # def breakDown(item):
    #     if len(item) == 1 and tryIt(item,r"[a-zA-Z]+")<=2:
    #         return item #we're done
    #     if len(item) == 1 and tryIt(item,r"[a-zA-Z]+")>2:
    #         return breakDown(item[0],item[1:]) #we've only just begun
    #     itty = []
    #     for it in item:
    #         if len(tryIt(it,r"[a-zA-Z]+"))>2:
    #             return breakDown(it[0],it[1:])
    #         itty.append(it)
    #     return itty
    
    # def breakDown(item):
    #     match len(tryIt(item,r"[a-zA-Z]+")):
    #         case 4:
    #             return [item[0],item[1],item[2:]]
    #         case 3:
    #             return [item[0],item[1:]]
    #         case 2:
    #             return item
    
    # for it in lhs:
    #     lengthy = len(re.findall((r"[a-zA-Z]+\d+",it)))
    #     if len(re.findall((r"[a-zA-Z]+\d+",it))) > 2:
    #         while
    #         s
    
    # Check for "XX### X X X XX##" format (i.e., 'fe45 p p p ti42')
    if len(parts) == 1 and parts[0][-1].isalpha():
        return parts[0], []
    
    # Check for "X#" format (i.e., 'fe54', 'h2'
    if len(parts) == 1:
        return re.findall("[a-zA-Z]+", parts[0])[0], [int(re.findall(r'\d+', parts[0])[0])]
    
    # Check for "X #1 #2" format (i.e., 'fe 52 54', 'he 3 4')
    if len(parts) == 3 and parts[0].isalpha() and all(part.isdigit() for part in parts[1:]):
        x = parts[0]
        start, end = int(parts[1]), int(parts[2])
        return x, list(range(start, end + 1))
    
    # If the format is not recognized, return an empty list for the values
    return "", []
    