#   /\_/\  
#  ( o.o ) 
#   > ^ < 

#SET UP
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression

#SETUP


import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random
import copy
import numpy as np
import sys

from scipy import stats 
from scipy.stats import norm
from scipy.stats import qmc

import itertools

import csv
import pandas as pd
import gc


#TOMIMOTO INITIALISE


# Elongation
# if st_d = num_stem --> Structured Elongation
# if st_d = 0 --> Stochastic Elongation

def mutInStemCells(num_stem, t, mu_0, Genom, st_d): #Num_of_stemcells, Num_of_division_in_elongation, mutation_rate, Genom_size, Num_of_stemcells - gamma

    ##stemCells = ''
    tCells = [] #list for history of stemcells in meristem
    stemCells = [[i,[]] for i in range(num_stem)] #[(ID number）,[list of mutations]]

    tCells = copy.deepcopy(stemCells)

    for k in range(t): #one cycle of division

        #one daughter cells
        ccells = [[0,0] for i in range(num_stem)]
        for i in range(num_stem):
            ##m = ""  
            m = random.random()
            if m > mu_0: #non mutaion
                ccells[i] = [stemCells[i][0], stemCells[i][1]]
            else: #mutation
                lst = copy.deepcopy(stemCells[i][1])
                lst.append(random.randint(1,Genom)) #site of Genom in which mut occure
                ccells[i] = [stemCells[i][0], lst]

        #the other daughter cells
        ccells2 = [[0,0] for i in range(num_stem)]
        for i in range(num_stem):
            m = random.random()
            if m > mu_0: #non mutaion
                ccells2[i] = [stemCells[i][0], stemCells[i][1]]
            else: #mutaion
                lst = copy.deepcopy(stemCells[i][1])
                lst.append(random.randint(1,Genom)) #site of Genom in which mut occure
                ccells2[i] = [stemCells[i][0], lst]

        #step1 + step2: choose cells of next meristem
        stemCells = random.sample(ccells + random.sample(ccells2, st_d), num_stem) #step1:random.sample(ccells2, st_d) by gammma(num_stem - st_d)

        tCells = tCells + copy.deepcopy(stemCells)

    return tCells

def mutInBrStemCells(num_stem, stemCells, t, mu_0, Genom, st_d):
    tCells = copy.deepcopy(stemCells)
    cells = copy.deepcopy(stemCells)  # to prevent reference of stemCells in for_loop

    ccells = [[0, 0] for _ in range(num_stem)]
    ccells2 = [[0, 0] for _ in range(num_stem)]

    for k in range(t):  # one cycle of div
        # one daughter cell
        for i in range(num_stem):
            m = random.random()
            if m > mu_0:  # non mutation
                ccells[i] = [cells[i][0], cells[i][1]]
            else:  # mutation
                lst = copy.deepcopy(cells[i][1])
                lst.append(random.randint(1, Genom))  # site of Genom in which mut occur
                ccells[i] = [cells[i][0], lst]

        # the other daughter cell
        for i in range(num_stem):
            m = random.random()
            if m > mu_0:  # non mutation
                ccells2[i] = [cells[i][0], cells[i][1]]
            else:  # mutation
                lst = copy.deepcopy(cells[i][1])
                lst.append(random.randint(1, Genom))  # site of Genom in which mut occur
                ccells2[i] = [cells[i][0], lst]

        # step1 + step2: choose cells of next meristem
        cells = random.sample(ccells + random.sample(ccells2, st_d), num_stem)  # step1:random.sample(ccells2, st_d) by gamma(num_stem - st_d)

        tCells = tCells + copy.deepcopy(cells)

    return tCells

# Branching
# bias_d = 0.5 ->BIASED
# bias_d = 10 -> UNBIASED

#for small sigma
def WNdist_f(x_value, mu_value, sig_value, n_value): # x, mean, var,

    f_func = 1/(sig_value*math.sqrt(2*math.pi))*np.sum([math.exp(-(x_value + 2*math.pi*k - mu_value)**2/(2*sig_value**2)) for k in range(-n_value,n_value+1)]) #norm.pdf(x, loc, scale)

    return f_func

#for large sigma
def WNdist_g(x_value, mu_value, sig_value, n_value): # x, mean, var,

    g_func = 1/(2*math.pi)*(1 + 2*np.sum([math.exp(-1/2*sig_value**2)**(k**2)*math.cos(k*(x_value-mu_value)) for k in range(1,n_value+1)])) #norm.pdf(x, loc, scale)

    return g_func


def WNdist(x_value, mu_value, sig_value):

    n_value = 5

    if sig_value < 3: #for small sigma
        dist_func = WNdist_f(x_value, mu_value, sig_value, n_value)

    elif sig_value >= 3: #for large sigma
        dist_func = WNdist_g(x_value, mu_value, sig_value, n_value)

    return dist_func


def weightList(num_stem, num_div, bias_d): # Num of stem cells, Num of division, sigma

    cell_dt = 2*math.pi/(num_stem*(2**num_div)) #distance between cells
    mean_num = random.uniform(0,2*math.pi) #center of dist

    wList = np.array([WNdist(i*cell_dt, mean_num, bias_d) for i in range(num_stem*(2**num_div))])

    return wList/np.sum(wList)

#Sampling of stem cells in Branching

def sample_mutations(stem_cells, num_stem, num_div, mu_0, Genom, bias_d): #　cells in apical meristem, Num of stem cells, Num of cell divisions in branching, Mutation rate, Genome size, sigma

    wList = weightList(num_stem, num_div, bias_d)

    cells = copy.deepcopy(stem_cells[-num_stem:])

    for k in range(0,num_div):

        ccells = [[0,0] for i in range(len(cells))]
        ccells2 = [[0,0] for i in range(len(cells))]

        #one daughter cells
        for i in range(len(cells)):
            m = random.random()
            if m > mu_0: #non mutaion
                ccells[i] = [copy.deepcopy(cells[i][0]), copy.deepcopy(cells[i][1])]
            else: #mutation
                lst = ''
                lst = copy.deepcopy(cells[i][1])
                lst.append(random.randint(1,Genom)) #site of Genom in which mut occure
                ccells[i] = [copy.deepcopy(cells[i][0]), lst]

        #the other daughter cells
        for i in range(len(cells)):
            m = random.random()
            if m > mu_0: #non mutaion
                ccells2[i] = [copy.deepcopy(cells[i][0]), copy.deepcopy(cells[i][1])]
            else: #mutaionI 
                lst = copy.deepcopy(cells[i][1])
                lst.append(random.randint(1,Genom)) #site of Genom in which mut occure
                ccells2[i] = [copy.deepcopy(cells[i][0]), lst]

        lst2 = [] #list for daughter cells
        for p,q in zip(ccells,ccells2):
            lst2.append(p)
            lst2.append(q)

        cells = copy.deepcopy(lst2) #daughter cells arranged in genealogical order.

        #check
        ##print(cells)

    sampleNum = np.random.choice([i for i in range(len(cells))], size=num_stem, replace=False, p=wList) #sampling of cells for axillary meristem
    sample = [cells[i] for i in sampleNum]

    return sample
def pickUpMut(allmutList):
    mutList = []
    for i in allmutList:
        for j in i[1]:
            mutList.append(j)

    reAllmutList = list(set(mutList)) #eliminate overlapping of mutation
    reAllmutList.sort()

    return reAllmutList

#making mutation matrix from history of meristems
def makeMutMatrix(List, mutList, num_stem): # List_of_meristems、list_of_all_mutated_site、Num_of_stemcells
    preMutMatrix = [[] for i in range(len(List))]
    for i,j in enumerate(List):
        for k in mutList:
            if k in j[1]:
                preMutMatrix[i].append(1)
            else:
                preMutMatrix[i].append(0)

    #check
    for i,j in zip(List, preMutMatrix):
        if len(i[1]) != sum(j):
            print("mutation occured at the same site")
            break

    #sum into meristems
    np_mer_mutMatrix = np.array([0 for i in range(int(len(List)/num_stem))])
    for i in range(int(len(preMutMatrix)/num_stem)): #int(len(preMutMatrix)/num_stem) = number of history of merstems
        mer_mutVector = np.zeros(len(mutList),dtype = int)
        for j in range(num_stem):
            arraylist = np.array(preMutMatrix[i*num_stem +j])
            mer_mutVector = mer_mutVector + arraylist

        if i == 0:
            np_mer_mutMatrix = mer_mutVector
        else:
            np_mer_mutMatrix = np.vstack((np_mer_mutMatrix, mer_mutVector))

    return np_mer_mutMatrix

#null regression
def reg1dim(x, y):
    a = np.dot(x, y)/ (x ** 2).sum()
    return a

#_for saving the data

#for saving chage the shape of list
def add_None_func(ori_List):
    max_length = max([len(i[0]) for i in ori_List]) #scaling of maximum length
    none_List = np.array([[None] for _ in range(len(ori_List[0]))]) #len(ori_List[0]): len(copy_List) due to the type of data

    pre_List = []
    #count = 0 #check
    for j_ in ori_List: #loop by NumTime
        #count += 1 #check
        copy_List = copy.deepcopy(j_)

        while len(copy_List[0]) < max_length: #len(i_[0]): number of mutations
            copy_List = np.hstack((copy_List ,none_List))

        pre_List.append(copy.deepcopy(copy_List))
        #print(count)#check

    return pre_List

#for Refrost data
def non_None(ori_List):
    non_matrix = []
    for i in ori_List:
        pre_matrix = []
        for j in i:
            pre_matrix.append([k for k in j if k is not None])
        non_matrix.append(np.array(copy.deepcopy(pre_matrix)))

    return non_matrix


#LATIN HYPERCUBE


#Constants
NumStem = 5  # Stem cells constant
nDiv = 1  # Elongation division constant (Re)
NumDiv = 7  # Branching division constant (Rb)
NumTime = 3  # 1 for loop per trial
nElim = 0 #threshold of detection

import numpy as np
from scipy.stats import qmc

def latin_hypercube_sampling(num_samples):
    # Define the parameter levels
    mu_0 = [6,60]
    GenSize = [500000000]  
    StD = [0, 5]  # Elongation param-> Stochastic = num_stem, Structured = 0
    biasVar = [0.5, 10]  # Branching param -> Unbiased= 0.5, Biased = 10
    tree_topologies = ["bS4", "ubS4", "bL4", "ubL4", "bS6", "ubS6", "bL6", "ubL6", "bS8", "ubS8", "bL8", "ubL8", "bS10", "ubS10", "bL10", "ubL10", "bS12", "ubS12", "bL12", "ubL12"]

    
    num_params = 4  
    sampler = qmc.LatinHypercube(d=num_params, optimization="random-cd")
    lhs_samples = sampler.random(n=num_samples)

    
    param_samples = []
    for i in range(num_samples):
        param_sample = {
            "mu_0": mu_0[int(lhs_samples[i, 0] * len(mu_0)) % len(mu_0)],  # Ensuring proper indexing
            "GenSize": GenSize[0],  # Fixed genome size
            "tree_topology": tree_topologies[int(lhs_samples[i, 1] * len(tree_topologies)) % len(tree_topologies)],
            "StD": StD[int(lhs_samples[i, 2] * len(StD)) % len(StD)],
            "biasVar": biasVar[int(lhs_samples[i, 3] * len(biasVar)) % len(biasVar)]
        }
        param_samples.append(param_sample)

    return param_samples



#Tree Topologies


# tree config
tree_topologies_dict = {
    "test":{
        "numBranch":8,
        "age": 200,
        "s10": 1,
        "b11": 36,"bb11": 130,"b12": 11,"bb12": 103,"b13": 11,"bb13": 98,"b14":93,
        "s40": 36, "b41": 84, "s41": 77, "b42": 80, "s42": 8, "b43": 70, "s44":64
    },
        
        
    "bS4": {
        "numBranch": 4,
        "age": 123,
        "s10": 10,
        "b11": 75, "bb11": 38, "b12": 38,
        "s40": 75, "b41": 38, "s41": 38
    },
    "ubS4": {
        "numBranch": 4,
        "age": 236,
        "s10": 10,
        "b11": 75, "bb11": 38, "b12": 75, "bb12": 38, "b13": 38,
        "s40": 38
    },
    "bL4": {
        "numBranch": 4,
        "age": 100,
        "s10": 10,
        "b11": 30, "bb11": 60, "b12": 60,
        "s40": 30, "b41": 60, "s41": 60
    },

    "ubL4": {
        "numBranch": 4,
        "age": 190,
        "s10": 10,
        "b11": 30, "bb11": 60, "b12": 30, "bb12": 60, "b13": 60,
        "s40": 60
    },
    "bS6": {
        "numBranch": 6,
        "age": 147,
        "s10": 10,
        "b11": 42, "bb11": 21, "b12": 42, "bb12": 21, "b13": 21,
        "s40": 42, "b41": 21, "s41": 42, "b42": 21, "s42": 42
    },
    "ubS6": {
        "numBranch": 6,
        "age": 262,
        "s10": 10,
        "b11": 42, "bb11": 21, "b12": 42, "bb12": 21, "b13": 42, "bb13": 21, "b14": 42, "bb14": 21, "b15": 21,
        "s40": 21
    },
    "bL6": {
        "numBranch": 6,
        "age": 124,
        "s10": 10,
        "b11": 19, "bb11": 38, "b12": 19, "bb12": 38, "b13": 38,
        "s40": 19, "b41": 38, "s41": 19, "b42": 38, "s42": 38
    },
    "ubL6": {
        "numBranch": 6,
        "age": 238,
        "s10": 10,
        "b11": 19, "bb11": 38, "b12": 19, "bb12": 38, "b13": 19, "bb13": 38, "b14": 19, "bb14": 38, "b15": 38,
        "s40": 38
    },
    "bS8": {
        "numBranch": 8,
        "age": 145,
        "s10": 10,
        "b11": 30, "bb11": 15, "b12": 30, "bb12": 15, "b13": 30, "bb13": 15, "b14": 15,
        "s40": 30, "b41": 15, "s41": 30, "b42": 15, "s42": 30, "b43": 15, "s43": 15
    },
    "ubS8": {
        "numBranch": 8,
        "age": 280,
        "s10": 10,
        "b11": 30, "bb11": 15, "b12": 30, "bb12": 15, "b13": 30, "bb13": 15, "b14": 30, "bb14": 15, "b15": 30, "bb15": 15, "b16": 30, "bb16": 15, "b17": 15,
        "s40": 15
    },
    "bL8": {
        "numBranch": 8,
        "age": 127,
        "s10": 10,
        "b11": 13, "bb11": 26, "b12": 13, "bb12": 26, "b13": 13, "bb13": 26, "b14": 26,
        "s40": 13, "b41": 26, "s41": 13, "b42": 26, "s42": 13, "b43": 26, "s43": 26
    },
    "ubL8": {
        "numBranch": 8,
        "age": 244,
        "s10": 10,
        "b11": 13, "bb11": 26, "b12": 13, "bb12": 26, "b13": 13, "bb13": 26, "b14": 13, "bb14": 26, "b15": 13, "bb15": 26, "b16": 13, "bb16": 26, "b17": 26,
        "s40": 26
    },
    "bS10": {
        "numBranch": 10,
        "age": 154,
        "s10": 10,
        "b11": 24, "bb11": 12, "b12": 24, "bb12": 12, "b13": 24, "bb13": 12, "b14": 24, "bb14": 12, "b15": 12,
        "s40": 24, "b41": 12, "s41": 24, "b42": 12, "s42": 24, "b43": 12, "s43": 24, "b44": 12, "s44": 12
    },
    "ubS10": {
        "numBranch": 10,
        "age": 310,
        "s10": 10,
        "b11": 24, "bb11": 12, "b12": 24, "bb12": 12, "b13": 24, "bb13": 12, "b14": 24, "bb14": 12, "b15": 24, "bb15": 12, "b16": 24, "bb16": 12, "b17": 24, "bb17": 12, "b18": 24, "bb18": 12, "b19": 12,
        "s40": 12
    },
    "bL10": {
        "numBranch": 10,
        "age": 142,
        "s10": 10,
        "b11": 11, "bb11": 22, "b12": 11, "bb12": 22, "b13": 11, "bb13": 22, "b14": 11, "bb14": 22, "b15": 22,
        "s40": 11, "b41": 22, "s41": 11, "b42": 22, "s42": 11, "b43": 22, "s43": 11, "b44": 22, "s44": 22
    },
    "ubL10": {
        "numBranch": 10,
        "age": 274,
        "s10": 10,
        "b11": 11, "bb11": 22, "b12": 11, "bb12": 22, "b13": 11, "bb13": 22, "b14": 11, "bb14": 22, "b15": 11, "bb15": 22, "b16": 11, "bb16": 22, "b17": 11, "bb17": 22, "b18": 11, "bb18": 22, "b19": 22,
        "s40": 22
    },
    "bS12": {
        "numBranch": 12,
        "age": 145,
        "s10": 10,
        "b11": 18, "bb11": 9, "b12": 18, "bb12": 9, "b13": 18, "bb13": 9, "b14": 18, "bb14": 9, "b15": 18, "bb15": 9, "b16": 9,
        "s40": 18, "b41": 9, "s41": 18, "b42": 9, "s42": 18, "b43": 9, "s43": 18, "b44": 9, "s44": 18, "b45": 9, "s45": 9
    },
    "ubS12": {
        "numBranch": 12,
        "age": 280,
        "s10": 10,
        "b11": 18, "bb11": 9, "b12": 18, "bb12": 9, "b13": 18, "bb13": 9, "b14": 18, "bb14": 9, "b15": 18, "bb15": 9, "b16": 18, "bb16": 9, "b17": 18, "bb17": 9, "b18": 18, "bb18": 9, "b19": 18, "bb19": 9, "b20": 18, "bb20": 9, "b21": 9,
        "s40": 9
    },
    "bL12": {
        "numBranch": 12,
        "age": 130,
        "s10": 10,
        "b11": 8, "bb11": 16, "b12": 8, "bb12": 16, "b13": 8, "bb13": 16, "b14": 8, "bb14": 16, "b15": 8, "bb15": 16, "b16": 16,
        "s40": 8, "b41": 16, "s41": 8, "b42": 16, "s42": 8, "b43": 16, "s43": 8, "b44": 16, "s44": 8, "b45": 16, "s45": 16
    },
    "ubL12": {
        "numBranch": 12,
        "age": 250,
        "s10": 10,
        "b11": 8, "bb11": 16, "b12": 8, "bb12": 16, "b13": 8, "bb13": 16, "b14": 8, "bb14": 16, "b15": 8, "bb15": 16, "b16": 8, "bb16": 16, "b17": 8, "bb17": 16, "b18": 8, "bb18": 16, "b19": 8, "bb19": 16, "b20": 8, "bb20": 16, "b21": 16,
        "s40": 16
    }
}



#Tree list/dict function


def create_tree_list_and_dict(tree_topology):
    
    keys = list(tree_topology.keys())  
    tree_list = []
    tree_dict = {
        "s10": tree_topology["s10"],
        "Ri": [],
        "Rt": [],
        "Li": [],
        "Lt": []
    }

    for key in keys:
        if key.startswith('s') or key.startswith('b'):
            tree_list.append(tree_topology[key])

    
    s40_index = keys.index('s40')
    right_keys = keys[1:s40_index]  
    for i, key in enumerate(right_keys):
        if key.startswith('b') and not key.startswith('bb'):
            if i == len(right_keys) - 1 or not right_keys[i + 1].startswith('bb'):
                tree_dict["Rt"].append(tree_topology[key])  # Last 'b' before s40
            else:
                tree_dict["Ri"].append(tree_topology[key])
        elif key.startswith('bb'):
            tree_dict["Rt"].append(tree_topology[key])

    
    left_keys = keys[s40_index:]  
    for i, key in enumerate(left_keys):
        if key.startswith('s') and not key.startswith('b'):
            if i == len(left_keys) - 1 or (i + 1 < len(left_keys) and not left_keys[i + 1].startswith('b')):
                tree_dict["Lt"].append(tree_topology[key])  # Last 'b' on the left side
            else:
                tree_dict["Li"].append(tree_topology[key])
        elif key.startswith('b'):
            tree_dict["Lt"].append(tree_topology[key])

    return tree_list, tree_dict, tree_topology["numBranch"], tree_topology["age"]


#Ave Meristem Time Function 


def calculate_ave_meristem_time(tree_dict, numBranch):
    # find branches
    root = tree_dict["s10"]
    right_internal = tree_dict["Ri"]
    right_terminal = tree_dict["Rt"]
    left_internal = tree_dict["Li"]
    left_terminal = tree_dict["Lt"]

    # root tip dist calc
    def calculate_root_to_tip(internal_branches, terminal_branches):
        total_distance = 0
        for i, terminal in enumerate(terminal_branches):
            # Sum up internal branches up to the corresponding terminal branch
            internal_sum = sum(internal_branches[:i + 1])
            total_distance += root + internal_sum + terminal
        return total_distance

    right_total_distance = calculate_root_to_tip(right_internal, right_terminal)
    left_total_distance = calculate_root_to_tip(left_internal, left_terminal)

    
    total_distance = right_total_distance + left_total_distance
    ave_meristem_time = total_distance / numBranch

    return ave_meristem_time


#SIMULATION FUNCTION!!!!!!!!! :0


def simulate_somatic_mutations(tree_dict, NumStem, NumTime, mu_0, GenSize, StD, NumDiv, nDiv, biasVar):
    NumTime_br_brmutMatrix = []
    µ_peryear = mu_0*10**(-10) * GenSize
    µ14 = (1 - (1 - µ_peryear) ** (1 / nDiv))
    µ13 = (1 - (1 - µ_peryear) ** (1 / nDiv))

    for num_time in range(NumTime):  
        # Simulation begin
        List10 = mutInStemCells(NumStem, tree_dict["s10"], µ14 / 2 + µ13 / 2, GenSize, StD)
        sampleCells10 = sample_mutations(List10, NumStem, NumDiv, µ14, GenSize, biasVar)
        print(List10)

        List_br = []
        List_left = []

        # R branch (list_br)
        if len(tree_dict["Ri"]) == 0:
            for i in range(len(tree_dict["Rt"])):
                List_br.append(mutInBrStemCells(NumStem, sampleCells10, tree_dict["Rt"][i], µ14, GenSize, StD))
        else:
            for i in range(len(tree_dict["Ri"])):
                if i == 0:
                    List_br.append(mutInBrStemCells(NumStem, sampleCells10, tree_dict["Ri"][i], µ14, GenSize, StD))
                else:
                    List_br.append(mutInBrStemCells(NumStem, List_br[-2][-NumStem:], tree_dict["Ri"][i], µ14, GenSize, StD))
                
                sampleCells_br = sample_mutations(List_br[-1], NumStem, NumDiv, µ14, GenSize, biasVar)
                List_br.append(mutInBrStemCells(NumStem, sampleCells_br, tree_dict["Rt"][i], µ14, GenSize, StD))
                
            # force last branch
            if len(List_br) > 1:
                #last_sampleCells_br = sample_mutations(List_br[-1], NumStem, NumDiv, µ14, GenSize, biasVar)
                #List_br.append(mutInBrStemCells(NumStem, last_sampleCells_br, tree_dict["Rt"][-1], µ14, GenSize, StD))
                List_br.append(mutInBrStemCells(NumStem, List_br[-2][-NumStem:], tree_dict["Rt"][-1], µ14, GenSize, StD))


        print(f"right branches: {len(List_br)}")
        sampleCells10_2 = sample_mutations(List10, NumStem, NumDiv, µ13, GenSize, biasVar)
        # L branch
        if len(tree_dict["Li"]) == 0:
            for i in range(len(tree_dict["Lt"])):
                List_left.append(mutInBrStemCells(NumStem, sampleCells10_2, tree_dict["Lt"][i], µ13, GenSize, StD))
        else:
            for i in range(len(tree_dict["Li"])):
                if i == 0:
                    List_left.append(mutInBrStemCells(NumStem, sampleCells10_2, tree_dict["Li"][i], µ13, GenSize, StD))
                else:
                    List_left.append(mutInBrStemCells(NumStem, List_left[-2][-NumStem:], tree_dict["Li"][i], µ13, GenSize, StD))
                
                sampleCells_left = sample_mutations(List_left[-1], NumStem, NumDiv, µ13, GenSize, biasVar)
                List_left.append(mutInBrStemCells(NumStem, sampleCells_left, tree_dict["Lt"][i], µ13, GenSize, StD))
            # force last branch
            if len(List_left) > 1:
                #last_sampleCells_left = sample_mutations(List_left[-1], NumStem, NumDiv, µ13, GenSize, biasVar)
                #List_left.append(mutInBrStemCells(NumStem, last_sampleCells_left, tree_dict["Lt"][-1], µ13, GenSize, StD))
                List_left.append(mutInBrStemCells(NumStem, List_left[-2][-NumStem:], tree_dict["Lt"][-1], µ13, GenSize, StD))
     
        print(f"left branches : {len(List_left)}")

        # combine L + R
        combined_list = [List10] + List_br + List_left 
        print(len(combined_list))

        flattened_combined_list = [item for sublist in combined_list for item in sublist]

        allMutations = pickUpMut(flattened_combined_list)
        print("all mut:", len(allMutations))

        # check in range 
        if len(List_br) > 0 and len(List_left) > 0:
            endallMutations = pickUpMut([item for sublist in List_br[1::2] for item in sublist] + 
                                        [item for sublist in List_left[1::2] for item in sublist] + 
                                        List_br[-1] + List_left[-1])
        elif len(List_br) > 0:
            endallMutations = pickUpMut([item for sublist in List_br[1::2] for item in sublist] + List_br[-1])
        elif len(List_left) > 0:
            endallMutations = pickUpMut([item for sublist in List_left[1::2] for item in sublist] + List_left[-1])
        else:
            endallMutations = []

        print(f"total mut {len(endallMutations)}")

        # Mutation Matrices
        ListTo_br_mut_matrices = []
        for i in range(0, len(List_br), 2):
            combined_list = List10 + [item for sublist in List_br[:i+2] for item in sublist]
            ListTo_br_mut_matrices.append(makeMutMatrix(combined_list, allMutations, NumStem))

        # force last branch
        if len(List_br) % 2 == 0 and len(List_br) > 1:
            combined_list = List10 + [item for sublist in List_br for item in sublist]
            ListTo_br_mut_matrices.append(makeMutMatrix(combined_list, allMutations, NumStem))

        ListTo_left_mut_matrices = []
        for i in range(0, len(List_left), 2):
            combined_list = List10 + [item for sublist in List_left[:i+2] for item in sublist]
            ListTo_left_mut_matrices.append(makeMutMatrix(combined_list, allMutations, NumStem))

        # Force Last branch
        if len(List_left) % 2 == 0 and len(List_left) > 1:
            combined_list = List10 + [item for sublist in List_left for item in sublist]
            ListTo_left_mut_matrices.append(makeMutMatrix(combined_list, allMutations, NumStem))

        pre_br_brmutMatrix = np.array([mat[-1, :] for mat in ListTo_br_mut_matrices + ListTo_left_mut_matrices])
        br_brmutMatrix = np.array([i for i in pre_br_brmutMatrix.T if sum(i) != 0]).T

        NumTime_br_brmutMatrix.append(br_brmutMatrix)


    if not NumTime_br_brmutMatrix:
        raise ValueError("NumTime_br_brmutMatrix is empty. Simulation FAILED :(")

    # save and load results
    result_sim = NumTime_br_brmutMatrix
    np.save('br_mutmatrix.npy', add_None_func(result_sim))
    test_matrixsim = non_None(np.load('br_mutmatrix.npy', allow_pickle=True))
    print(np.array(test_matrixsim[0]).shape[0])

    return test_matrixsim


#Check for Bugs!

def is_br_mutmatrix_empty(resultList_2):
    for matrix in resultList_2:
        if len(matrix) == 0 or matrix.shape[0] == 0:
            return True
    return False



#Count Mutations Function!

def mut_dist_func(resultList_2):
    try:
        numBranch = len(resultList_2[0])
        print(f"Number of branches (numBranch): {numBranch}")

        mutShapeTemplate = np.array(list(itertools.product([0, 1], repeat=numBranch)))
        print(f"mutShapeTemplate shape: {mutShapeTemplate.shape}")

        NumStemList = []
        Nummut_freqTimeFreq = []  
        Numallmut_freqDifMatrix = []  

        # The pattern of expansion of somatic mutations
        NumAveallMutDist = []
        NumAveSDallMutDist = []  # for SD

        NumAveallMutFreqDist = []  # freq

        for NumStem in [5]:
            NumTime = len(resultList_2)

            for num_time in range(NumTime):
                try:
                    br_brmutMatrix = resultList_2[num_time]
                    print(f"br_brmutMatrix shape: {br_brmutMatrix.shape}")

                    # The distribution pattern of somatic mutations
                    # focusing on mutations at the tip of branches
                    mutShapeMatrix = [0 for i in range(2 ** numBranch)]
                    br_brmut01Matrix = np.nan_to_num(br_brmutMatrix / br_brmutMatrix)  # convert to 01 (changing Nan, resulting from 0/0, into 0), row: branch, column: mutation

                    for icut, i in enumerate(mutShapeTemplate):  # mutShapeTemplate is a list for all possible combinations
                        for k in range(len(br_brmutMatrix[0])):
                            if np.all(br_brmut01Matrix[:, k] == i):
                                mutShapeMatrix[icut] += 1

                    if num_time == 0:
                        allMutDist = np.array(mutShapeMatrix)
                        totalMutDist = [mutShapeMatrix]  # for SD
                    else:
                        allMutDist = allMutDist + np.array(mutShapeMatrix)
                        totalMutDist.append(mutShapeMatrix)  # for SD

                    # The distribution pattern of somatic mutations
                    # freq_(枝末端間の変異細胞数の average)
                    mutShapeFreqMatrix = [[0 for i in range(sum(j) * NumStem)] for j in mutShapeTemplate]  # for different Ave_mutations, sum(j)*NumStem: number of average variations

                    for icut, i in enumerate(mutShapeTemplate):  # mutShapeTemplate is a list for all possible combinations
                        for j in br_brmutMatrix.T:
                            if sum(j) != 0 and np.all(np.nan_to_num(j / j) == i):  # eliminate [0,0,...,0] all 0 list
                                mutShapeFreqMatrix[icut][sum(j) - 1] += 1  # sum(j): sum of mutated cells, sum(j)-1: index of list

                    if num_time == 0:
                        allMutFreqDist = copy.deepcopy(mutShapeFreqMatrix)
                    else:
                        for icut, i in enumerate(copy.deepcopy(mutShapeFreqMatrix)):
                            for jcut, j in enumerate(i):
                                allMutFreqDist[icut][jcut] += j

                except IndexError:
                    print(f"IndexError encountered at NumTime {num_time}. Skipping this trial.")
                    continue

            # Distribution of mutations
            AveallMutDist = allMutDist / NumTime  # The pattern of expansion of somatic mutations
            NumAveallMutDist.append(AveallMutDist)  # for different NumStem

            # Standard deviation of dist
            sumtotalMutDist = np.array([float(0) for _ in mutShapeTemplate])
            for k in totalMutDist:
                sumtotalMutDist += (np.array(k) - AveallMutDist) ** 2  # branch-wise
            NumAveSDallMutDist.append((sumtotalMutDist / len(totalMutDist)) ** (1 / 2))  # standard deviation for each branch

            # freq_pattern_Distribution of mutations
            AveallMutFreqDist = [[i / NumTime for i in j] for j in allMutFreqDist]  # The pattern of expansion of somatic mutations with freq
            NumAveallMutFreqDist.append(AveallMutFreqDist)

            NumStemList.append(NumStem)

        # Calculate shared, unique, and total mutations
        shared_mutations = {}
        unique_mutations = [0 for _ in range(numBranch)]
        mutations_counts = [0 for _ in range(numBranch)]

        mutDist_TemFreqList = [[i, j] for i, j in zip(mutShapeTemplate.tolist(), NumAveallMutFreqDist[0])]
        mutDist_TemFreqList_sort_key = sorted(mutDist_TemFreqList, key=lambda x: sum(x[0]))  # sort by template

        mutDist_TemFreqList_sort_key_short = [i for i in mutDist_TemFreqList_sort_key if sum(i[1]) >= 0.1]  # removing mutation with freq < 0.1

        for pattern_freq in mutDist_TemFreqList_sort_key_short:
            pattern, freq_list = pattern_freq[0], pattern_freq[1]
            # Identify shared mutations
            if sum(pattern) > 1:  # Mutation is shared between branches
                branches = tuple([i + 1 for i, x in enumerate(pattern) if x == 1])
                shared_mutations[branches] = shared_mutations.get(branches, 0) + sum(freq_list)

            # Count unique mutations
            if sum(pattern) == 1:  # Mutation is unique to one branch
                branch = pattern.index(1)
                unique_mutations[branch] += sum(freq_list)

            # Count total mutations for each branch
            for i, has_mutation in enumerate(pattern):
                if has_mutation == 1:
                    mutations_counts[i] += sum(freq_list)

        print(f"Unique mutations: {unique_mutations}")
        print(f"Mutation counts: {mutations_counts}")

        return NumAveallMutFreqDist, NumAveSDallMutDist, shared_mutations, unique_mutations, mutations_counts

    except IndexError:
        print("IndexError encountered in mut_dist_func. Skipping this function.")
        return [], [], {}, [], []


#Make dist matrices function


def gen_matrices(tree_dict, unique_mutations, shared_mutations, mutations_counts):
    numBranch = len(mutations_counts)
    print(f"Number of branches (numBranch): {numBranch}")

    genetic_distance_matrix = np.zeros((numBranch, numBranch))

    # pairwise dist calc
    for i in range(numBranch):
        for j in range(numBranch):
            if i == j:
                genetic_distance_matrix[i][j] = 0
            else:
                shared = shared_mutations.get((min(i+1, j+1), max(i+1, j+1)), 0)
                genetic_distance_matrix[i][j] = mutations_counts[i] + mutations_counts[j] - shared

    # branch length extract
    root_length = tree_dict['s10']
    right_internal = tree_dict['Ri']
    right_terminal = tree_dict['Rt']
    left_internal = tree_dict['Li']
    left_terminal = tree_dict['Lt']

    # root tip dist R
    right_root_to_tip = []
    for i in range(len(right_terminal)):
        distance = root_length + sum(right_internal[:i+1]) + right_terminal[i]
        right_root_to_tip.append(distance)

    # root tip dist L
    left_root_to_tip = []
    for i in range(len(left_terminal)):
        distance = root_length + sum(left_internal[:i+1]) + left_terminal[i]
        left_root_to_tip.append(distance)

    root_to_tip_distances = right_root_to_tip + left_root_to_tip

    # PHYSICAL DIST !
    physical_distance_matrix = np.zeros((numBranch, numBranch))

    # pairwise physical
    for i in range(numBranch):
        for j in range(numBranch):
            if i == j:
                physical_distance_matrix[i][j] = 0
            else:
                physical_distance_matrix[i][j] = abs(root_to_tip_distances[i] - root_to_tip_distances[j])

    # make sure length match
    if len(root_to_tip_distances) != numBranch:
        raise ValueError("The length of root_to_tip_distances does not match numBranch")

    print(f"Root to tip distances: {root_to_tip_distances}")
    print(f"Length of root_to_tip_distances: {len(root_to_tip_distances)}")

    return genetic_distance_matrix, physical_distance_matrix


#Calc Variants + mut_rate function


def calc_variants(genetic_matrix, physical_matrix, tree_dict, GenSize):
    physical_distances = []
    genetic_distances = []

    # upper triangle 
    for i in range(len(physical_matrix)):
        for j in range(i + 1, len(physical_matrix)):
            physical_distances.append(physical_matrix[i, j])
            genetic_distances.append(genetic_matrix[i, j])

    physical_distances = np.array(physical_distances).reshape(-1, 1)
    genetic_distances = np.array(genetic_distances)

    # linear regression
    regressor = LinearRegression()
    regressor.fit(physical_distances, genetic_distances)
    predicted_genetic_distances = regressor.predict(physical_distances)

    
    root_length = tree_dict['s10']
    right_internal = tree_dict['Ri']
    right_terminal = tree_dict['Rt']
    left_internal = tree_dict['Li']
    left_terminal = tree_dict['Lt']

    # root tip R
    right_root_to_tip = []
    for i in range(len(right_terminal)):
        distance = root_length + sum(right_internal[:i+1]) + right_terminal[i]
        right_root_to_tip.append(distance)

    # root tip L
    left_root_to_tip = []
    for i in range(len(left_terminal)):
        distance = root_length + sum(left_internal[:i+1]) + left_terminal[i]
        left_root_to_tip.append(distance)

    root_to_tip_distances = right_root_to_tip + left_root_to_tip

    # Height
    max_root_to_tip_distance = max(root_to_tip_distances)
    sum_root_tip = sum(root_to_tip_distances)
    
    #TOTAL length sampled in m
    sum_samp = root_length + sum(right_internal) + sum(right_terminal) + sum(left_internal) + sum(left_terminal)
    length_samp = sum_samp/10 


    # find total variants
    #variants = regressor.predict(np.array([[max_root_to_tip_distance]]))[0]
    variants = regressor.predict(np.array([[length_samp]]))[0]

    # MUTATION PER NUCLEOTIDE!!!!!!!! 
    output_mut_rate = variants / GenSize
  
    regression_equation = f"y = {regressor.coef_[0]} * x + {regressor.intercept_}"
    return variants, regression_equation, max_root_to_tip_distance, output_mut_rate



#Model definition function


# which model
def determine_model(StD, biasVar):
    if StD == NumStem and biasVar == 0.5:
        return "Model #1"
    elif StD == NumStem and biasVar == 10:
        return "Model #2"
    elif StD == 0 and biasVar == 0.5:
        return "Model #3"
    elif StD == 0 and biasVar == 10:
        return "Model #4"
    return "Unknown"


#Output Csv Function

import sys
def write_samples_to_csv(param_samples, job_id, filename="simulation_results100.csv"):
    filename = f"simulation_results_{job_id}.csv"  
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Sample_Number", "GenSize", "tree_topology", "StD", "biasVar", "Model", "AveMeristemTime", "Variants", "Regression_Equation", "Unique_Mutations", "Shared_Mutations", "Age", "Height_m", "input_mut_rate", "output_mut_rate"])
        
        for i, sample in enumerate(param_samples):
            sample_number = i + 1
            GenSize = sample["GenSize"]
            tree_topology = sample["tree_topology"]
            StD = sample["StD"]
            biasVar = sample["biasVar"]
            model = determine_model(StD, biasVar)
            mu_0 = sample["mu_0"]
            input_mut_rate = mu_0 * 10**(-10)

            tree_list, tree_dict, numBranch, age = create_tree_list_and_dict(tree_topologies_dict[tree_topology])
            ave_meristem_time = calculate_ave_meristem_time(tree_dict, numBranch)

            resultList_2 = simulate_somatic_mutations(tree_dict, NumStem, NumTime, mu_0, GenSize, StD, NumDiv, nDiv, biasVar)

            # check br_mutmatrix
            if is_br_mutmatrix_empty(resultList_2):
                print(f"no mutations detected!")
                continue  

            _, _, shared_mutations, unique_mutations, mutations_counts = mut_dist_func(resultList_2)
            genetic_matrix, physical_matrix = gen_matrices(tree_dict, unique_mutations, shared_mutations, mutations_counts)
            variants, regression_equation, max_root_to_tip_distance, output_mut_rate = calc_variants(genetic_matrix, physical_matrix, tree_dict, GenSize)

            height_cm = age * 10
            height_m = height_cm / 100

            writer.writerow([sample_number, GenSize, tree_topology, StD, biasVar, model, ave_meristem_time, variants, regression_equation, unique_mutations,shared_mutations, age, height_m, input_mut_rate, output_mut_rate])
            gc.collect()  # clean

# RUN IT :)
if __name__ == "__main__":
    job_id = sys.argv[1]
    num_samples = 100
    param_samples = latin_hypercube_sampling(num_samples)
    write_samples_to_csv(param_samples, job_id)
