# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 22:23:06 2024

@author: hp
"""

import random
import json
import os.path as osp
from collections import defaultdict
from functools import reduce
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from utils import sketch_based_greedy_RTlL, order_based_SBG_RTlL, build_upper_bound_label, calculate_candidate_edges, generate_user_groups
from utils import read_temporary_graph_data, read_graph_from_edgefile
from utils import draw_networkx, draw_evaluation
import time
import copy
import random
import heapq
from typing import List
from functools import reduce
import numpy as np
from tqdm import tqdm

import networkx as nx


T = 6
R=200

bit_operationdic={
         0: 0b00000001,
         1: 0b00000010,
         2: 0b00000100,
         3: 0b00001000,
         4: 0b00010000,
         5: 0b00100000,
         6: 0b01000000,
         7: 0b10000000}

def compute_montesimu_spread(graph: nx.Graph, users: List[int], mask: List[int] = []):
    """Compute independent cascade in the graph

    Args:
        graph (nx.Graph): graph object
        users (List[int]): a set of user nodes
        mask (List[int], optional): a set of unparticipation cascade mask nodes. default is []. 

    Returns:
        spread, reached_node (int, List[int]): the number of vertexes reached by users in graph and the set of nodes reached by users.
    """
    spread=[]
    popo =[0.1,0.05,0.01]
    for i in range(10000):
        np.random.seed(i)
        new_active, active = users[:], users[:]
            
        # for each newly activated nodes, find its neighbors that becomes activated
        while new_active:
            activated_nodes = []
            for node in new_active:
                for nodenei in list(graph.neighbors(node)):                    
                    # if np.random.uniform(0,1) < 1/graph.degree(nodenei):
                    # if np.random.uniform(0,1) < 0.1:
                    if np.random.uniform(0,1) < random.choice(popo):
                       activated_nodes.append(nodenei)                                                 
            new_active = list(set(activated_nodes) - set(active) - set(mask))
            active += new_active  
        spread.append(len(active))
    return np.mean(spread)

def subgraph_getfirst_descendants(graphs: List[nx.Graph], users: int):
    """Compute independent cascade in the graph

    Args:
        graph (nx.Graph): graph object
        users (List[int]): a set of user nodes
        mask (List[int], optional): a set of unparticipation cascade mask nodes. default is []. 

    Returns:
        spread, reached_node (int, List[int]): the number of vertexes reached by users in graph and the set of nodes reached by users.
    """
    descendants =[]
    for j in range(0,R):
        visited = set()
        if(users in graphs[j]):
            visited.add(users)
            Queue1 = []
            Queue1.append(users)
            while(len(Queue1)>0):
                nodecur=Queue1.pop(0)
                neighbors=list(graphs[j].neighbors(nodecur))
                for neighbor in neighbors:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        Queue1.append(neighbor)  
        descendants.append(visited)
    # for j in range(0,R):
    #     if users in graphs[j].nodes:
    #         tmp = list(nx.bfs_successors(graphs[j], users))
    #         # tmp = set(tmp)
    #     else:
    #         tmp = set() 
    #     descendants.append(tmp)                      
    return descendants

def subgraph_get_descendants(graphs: List[nx.Graph], users: int, K: int):
    """Compute independent cascade in the graph

    Args:
        graph (nx.Graph): graph object
        users (List[int]): a set of user nodes
        mask (List[int], optional): a set of unparticipation cascade mask nodes. default is []. 

    Returns:
        spread, reached_node (int, List[int]): the number of vertexes reached by users in graph and the set of nodes reached by users.
    """
    nodespread=[]
    descendants =[]
    
    for j in range(0,R):
        visited = set()
        if(users in graphs[j]):
             if users not in Seedsetdictionary[K][j]:
                visited.add(users)
                Queue1 = []
                Queue1.append(users)
                while(len(Queue1)>0):
                    nodecur=Queue1.pop(0)
                    neighbors=list(graphs[j].neighbors(nodecur))
                    for neighbor in neighbors:
                        if neighbor not in visited:
                            if neighbor not in Seedsetdictionary[K][j]:
                                visited.add(neighbor)
                                Queue1.append(neighbor) 
            # visited = visited - Seedsetdictionary[K][j]
        descendants.append(visited)
        nodespread.append(len(visited)) 
    # for j in range(0,R):
    #     if user in graphs[j].nodes:
    #         tmp = nx.descendants(graphs[j], user) - Seedsetdictionary[K][j]
    #     else:
    #         tmp = set() 
    #     descendants.append(tmp)
    #     nodespread.append(len(tmp))                      
    return descendants, np.mean(nodespread)



def generate_snapshots(graph: nx.Graph, r: int, seed: int = 42):
    """Generate r random sketch graph by removing each edges with probability 1-P(u,v), which defined as 1/degree(v). 

    Args:
        graph (nx.Graph): graph object
        r (int): the number of snapshots generated
        seed (int): the random seed of numpy

    Returns:
        snapshots (List[nx.Graph]): r number sketch subgraph
    """
    popo =[0.1,0.05,0.01]
    np.random.seed(seed)
    snapshots = []
    for _ in range(r):
        # select_edges = [edge for edge in graph.edges if np.random.uniform(0, 1) < 0.1]
        # select_edges = [edge for edge in graph.edges if np.random.uniform(0, 1) < 1/graph.degree(edge[1])]
        select_edges = [edge for edge in graph.edges if np.random.uniform(0, 1) < random.choice(popo)]        
        snapshots.append(graph.edge_subgraph(select_edges))

    return snapshots


def HCS_influence_compute(graphs: List[nx.Graph], users: int):
    """The Forward Influence Sketch method

    Args:
        graphs (List[nx.Graph]): a set of snapshot graph
        users (List[int]): a group of user nodes
        reconneted_edge (tuple[int, int]): the reconneted edge

    Returns:
        spread (float): the mean of additional spread of vertexes reached by users in all sketch subgraph.
    """
   # T=12
    R=200
  
    Spreadcompu=[0]*T
    
    
    global e_dic    
    
    nl_dictmp={}
    nt_dictmp={}
    
   
    for j in range(0,R):
         # for node in graphs[j].nodes:
         #     nl_dictmp[node]=n_dic[j][node]
          if (users in graphs[j]):
              # if (n_dic[j][users]):
                # for i in range(0,T):
                    # if (n_dic[j][users] & bit_operationdic[i]):
                    # Spreadcompu[i]=Spreadcompu[i]+1
                visited = set()
                visited.add(users)
                nl_dictmp[users]=0b00111111
                Queue1 = []
                Queue2 = []   
                Queue1.append(users)
                Queue2.append(nl_dictmp[users])
                while(len(Queue1)>0):
                    nodecur=Queue1.pop(0)
                    bitcur=Queue2.pop(0)
                    neighbors=list(graphs[j].neighbors(nodecur))
                    for neighbor in neighbors:
                        if neighbor not in visited:
                              nl_dictmp[neighbor]=0b00111111
                              # d = e_dic[j][(nodecur,neighbor)]
                        b=bitcur & e_dic[j][(nodecur,neighbor)] & nl_dictmp[neighbor]
                        
                        if (b):
                            for i in range(0,T):
                                if (b & bit_operationdic[i]):
                                    Spreadcompu[i]=Spreadcompu[i]+1
                            visited.add(neighbor)
                            # if
                            nl_dictmp[neighbor]=nl_dictmp[neighbor] ^ b
                            Queue1.append(neighbor)
                            Queue2.append(b)
                            
    for i in range(0,T):
        Spreadcompu[i]=Spreadcompu[i]/R+1
     
    
    return Spreadcompu

datasets = ['EmailEuCore', 'MathOverflow', 'AskUbuntu', 'StackOverflow']
Seedsize = [5,10,15,20,25]  
# datasets = ['Superuser']
# Seedsize = [20]   
t=60
for dataset in datasets:      
    snapshots = []
    subgraphs=[]
    Graphtotalnodes=set()
    e_dic=[]
    GraphHC=[]
    timelapse = []
    Totalnodes = []
    for i in range(1,T+1):
        pred_graph = read_graph_from_edgefile(f'data/SEALDataset/{dataset}/T{t}_pred_edge{i}.pt')
        Graphtotalnodes=Graphtotalnodes.union(set(pred_graph.nodes))
        snapshots.append(pred_graph)
    Totalnodes = list(Graphtotalnodes)
    for w in range(0,T):
        tmpsubgraph = generate_snapshots(snapshots[w], R, 42)
        subgraphs.append(tmpsubgraph)
    
    for j in range(0,R):    
        GraphHCtmp=nx.DiGraph()
        e_dictmp={}
        n_dictmp={}
        for i in range(0,T):
            for edge in subgraphs[i][j].edges:  
                if(GraphHCtmp.has_edge(*edge)):                
                    e_dictmp[edge]=e_dictmp[edge]|bit_operationdic[i]
                else:
                    GraphHCtmp.add_edge(*edge)
                    e_dictmp[edge]=bit_operationdic[i]
        GraphHC.append(GraphHCtmp)
        e_dic.append(e_dictmp)
    for SZ in Seedsize:     
        start_time = time.time() 
        SEED=[]
        TOTALSPREAD=[]
        Spreadtotal=[HCS_influence_compute(GraphHC,node) for node in Totalnodes]
        upperbound = list(zip(*Spreadtotal))    
    
        for snapshotindex in range(0,T):    
            Seedsetdictionary = []
            # start_time = time.time()
            Q = sorted(zip(Totalnodes,upperbound[snapshotindex]), key=lambda x: x[1],reverse=True)
            S, spread, SPREAD = [Q[0][0]], Q[0][1], [Q[0][1]] 
            # teststarttime = time.time()
            Seedsetdictionary.append(subgraph_getfirst_descendants(subgraphs[snapshotindex],Q[0][0]))  
        
            Q= Q[1:]
            for i in range(1, SZ):    
        
              check = False
                
              while not check:
        
                    # node_lookup += 1
        
                    current = Q[0][0]
        
                    [tmpchildren,tmpspread] = subgraph_get_descendants(subgraphs[snapshotindex], current, i-1)
                    Q[0] = (current,tmpspread)
        
                    Q = sorted(Q, key = lambda x: x[1], reverse = True)
        
                    check = (Q[0][0] == current)
        
              Seedsetdictionary.append(tmpchildren)
              for j in range(0,R):
                    Seedsetdictionary[i][j] = Seedsetdictionary[i][j].union(Seedsetdictionary[i-1][j])
              spread += Q[0][1]
              S.append(Q[0][0])
              # SPREAD.append(spread)
              # LOOKUPS.append(node_lookup)
              Q = Q[1:] 
            SEED.append(S)
            TOTALSPREAD.append(spread)
            # LOOKUPNUM.append(LOOKUPS)
        timelapse = time.time() - start_time
        evaspread=[]
        for i in range(0,T):
            pred_graph = read_graph_from_edgefile(f'data/SEALDataset/{dataset}/T{t}_pred_edge{i+1}.pt')
            evaspread.append(compute_montesimu_spread(pred_graph, SEED[i]))
        # data={'Seedsize':SZ,'Runtime':timelapse,'Totalspread0':evaspread[0],'Totalspread1':evaspread[1],'Totalspread2':evaspread[2],'Totalspread3':evaspread[3],'Totalspread4':evaspread[4]}
        # data={'Seedsize':SZ,'Runtime':timelapse,'Totalspread0':evaspread[0],'Totalspread1':evaspread[1],'Totalspread2':evaspread[2],'Totalspread3':evaspread[3],'Totalspread4':evaspread[4],'Totalspread5':evaspread[5]}
        # data={'Seedsize':SZ,'Runtime':timelapse,'Totalspread0':evaspread[0],'Totalspread1':evaspread[1],'Totalspread2':evaspread[2],'Totalspread3':evaspread[3]}
        # data={'Seedsize':SZ,'Runtime':timelapse,'Totalspread0':evaspread[0],'Totalspread1':evaspread[1],'Totalspread2':evaspread[2]}
        # data={'Seedsize':SZ,'Runtime':timelapse,'Totalspread0':evaspread[0],'Totalspread1':evaspread[1]}
        data={'Seedsize':SZ,'Runtime':timelapse,'Totalspread0':evaspread[0],'Totalspread1':evaspread[1],'Totalspread2':evaspread[2],'Totalspread3':evaspread[3],'Totalspread4':evaspread[4],'Totalspread5':evaspread[5]}
        df = pd.DataFrame(data, index=[0])    
        df.to_csv(f'CSINT{dataset}SZ{SZ}N{T}Random.csv', index=False) 
