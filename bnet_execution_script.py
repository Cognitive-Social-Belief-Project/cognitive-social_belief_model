##------------------------------------------------------------------------------
##
##    Imports
##
##------------------------------------------------------------------------------
import matplotlib
from networkx.classes.function import subgraph
# from sympy.functions.elementary.complexes import re
matplotlib.use('Agg') ### remember to comment this back in#####
import networkx as nx
import copy
import random as rnd
import scipy.stats as stats
import scipy.spatial.distance as scydist
import scipy.optimize as opt
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import bisect
import sys
import pickle
import operator
# from sklearn.cluster import DBSCAN
import math
import scipy.cluster.hierarchy as hier
import os
import subprocess as sub
from collections import Counter
from collections import defaultdict
from audiodev import test
from multiprocessing import Process, Queue
from datetime import datetime
import itertools
import shutil
import scipy.ndimage
# from mpmath import mp
# import heapq
# import randht
# import plvar
# import plfit
# import plplot
# import plpva

##------------------------------------------------------------------------------
##
##    Classes
##
##------------------------------------------------------------------------------

                
                
##------------------------------------------------------------------------------
##
##    Functions
##
##------------------------------------------------------------------------------
def Normalize(points, minval=None, maxval=None):

    if minval==None:
        minval = np.min(points)
    if maxval==None:
        maxval = np.max(points)
    points = np.array(points)

    return (points.astype(float) - minval) / (maxval - minval)

def relabel_nodes(nxGraph):
    """
    """
    
    mapping=dict(zip(nxGraph.nodes(),range(0,len(nxGraph))))
    relabelled_graph = nx.relabel_nodes(nxGraph,mapping)
    
    return relabelled_graph

def save_object(obj, filename):
    """
    Save an object to file for later use.
    """
    
    file = open(filename, 'wb')
    pickle.dump(obj, file)
    file.close()

def load_object(filename):
    """
    Load a previously saved object.
    """
    
    file = open(filename, 'rb')
    return pickle.load(file)

def cycle_search(graph, start, target, result, cycles):
    """
    A depth first cycle search for all cycles that pass through target
    Changes the cycles list to produce output
    Cycles list is empty if none are found
    
    graph
        uses a networkx graph object
    """
    if target == start:
        cycles.append(copy.copy(result + [start]))
        
    if start in result:
        return False
    
    result.append(start)
    
    for neighbor in self.graph.neighbors(start):
        cycle_search(neighbor, target, result, cycles)
         
    # No cycle was found
    result.pop()
    return False

def list_triangles(G):
    """
    G
        a networkx graph object
        
    Returns a list of triangles for the assoc_net
    """
    
    # Get a sorted list of nodes from greatest to least
    sorted_node_list, sorted_degrees = zip(*dict_sorter(G.degree(), True))
    triangles = []
    mem = [ [] for i in G.nodes()]
    for i, node in enumerate(sorted_node_list): #low i is high degree
        # Loop over neighbors of node
        for j in G.neighbors(node):
            # Only run if the edge has not been checked
            if sorted_node_list.index(j) > i:
                for w in list(set(mem[j]).intersection(mem[node])):
                    triangles.append((w,j,node))          

                mem[j].append(node)
    
    return triangles

def dict_sorter(x, reverse=False):
    """
    Return a list of tuple where first element is node and second is value.
    Sorts dict based on value not key, keeps key in tuple.
    """
    
    # return copy.deepcopy(sorted(x.iteritems(), key=operator.itemgetter(1), \
    #                             reverse=reverse))
    return sorted(x.iteritems(), key=operator.itemgetter(1), reverse=reverse)
    
def boltzmann_prob(delta_coherence, I):
    """
    """
    
    if I != 0.0:
        return np.exp(delta_coherence / I)
    else:
        return 0

def boltzmann_prob_mem(delta_coherence, memory, I):
    """
    """
    
    return np.exp(delta_coherence / memory / I)

def weighted_choice(nodes, weights):
    """
    Returns a random weighted selection
    """
    
    assert(len(nodes) == len(weights))
    cumdist = list(accumulate(weights))
    pick = rnd.random() * cumdist[-1]
    return nodes[bisect.bisect(cumdist, pick)]

def accumulate(iterable):
    """
    Return running totals -- P3 -> P2.7
    accumulate([1,2,3,4,5]) --> 1 3 6 10 15
    -Bumbed off of python site cus 2.7 doesn't have accumulate or operator
    """

    it = iter(iterable)
    total = next(it)
    yield total
    for element in it:
        total += element
        yield total
        
def two_cliq_a2a(N, Xo):
    """
    Creates a two clique adjacency matrix where each clique is all-to-all
    connected, but there are no connections between the two cliques.
    
    N
        set the number of nodes
    Xo
        set the proportion of nodes in each clique: float [0,1]
        
    Returns an NxN numpy array where the first floor(N*Xo) i columns represent 
    the first clique.
    
    """
    
    # Create adjacency matrix
    adj_matrix = np.zeros((N,N))
    adj_size = np.shape(adj_matrix)
    
    # Determine structure of matrix
    c1_range = range(int(math.ceil(N*Xo)))
    c2_range = range(max(c1_range)+1,N)     
    
    # Assign connections
    for row in range(adj_size[0]):
        for col in range(adj_size[1]):
            if (row != col) and ( ((row in c1_range) and (col in c1_range)) \
                or ( (row in c2_range) and (col in c2_range))):
                adj_matrix[row][col] = 1
    
    return adj_matrix

def connect_cliques(adj_matrix, connection_prob):
    """
    This function connects the nodes of each either clique. Each node has a 
    connection probability to build a connection to a node of another clique.
    
    adj_matrix
        the adjacency matrix to be operated on. Must be a square, symmetric 
        matrix.
    connection_prob
        the probability of each node to build a connection with another node
        
    Returns nothing. The adj_matrix is modified in place.
    
    """
    
    matrix_size = np.shape(adj_matrix)
    for i in range(matrix_size[0]):
        for j in range(matrix_size[1]):
            if (i != j) and (adj_matrix[i][j] != 1) and (adj_matrix[j][i] != 1):
                rand_val = rnd.random()
                if rand_val <= connection_prob:
                    adj_matrix[i][j] = 1
                    adj_matrix[j][i] = 1

def log_linspace(start, stop, number):
    """
    """
    
    lstart = np.log(start)
    lstop = np.log(stop)
    ldata = np.linspace(lstart, lstop, number)
    
    return np.exp(ldata)

def readcol(filename,format='s',delimiter=None,col=None,skiplines=[]):
    """
    readcol(filename,format='s',delimiter=None,col=None,skiplines=[])

    Reads columns from a file.
    
    filename 
        filename of input file
    format 
        string of format characters (auto formats output)
        s - string
        i - int
        f - float
        Letters must be chosen in the format of a string for each
        column to be output. The default setting for format is 's'.
    delimiter 
        char used to separate columns in a line of text. If None is
        chosen then the full line of a file will be read in.
    col 
        can optionally list the columns you want to return (list of ints).
        By default col=None which means that all columns will be read if
        a delimiter is chosen, or only one if no delimiter is chosen.
    skiplines 
        list of lines to skip in reading.
    
    Returns lists. When multiple columns are chosen the function returns
    each one in the order provided in the col argument. If no col argument was
    given, yet all the columns were read, then a single list containing all the 
    columns is returned.
    
        example: firstcol, secondcol, thrdcol = readcol(file,'sif',',',[1,2,3])
    or
        example cols = readcol(file,'s',',')
    
    """

    # Reject bad input
    if  (delimiter == None) and (col != None): 
        if (len(format) == len(col)):
            sys.exit("Must have delimiter for multiple columns.")

    # Open file and read lines into variables.
    if os.path.exists(filename):
        infile = open(filename,'r')
    else:
        sys.exit("The file named: " + filename + " does not exist.")
    
    # Default read all columns as a single format
    # This requires a delimiter to be defined, col set to None
    if (delimiter != None) and (col == None) and (len(format) == 1):
        data = []
        for i,line in enumerate(infile.readlines()):
            if (i+1 in skiplines) == False: 
                fields = line.strip().split(delimiter)
                fields = [field for field in fields if field != ""]
                row = []
                for str in fields:
                    row.append(str)
                data.append(row)
        # Put row values into their respective columns
        columns = zip(*data)
        
        if len(columns) == 1:   
            # Format data 
            col_list = list(columns[0])
            type = format[0]
            if type == 'i':
                col_list = [int(val) for val in col_list]
            elif type == 'f':
                col_list = [int(val) for val in col_list]
            elif type != 's':
                sys.exit("Warning: Unrecognized type " + type + " chosen!")    
            return col_list
        else:
            # Format data
            col_list = [list(tupl) for tupl in columns]
            type = format[0]
            if type == 'i':
                col_list = [map(int,col) for col in col_list]
            elif type == 'f':
                col_list = [map(float,col) for col in col_list]
            elif type != 's':
                sys.exit("Warning: Unrecognized type " + type + " chosen!")    
            return col_list 
        
    # Read a single column file
    # This requires the delimiter to be set to None, and col to be set to None
    # Only the first format character is used for formating, all others are 
    # ignored
    elif (delimiter == None) and (col == None):
        data = []
        for i,line in enumerate(infile.readlines()):
            if (i+1 in skiplines) == False: 
                field = line.strip()
                data.append(field)
                
        # Format data 
        type = format[0]
        if type == 'i':
            data = map(int,data)
        elif type == 'f':
            data = map(float,data)
        elif type != 's':
            sys.exit("Warning: Unrecognized type " + type + " chosen!")    
        return data
    
    # Read multicolumn file with different formats for the first N columns
    # where N is the number of format options chosen
    # This requires a delimiter to be set and col to be set to None
    elif (delimiter != None) and (col == None) and (len(format) > 1):
        data = [[] for dummy in range(len(format))]
        for i,line in enumerate(infile.readlines()):
            if (i+1 in skiplines) == False: 
                fields = line.strip().split(delimiter)
                fields = [field for field in fields if field != ""]
                for j,val in enumerate(fields):
                    try:
                        data[j].append(val)
                    except IndexError:
                        pass
    
    # Read multicolumn file with different formats for specific columns
    # This requires a delimiter to be set and col to be set
    # Col must be the same length as format 
    # *there is no difference between reading a single column with col set
    #  and having default col with no delimiter set 
    elif (delimiter != None) and (col != None):
        data = [[] for dummy in range(len(format))]
        assert(len(col) == len(format))
        for i,line in enumerate(infile.readlines()):
            if (i+1 in skiplines) == False: 
                fields = line.strip().split(delimiter)
                fields = [field for field in fields if field != ""]
                for j in range(len(col)):
                    try:
                        data[j].append(fields[col[j]-1])
                    except IndexError:
                        pass
    
    # Read 
    else:
        sys.exit("Error: Inappropriate input provide.")
                    
    infile.close()
        
    # Format data
    for i,type in enumerate(format):
        if type == 'i':
            data[i] = map(int,data[i])
        elif type == 'f':
            data[i] = map(float,data[i])
        elif type != 's':
            sys.exit("Warning: Unrecognized type " + type + " chosen!")
    
    # Convert to tuple 
    if len(data) == 1:
        columns = tuple(data[0])
    else: 
        columns = tuple(data)
    
    return columns
        
def read_LFR_Benchmark(edge_file, community_file):
    """
    Reads a LFR style output file into a networkx graph object
    """

    col1, col2 = readcol(edge_file,'ii',delimiter='\t')
    col1 = np.array(col1) - 1 # Set first node ID to 0
    col2 = np.array(col2) - 1 # Set first node ID to 0
    edge_list = zip(col1, col2)
    nodes, clusters = readcol(community_file, 'ii',delimiter='\t')
    nodes = np.array(nodes) - 1 # Set first node ID to 0
    community_list = zip(nodes, clusters)
    
    # Create graph
    LFR_graph = nx.Graph(edge_list)
    
    # Assigne community values
    for node in LFR_graph:
        LFR_graph.node[node]['community'] = community_list[node][1]
    
    return LFR_graph

def get_nodes_list(graph, community, community_key):
    """
    """

    return [ node for node in graph.nodes() if graph.node[node][community_key]==community ]

def relabel_graph(graph, assignment_dict):
    """
    """
    
    # Get dictionary of key-communities values-nodes
    # Get community list
    communities = []
    for node in graph.nodes():
        communities.append(graph.node[node]['community'])

    # Remove copies
    communities = set(communities)

    # Community nodes lists
    community_node_lists = {}
    for community in communities:
        community_node_lists[community] = get_nodes_list(graph, community, 'community')
        
    # Loop through zealot communities and great new labels
    new_labels = {}
    node_id = 0
    for key in assignment_dict.keys():
        if assignment_dict[key][0] == 'z':
            for node in community_node_lists[key]:
                new_labels[node] = node_id
                node_id += 1
    # Loop through non-zealots and add to labels
    for key in assignment_dict.keys():
        if assignment_dict[key][0] != 'z':
            for node in community_node_lists[key]:
                new_labels[node] = node_id
                node_id += 1
                
    # Do reassignment
    return nx.relabel_nodes(graph, new_labels, copy=True)

def community_assignment(graph, assignment_dict):
    """
    assignment_dict has the following structure
    
    key = community in qoutes: '1', '2', '3', ... etc
    value = first char is 'b' or 'z' for belief or zealot
            second char is index # in set: 'b1', 'z2', ... etc
    """
    
    assignment_list = []
    for node in graph.nodes():
        value = assignment_dict[graph.node[node]['community']]
        type = value[0]
        subset = value[1:]
        assignment_list.append((node, type, subset))
        
    return assignment_list
    
def write_LFR_param_file(size, avg_deg, max_deg, exp_deg, exp_com, mix, \
                       min_com=None, max_com=None):
    """
    """
    
    lines = []
    lines.append(str(size) + '\n')
    lines.append(str(avg_deg) + '\n')
    lines.append(str(max_deg) + '\n')
    lines.append(str(exp_deg) + '\n')
    lines.append(str(exp_com) + '\n')
    lines.append(str(mix) + '\n')
    if min_com != None:
        lines.append(str(min_com) + '\n')
    if max_com != None:
        lines.append(str(max_com) + '\n')
        
    outfile = open('parameters.dat', 'w')
    outfile.writelines(lines)
    outfile.close()

def generate_LFR_graph(size, avg_deg, max_deg, exp_deg, exp_com, mix, \
                       min_com=None, max_com=None):
    """
    """
    
    # Create parameter file
    write_LFR_param_file(size, avg_deg, max_deg, exp_deg, exp_com, mix,\
                         min_com, max_com)
    
    # Run program
    os.popen("./benchmark")
    
    # Read output
    LFR_graph = read_LFR_Benchmark('network.dat', 'community.dat')
    
    return LFR_graph
        
def fit_power_law(filename, cluster_size_file, ilist_file, piter=0):
    """
    cutoff is a value at which the data is cutoff it is
    exceeds (for estimation and fitting purposes).
     
    """
    
    cluster_sizes = load_object(cluster_size_file)
    ilist = load_object(ilist_file)
    
    pvalues = []
    alphas = []
    xmins = []
    Ls = []
    # Loop through each I
    for i, ival in enumerate(ilist):
        pvalue = None
        
        if len(cluster_sizes[i]) < 120:
            alpha, xmin, L = plfit.plfit(cluster_sizes[i], 'finite', 'range',[1.01,4.00,0.01], 'nowarn', True)
        else:
            alpha, xmin, L = plfit.plfit(cluster_sizes[i], 'nosmall','range',[1.01,4.00,0.01], 'nowarn', True)
        
        if alpha != 'Not a Number':
            # P-value of optimal parameters
            if piter != 0:
                pvalue, gof = plpva.plpva(cluster_sizes[i],xmin,'silent','reps', piter)
        
            if piter != 0:
                pvalues.append(pvalue)
        
            alphas.append(alpha)
            xmins.append(xmin)
            Ls.append(L)
            plot_powerlaw_fit(filename, cluster_sizes[i], alpha, xmin, ival, pvalue)
        
    return alphas, xmins, Ls, pvalues

def plot_powerlaw_fit(filename, data, alpha, xmin, ival, pvalue=None):
    """
    """
    
    cdf, bins, patches = plt.hist(data, np.logspace(np.log10(min(data)), \
                            np.log10(max(data))), normed=True, cumulative=-1)                             
    xbins = (bins[:-1] + bins[1:]) / 2.
    mcdf = norm_powlaw_cdf(xbins, alpha, xmin)
    
    plt.clf()
    plt.loglog(xbins, cdf, linestyle='None', marker='o', color='blue', label='Model')
    plt.loglog(xbins, mcdf, linestyle='-', marker='None', color='red', label='Best fit')
    plt.xlabel(r'$S_{g} / N$')
    plt.ylabel(r'$P(X \leq x)$')
    if pvalue == None:
        plt.title(r'$I=$' + str(round(ival,3)) + r' $\alpha=$' + str(round(alpha,2)) + \
                  r' $x_{min}=$' + str(round(xmin,2)) + r' $p=$' + 'na')
    else:
        plt.title( r'$I=$' + str(round(ival,3)) + r' $\alpha=$' + \
                  str(round(alpha,2)) + r' $x_{min}=$' + str(round(xmin,2)) + r' $p=$' + str(round(pvalue,2)))
    plt.legend()
    plt.savefig(filename + '-' + str(ival) + '-powlaw-fits.pdf', dpi=1200,bbox_inches='tight')
    plt.clf()
    plt.close()
        
def cont_powlaw_gen(n, alpha, xmin):
    """
    """
    
    # Does not work for alpha less than 1.0
    assert(alpha >= 1.0)
    
    return xmin * np.power(1 - stats.uniform.rvs(loc=0, scale=1, size=n), -1. / (alpha - 1.))
        
def ks_stat(xmin, data):
    """
    """
    
    alpha = cont_powlaw_mle(data, xmin)
    
    # Remove data
    cut_data = [ val for val in data if val > xmin ]
    if len(cut_data) <= 3:
        return 100000
    # cut_data = data # removed cut data stuff
    
    # Calculate KS statistic
    D, p = stats.kstest(cut_data, powlaw_cdf, args=(alpha, xmin))
    
    return D

def powlaw_cdf(x, alpha=1, xmin=1):
    """
    """

    return np.power(x / xmin, -alpha + 1)

def norm_powlaw_cdf(x, alpha=1, xmin=1):
    """
    """
    
    xo = powlaw_cdf(x[0], alpha, xmin)
    return powlaw_cdf(x, alpha, xmin) / xo
    
def cont_powlaw_mle(data, xmin=1):
    """
    """
    
    cut_data = [ val for val in data if val >= xmin ]
    
    return 1 + len(data) * np.power(np.sum(np.log(np.array(cut_data) / float(xmin))), -1)

def discrete_powlaw_mle(data, xmin=1.):
    """
    """
    
    cut_data = [ val for val in data if val >= xmin ]
    
    return 1 + len(data) * np.power(np.sum(np.log(np.array(cut_data) / (xmin - 0.5))), -1)

def alpha_mle_std_err(est, data):
    """
    """
    
    return (est - 1) / np.sqrt(len(data))
    
def nested_run_completion(current, max, pid=1):
    """
    
    current
        list of current values in iteration from inner most to outer
    
    max
        list of maximum values for iteration from inner most to outer
    
    """
    
    return np.sum([ val * np.prod([ max[j] for j in xrange(i) ]) \
                   for i, val in enumerate(current) ]) / np.prod(max)
                   
def update_progress(progress):
    """
    Progress bar
    """
    
    sys.stdout.write('\r') 
    sys.stdout.write("[%-50s] %d%%" % ('='*int(progress*100 / 2.), progress*100))
    sys.stdout.flush()

def mean_and_error(list_of_domains):
    """
    """
    
    means = np.array([ np.mean(domain) for domain in list_of_domains ])
    std_err = np.array([ np.std(domain) / np.sqrt(float(len(domain))) * 1.96 for domain in list_of_domains ])
    
    return means, std_err

def pre_cache_coherence(num_nodes, filename='coh_file.dat', formulation='strong'):
    """
    caches all the association combinations as a dictionary whose keys are 
    the combination (in vector form)
    
    """

    num_assoc = (num_nodes * (num_nodes - 1)) / 2
    combinations = itertools.product([-1,1],repeat=num_assoc)
    sequence = [ (combo, default_map(combo, num_nodes, num_assoc, formulation))\
                 for combo in combinations ]
    cached_coh = { key: value for (key,value) in sequence }
    translator = { key: value for (key,value) in zip(range(num_assoc), nx.fast_gnp_random_graph(num_nodes, 1.0).edges()) }
    rev_translator = { key: value for (key,value) in zip(nx.fast_gnp_random_graph(num_nodes, 1.0).edges(), range(num_assoc) ) }
    
    save_object((cached_coh, translator, rev_translator), filename)
    
def default_map(assoc, num_nodes, num_assoc, formulation='strong'):
    """
    """
    
    # Initialize graph
    frame = nx.fast_gnp_random_graph(num_nodes, 1.0)
    for i, edge in enumerate(frame.edges()):
        frame[edge[0]][edge[1]]['tie'] = assoc[i]
    
    # Determine coherence
    return calc_avg_coherence(frame, formulation)

def calc_avg_coherence(G, formulation='strong'):
    """
    Calculate the coherence (-energy) of graph G.
    
    G
        a networkx object with a 'tie' property that is either True/False
        if G has weights, then they will be used to determine energy
        
    formulation
        determines whether the 'strong' or 'weak' interpretation of ties
        is used. If strong, then --- configuration is unstable. If weak,
        then --- configuration is stable. 
    
    """

    # Get a list of triangles for the graph
    triangles = list_triangles(G)

    # Compute the coherence contribution for each triangle
    coherence = 0.0
    if formulation == 'strong':
        for triangle in triangles:
            product = 1.0
            product *= G[triangle[0]][triangle[1]]['tie']
            product *= G[triangle[1]][triangle[2]]['tie']
            product *= G[triangle[0]][triangle[2]]['tie']
            coherence += product
            
    elif formulation == 'weak':
        for triangle in triangles:
            product = 1.0

            # Check if triangle has --- configuration
            if G[triangle[0]][triangle[1]]['tie'] == \
            G[triangle[1]][triangle[2]]['tie'] == \
            G[triangle[0]][triangle[2]]['tie'] == -1:
                product *= 1.
                product *= 1.
                product *= 1.
            else:
                product *= G[triangle[0]][triangle[1]]['tie']
                product *= G[triangle[1]][triangle[2]]['tie']
                product *= G[triangle[0]][triangle[2]]['tie']
            coherence += product

    # Return average coherence contribution per triangle
    if len(triangles) == 0:
        return -1.0 # Assert that a triangleless network is unstable
    return coherence / float(len(triangles))

def portion_list(nHigh, nLow, ratioHigh, ratioLow):
    """
    """
    
    assert(ratioHigh + ratioLow == 1.)
    
    total = nHigh + nLow
    ratio_list = []
    for i in xrange(nHigh):
        ratio_list.append(1. / nHigh * ratioHigh)
        
    for i in xrange(nLow):
        ratio_list.append(1. / nLow * ratioLow)
    
    return ratio_list

def writecol(filename,delimiter,firstline,*args):
    """
    writecol(filename,delimiter,firstline,*args)

    Writes to file a list of columns provided as arguments to the function.
    If input is provided for firstline that is not "", then that string
    is made the first line of the file. Columns must be of same length 
    
    filename
        file to be read
    delimiter
        character or characters used as a delimiter between columns
    firstline
        header for file, if set to '', then none is written
    *args
        lists of columns to be written to the text file
    
    """
    
    col = [arg for arg in args]
    
    # Make sure columns are of the same length
    for x in col:
        assert(len(x) == len(col[0]))
    
    lines = []
    if firstline != "":
        lines.append(firstline + '\n')
    
    col_num = range(0,len(col))
    end = col_num[-1]
    for i in range(0,len(col[0])):
        line = ''
        for j in col_num:
            if j == end:
                line += str(col[j][i])
            else:
                line += str(col[j][i]) + delimiter
        line += '\n' 
        lines.append(line)
    
    outfile = open(filename,'w')
    outfile.writelines(lines)
    outfile.close()

def write_csv_cohfiles(coherence_filelist):
    """
    """
    
    # !!! This converts to an edge based representation which is more informative
    # than an object based representation.
    # !!! I have also converted to energies. I am not doing things in 
    # coherence anymore.
    for coherence_file in coherence_filelist:
        coherence_dict, translator, rev_translator = load_object(coherence_file)
        num_associations = max(translator.keys()) + 1
        
        # Write data to coherence file
        edge_columns = zip(*coherence_dict.keys())
        energy_column = coherence_dict.values()
        # !!! Convert coherence to energy: !!!
        energy_column = [ -1.0*energy for energy in energy_column]
        write_cols = edge_columns + [energy_column]
        writecol('coh' + str(num_associations) + '_energy.cdat', ',', str(num_associations), *write_cols)
        
        # # Write data to translator file
        # assoc_column = translator.keys()
        # obj_columns = zip(*translator.values())
        # write_cols = [assoc_column] + obj_columns
        # writecol('coh' + str(num_associations) + '_trans.cdat', ',', '', *write_cols)
        
    print 'Translated cohfiles.'
    
def write_voter_energyfile(num_nodes = 5):
    """
    """
    
    num_assoc = (num_nodes * (num_nodes - 1)) / 2
    combinations = list(itertools.product([-1,1],repeat=num_assoc))
    edge_columns = zip(*combinations)

    energy_column = [ 0.0 for combo in combinations]
    write_cols = edge_columns + [energy_column]
    writecol('voter_' + str(num_assoc) + '_energy.cdat', ',', str(num_assoc), *write_cols)
    
def write_cparam_file(filename, model_parameters):
    """
    """
    
    # Text to be added to final file:
    lines = ''
    
    # General parameters; must have all of these
    out_line = str(model_parameters['outfile']) + '\n'
    sizes_line = str(model_parameters['sizesfile']) + '\n'
    oenergy_line = str(model_parameters['out_energy_file']) + '\n'
    I_line = str(model_parameters['I']) + '\n'
    T_line = str(model_parameters['T']) + '\n'
    J_line = str(model_parameters['J']) + '\n'
    runs_line = str(model_parameters['runs']) + '\n'
    coh_line = model_parameters['energy_file'] + '\n'
    update_line = model_parameters['update'] + '\n'
    points_line = str(model_parameters['npoints']) + '\n'
    iter_line = str(model_parameters['niterations']) + '\n'
    graph_line = model_parameters['graph'] + '\n'
    random_line = str(int(model_parameters['random'])) + '\n'
    zealot_line = str(int(model_parameters['zealots'])) + '\n'
    sets_line = str(int(model_parameters['num_sets'])) + '\n'
    lines += out_line + sizes_line + oenergy_line + I_line + \
                T_line + J_line + points_line + iter_line + runs_line + coh_line + update_line + \
                graph_line + sets_line + random_line
    
    # Add lines if initial belief-sys config is not random
    if model_parameters['random'] == False:
        # Add initial belief settings
        num_vectors_line = str(model_parameters['num_belief_types']) + '\n'
        belief_types_lines = ''
        for set in xrange(model_parameters['num_sets']):
            for beliefs in model_parameters['belief_types'][set]:
                belief_types_lines += ','.join(map(str, beliefs))
                belief_types_lines += '\n'

        belief_ratios_line = ','.join(map(str,model_parameters['belief_ratios'])) + '\n'
        
        lines += num_vectors_line + belief_types_lines + belief_ratios_line
    
    lines += zealot_line
    # Add lines if zealots are active
    if model_parameters['zealots'] == True:
        # Add initial zealot settings
        num_zealot_types_line = str(model_parameters['num_zealot_types']) + '\n'
        zealot_types_lines = ''
        for set in xrange(model_parameters['num_sets']):
            for zealots in model_parameters['zealot_types'][set]:
                zealot_types_lines += ','.join(map(str,zealots))
                zealot_types_lines += '\n'
            
        num_zealots_line = ','.join(map(str,model_parameters['zealot_ratios'])) + '\n'
        
        lines += num_zealot_types_line + zealot_types_lines + num_zealots_line
    
    # write graph size line
    lines += str(model_parameters['N']) + '\n'
    
    # write seed to file
    if model_parameters['seed'] == None:
        lines += '0\n'
    else:
        lines += '1\n'
        lines += str(model_parameters['seed']) + '\n'
    
    # Write assignment line
    if model_parameters['assignment']:
        lines += '1\n'
    else:
        lines += '0\n'
        
    if model_parameters['assignment'] == True:
        assignment_lines = ''
        for node in model_parameters['assignments']:
            assignment_lines += str(node[0]) + ',' + str(node[1]) + ',' + str(node[2]) + '\n'

        lines += assignment_lines

    # Write data lines
    if model_parameters['timeseries']:
        lines += '1\n'
    else:
        lines += '0\n'

    if model_parameters['finalstate']:
        lines += '1\n'
    else:
        lines += '0\n'

    if model_parameters['ts_sizes']:
        lines += '1\n'
    else:
        lines += '0\n'

    if model_parameters['AppendOut']:
        lines += '1\n'
    else:
        lines += '0\n'

    outfile = open(filename,'w')
    outfile.writelines(lines)
    outfile.close()

    print 'Parameter file: Done'
    
def write_graph_file(filename, nx_graph):
    """
    """

    num_nodes = nx.number_of_nodes(nx_graph)
    num_edges = nx.number_of_edges(nx_graph)
    lines = str(num_nodes) + ',' + str(num_edges) +'\n'
    for edge in nx_graph.edges():
        lines += str(edge[0]) + ',' + str(edge[1]) + '\n'
    
    outfile = open(filename,'w')
    outfile.writelines(lines)
    outfile.close()
    
    print 'Graph File: Done'
    
def generate_beliefsys_sets(coh_file,zealots, 
                            zealot_energy_ranges, belief_energy_ranges, num_sets):
    """
    """
    
    # Read in coherence values
    coherence_dict, translator, rev_translator = load_object(coh_file)

    belief_options = tuple([ tuple([ tuple(vect) for vect, val in coherence_dict.items()
     if (-1*val)>=belief_energy_ranges[i][0] 
     and (-1*val)<=belief_energy_ranges[i][1]]) for i in xrange(len(belief_energy_ranges))])
    
    if zealots == True:
        zealot_options = tuple([ tuple([ tuple(vect) for vect, val in coherence_dict.items()
                            if (-1*val)>=zealot_energy_ranges[i][0]
                            and (-1*val)<=zealot_energy_ranges[i][1]])
                          for i in xrange(len(zealot_energy_ranges))])
    end_counter = 0
    set_counter = 0
    zealot_system_sets = []
    belief_system_sets = []
    while (end_counter <= 100 and set_counter < num_sets):
        belief_selection = []
        for boption in belief_options: 
            belief_selection.append(rnd.choice(boption))

        zealot_selection = []
        if zealots == True:
            # If immutables are on, makes zealots
            for zoption in zealot_options:
                zealot_selection.append(rnd.choice(zoption))
            # If zealots are not duplicates, then check that they don't conflict with init-beliefs
            if len(zealot_selection) == len(set(tuple(zealot_selection))):
                conflict = False
                for belief in belief_selection:
                    # Adds to end-counter and retarts selection if a zealot and initial belief are the same
                    if belief in zealot_selection:
                        end_counter += 1
                        conflict = True
                        break
                if conflict == True:
                    continue
            
            else: # dont allow duplicate zealots
                end_counter += 1
                continue
                    
        # dont allow duplicate beliefs
        if len(belief_selection) != len(set(belief_selection)):
            end_counter += 1
            continue

        # Checks if init-beliefs are in one of the existing sets already
        if belief_selection in belief_system_sets:

            # If immutablility is on, we allow for the dublicate. But have to do additional checks
            if zealots == True:
                
                # CHeck if zealots have already been chosen in a set
                if zealot_selection in zealot_system_sets:

                    # If both the zealots and init beliefs found in the sets are both in the same set, we reset loop
                    for i in xrange(len(belief_system_sets)):
                        conflict = False
                        if (belief_selection == belief_system_sets[i]) and (zealot_selection == zealot_system_sets[i]):
                            end_counter += 1
                            conflict = True
                            break
                    
                    if conflict == True:
                        continue
                    # If they are in different sets, we continue and add them as a new set
                    else:
                        belief_system_sets.append(belief_selection)
                        set_counter += 1
                        end_counter = 0
                        zealot_system_sets.append(zealot_selection)
                        continue
                    
                # if the zealots are new, then we can add this configuration of zealots and init-beliefs as a set
                else:
                    belief_system_sets.append(belief_selection)
                    set_counter += 1
                    end_counter = 0
                    zealot_system_sets.append(zealot_selection)
                    continue
            
            # If immutability is off, we reset selection and add to end-counter
            end_counter += 1
            continue
        
        # if init-belief config is not in a set already, and no zealots overlap, then we can add it on the set
        else:
            belief_system_sets.append(belief_selection)
            set_counter += 1
            end_counter = 0
            zealot_system_sets.append(zealot_selection)
            continue


    return belief_system_sets, zealot_system_sets, len(belief_system_sets)

def run_cpp(filename, param_file):
    """
    """
    
    command = './' + filename
    
    try:
        sub.call([command, param_file])
    except:
        print "Failure in, or inability to run: " + filename
        sys.exit()
        
def read_sizes_output(filename, num_points):
    """
    Read the sizes data file from the cpp program
    """
    
    try:
        sets = [ [ ] ]
        with open(filename) as inFile:
            point_coint = 0
            for line in inFile:

                line_list = line.strip().split(',')
                if "NEWPOINT" == line_list[1]:
                    if point_coint == num_points:
                        
                        sets.append([[]])
                        point_coint = 1
                    else:
                        point_coint += 1
                        sets[-1].append([])
                else:
                    belief = line_list[:-1]
                    sets[-1][-1].append((int(line_list[-1]), belief))
                    
    except IOError:
        print "Could not open: " + filename
        sys.exit()
    
    return sets

def read_energies(filename, num_points, num_sets):
    """
    """
    
    # Read all data from file
    all_global_energies, all_indiv_energies, all_total_energies = readcol(filename, 'fff', ',')
    
    # reshape the data into different sets
    global_energy_sets = [ all_global_energies[i*num_points:(i+1)*num_points] for i in xrange(num_sets)]
    indiv_energy_sets = [ all_indiv_energies[i*num_points:(i+1)*num_points] for i in xrange(num_sets)]
    total_energy_sets = [ all_total_energies[i*num_points:(i+1)*num_points] for i in xrange(num_sets)]
    
    return global_energy_sets, indiv_energy_sets, total_energy_sets

def read_cpp_output(filename, num_points, num_nodes, num_sets):
    """
    returns a list of datapoints for each set
    
    """
    
    # Attempt to open file
    try:
        data_file = open(filename, 'r')
    except IOError:
        print "could not open: " + filename
        sys.exit()
    
    # Loop through each line
    sets = []
    for i in xrange(num_sets):
        set = []
        for j in xrange(num_points):
            system_state = []
            for k in xrange(num_nodes):
                line = data_file.readline().strip().split(',')
                system_state.append( (line[:-1], float(line[-1])) )
                
            set.append(system_state)
        
        sets.append(set)
    
    data_file.close()    
    
    return sets

def clear_data_files(model_parameters):
    """
    """
    try:
        file = open(model_parameters['outfile'], 'w')
        file.write("")
        file.close()
    except:
        pass
    try:
        file = open(model_parameters['sizesfile'], 'w')
        file.write("")
        file.close()
    except:
        pass
    try:
        file = open(model_parameters['out_energy_file'], 'w')
        file.write("")
        file.close()
    except:
        pass

def bin_data(time_series, bin_width):

    return time_series[:(time_series.size // bin_width) * bin_width].reshape(-1, bin_width).sum(axis=1)

def difference_filter(time_series):
    return time_series[1:] - time_series[:-1]

def ratio_filter(timeseries):
    return np.abs(timeseries[1:] / timeseries[:-1])

def find_convergence(parameter_obj_file):

    # Load object file
    model_parameters = load_object(parameter_obj_file)

    # Time
    tmod =  model_parameters['niterations'] / model_parameters['npoints']
    trange = [ tmod * i for i in range(model_parameters['npoints']) ]
    
    # Get data
    global_energy_sets, indiv_energy_sets, total_energy_sets = read_energies(model_parameters['out_energy_file'], model_parameters['npoints'], model_parameters['num_sets'])

    plt.clf()
    # For each set, find the convergence
    # mp.dps = 50
    # mp.prec = 53
    global_mean_energy = [ np.mean(time_point) for time_point in zip(*global_energy_sets) ]
    bin_energies = bin_data( np.array(global_mean_energy), 20) / 20.
    energies = np.cumsum(bin_energies)
    energies = ratio_filter(np.cumsum(bin_energies))
    plt.axhline(1.0, lw=2, color='black', ls='--')
    # global_shanks = mp.shanks(energies)[-1][-1]
    # with mp.extraprec(2 * mp.prec): # levin needs a high working precision
        # levin_model = mp.levin(method = "levin", variant = "u")
        # bestLevin, errLevin = levin_model.update_psum(energies)

    plt.plot(energies, marker='o')
    # plt.axhline(global_shanks, lw=2, ls='--',color='green')
    # plt.axhline(bestLevin, lw=2, ls='--',color='red')
    # plt.yscale('log')
    # plt.xscale('log')

    # S = [1. / n / n for n in xrange(1,len(energies))]
    # plt.plot(S)
    plt.savefig('test-conv-I2-ratio_plot.png')

    # for i, global_energy_set in enumerate(global_energy_sets):
    #     bin_energies = mp.one * bin_data(np.array(global_energy_set), 100) / 100.
    #     energies = ratio_filter(bin_energies)
    #     cum_energies = np.cumsum(bin_energies)
    #     # rates = difference_filter(np.array(bin_energies))
    #     # drates = difference_filter(rates)
        
    #     # print energies
    #     with mp.extraprec(2 * mp.prec): # levin needs a high working precision
    #         levin_model = mp.levin(method = "levin", variant = "u")
    #         bestLevin, errLevin = levin_model.update_psum(cum_energies)

    #     global_shanks = mp.shanks(cum_energies)[-1][-1]
        
        # plt.plot(energies, marker='o')
    #     # plt.plot(global_energy_set, marker='o')
    #     # plt.axhline(global_shanks, lw=2)
    #     # plt.yscale('log')
    #     # plt.xscale('log')
    #     plt.axhline(bestLevin, lw=2, color='red')

    # # plt.ylim(0, 9000)
    # plt.show()

def plot_energies(prefix, parameter_obj_file, xdim=300, ydim=300):
    """
    """
    
    # Load object file
    model_parameters = load_object(parameter_obj_file)

    # Time
    tmod =  model_parameters['niterations'] / model_parameters['npoints']
    trange = np.array([ tmod * i for i in range(model_parameters['npoints']) ]) / float(model_parameters['N'])
    
    # Get data
    global_energy_sets, indiv_energy_sets, total_energy_sets = read_energies(model_parameters['out_energy_file'], model_parameters['npoints'], model_parameters['num_sets'])

    # Find mean of data
    global_mean_energy = [ np.mean(time_point) for time_point in zip(*global_energy_sets) ]
    indiv_mean_energy = [ np.mean(time_point) for time_point in zip(*indiv_energy_sets) ]
    total_mean_energy = [ np.mean(time_point) for time_point in zip(*total_energy_sets) ]

    # Plotting
    plt.clf()

    f, (ax, bx, cx) = plt.subplots(3, sharex=False)
    f.subplots_adjust(hspace=0.1)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

    for eg in global_energy_sets:
        ax.plot(trange, eg, ls='-', marker=None, color='0.75')
    for ei in indiv_energy_sets:
        bx.plot(trange, ei, ls='-', marker=None, color='0.75')
    for ee in total_energy_sets:
        cx.plot(trange, ee, ls='-', marker=None, color='0.75')
        
    bx.plot(trange, indiv_mean_energy, color='black', ls='-', linewidth=2.0, marker=None)
    ax.plot(trange, global_mean_energy, color='black', ls='-', linewidth=2.0, marker=None)
    cx.plot(trange, total_mean_energy, color='black', ls='-', linewidth=2.0, marker=None)

    ax.set_ylim(0.0, ax.get_ylim()[1])
    bx.set_ylim(-1,1.)
    bx.set_xlim(0, max(trange))
    ax.set_xlim(0, max(trange))
    cx.set_xlim(0, max(trange))
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.yaxis.get_major_formatter().set_powerlimits((4, 4))
    ax.xaxis.get_major_formatter().set_powerlimits((-40, 40))
    cx.xaxis.get_major_formatter().set_powerlimits((-40, 40))
    bx.xaxis.get_major_formatter().set_powerlimits((-4, 4))
    ax.set_ylabel(r'$E^{(g)}$')
    bx.set_ylabel(r'$\langle E^{(i)} \rangle$')
    cx.set_xlabel(r'$Time/N$')
    cx.set_ylabel(r'$H$')
    
    plt.savefig(prefix + '-energies.pdf' ,dpi=1200,bbox_inches='tight')
    plt.close()

def energies_plot_pubversion(prefix, parameter_obj_file, xdim=300, ydim=300):
    """
    """
    
    # Load object file
    model_parameters = load_object(parameter_obj_file)

    # Time
    tmod =  model_parameters['niterations'] / model_parameters['npoints']
    trange = np.array([ tmod * i for i in range(model_parameters['npoints']) ]) / float(model_parameters['N'])
    
    # Get data
    global_energy_sets, indiv_energy_sets, total_energy_sets = read_energies(model_parameters['out_energy_file'], model_parameters['npoints'], model_parameters['num_sets'])

    # Find mean of data
    global_mean_energy = [ np.mean(time_point) for time_point in zip(*global_energy_sets) ]
    indiv_mean_energy = [ np.mean(time_point) for time_point in zip(*indiv_energy_sets) ]

    # Plotting
    plt.clf()

    f, (ax, bx) = plt.subplots(2, sharex=False)
    f.subplots_adjust(hspace=0.1)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

    # for eg in global_energy_sets:
    #     ax.plot(trange, eg, ls='-', marker=None, color='0.75')
    # for ei in indiv_energy_sets:
    #     bx.plot(trange, ei, ls='-', marker=None, color='0.75')
        
    bx.plot(trange, indiv_mean_energy, color='black', ls='-', linewidth=2.0, marker=None)

        
    ax.plot(trange, global_mean_energy, color='black', ls='-', linewidth=2.0, marker=None)

    
    ax.set_ylim(0.0, ax.get_ylim()[1])
    bx.set_ylim(-1,1.)
    bx.set_xlim(0, max(trange))
    ax.set_xlim(0, max(trange))
    # cx.set_xlim(0, max(trange))
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.yaxis.get_major_formatter().set_powerlimits((4, 4))
    ax.xaxis.get_major_formatter().set_powerlimits((-40, 40))
    bx.xaxis.get_major_formatter().set_powerlimits((-4, 4))
    ax.set_ylabel(r'$E^{(s)}$')
    bx.set_ylabel(r'$\langle E^{(i)} \rangle$')
    bx.set_xlabel(r'$Time$')
    plt.tight_layout()
    plt.savefig(prefix + '-energies.pdf' ,dpi=1200,bbox_inches='tight')
    plt.close()

def sizes_plot_pubversion(prefix, parameter_obj_file):
    """
    """

    # Load object file
    model_parameters = load_object(parameter_obj_file)

    # Time
    tmod =  model_parameters['niterations'] / model_parameters['npoints']
    trange = np.array([ tmod * i for i in range(model_parameters['npoints']) ]) / float(model_parameters['N'])
    
    # Get data
    sets = read_sizes_output(model_parameters['sizesfile'], model_parameters['npoints'])

    # Calculate average over all sets
    largest_groups = []
    second_groups = []
    for j, size_set in enumerate(sets):
        # Get zealot time-series
        largest = [ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in size_set ]
        try:
            second = [ sorted(point, reverse=True)[1][0] / float(model_parameters['N']) for point in size_set ]
        except:
            second = [ 0.0 for point in size_set ]
        
        largest_groups.append(largest)
        second_groups.append(second)
    
    largest_groups = zip(*largest_groups)
    second_groups = zip(*second_groups)

    sizes_first = [ np.mean(time_point) for time_point in largest_groups ]
    sizes_second = [ np.mean(time_point) for time_point in second_groups ]

    # Plotting
    plt.clf()
    plt.plot(trange, sizes_first, color='black', ls='-', linewidth=2.0, marker=None, label='$S_1$')
    plt.plot(trange, sizes_second, color='red', ls='--', linewidth=2.0, marker=None, label='$S_2$')    
    plt.ylim(0.0, 1.0)
    plt.ylabel(r'$\langle S_g/N \rangle$')
    plt.xlabel(r'$Time$')
    plt.legend(loc=7, frameon=False)
    plt.tight_layout()
    plt.savefig(prefix + '-sizes.pdf' ,dpi=1200,bbox_inches='tight')
    plt.close()

def unstable_collapes_1col(prefix, parameter_obj_file):
    """
    """

    # Load object file
    model_parameters = load_object(parameter_obj_file)

    # Time
    tmod =  model_parameters['niterations'] / model_parameters['npoints']
    trange = np.array([ tmod * i for i in range(model_parameters['npoints']) ]) / float(model_parameters['N'])
    
    # Get data
    sets = read_sizes_output(model_parameters['sizesfile'], model_parameters['npoints'])    

    # Calculate average over all sets
    mainstream_groups = []
    winner_groups = []
    for j, size_set in enumerate(sets):

        winning_belief = map(int, sorted(size_set[-1])[-1][1])
        # Get zealot time-series
        mainstream = [ find_belief_size(point, model_parameters['belief_types'][j][0]) / float(model_parameters['N']) for point in size_set ]
        winner = [ find_belief_size(point, winning_belief) / float(model_parameters['N']) for point in size_set ]
        mainstream_groups.append(mainstream)
        winner_groups.append(winner)
    
    mainstream_groups = zip(*mainstream_groups)
    winner_groups = zip(*winner_groups)

    sizes_main = [ np.mean(time_point) for time_point in mainstream_groups ]
    sizes_winner = [ np.mean(time_point) for time_point in winner_groups ]

    # Energies
    global_energy_sets, indiv_energy_sets, total_energy_sets = read_energies(model_parameters['out_energy_file'], model_parameters['npoints'], model_parameters['num_sets'])

    # Find mean of data
    global_mean_energy = [ np.mean(time_point) for time_point in zip(*global_energy_sets) ]
    indiv_mean_energy = [ np.mean(time_point) for time_point in zip(*indiv_energy_sets) ]

    # Plotting
    plt.clf()
    f, (ax, bx, cx) = plt.subplots(3, sharex=True, figsize=(8, 10))
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

    ax.plot(trange, global_mean_energy, color='black', ls='-', linewidth=2.0, marker=None)
    bx.plot(trange, indiv_mean_energy, color='black', ls='-', linewidth=2.0, marker=None)

    ax.set_ylim(0.0, ax.get_ylim()[1])
    bx.set_ylim(-1,1.)
    bx.set_xlim(0, max(trange))
    ax.set_xlim(0, max(trange))
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.yaxis.get_major_formatter().set_powerlimits((4, 4))
    ax.xaxis.get_major_formatter().set_powerlimits((-40, 40))
    bx.xaxis.get_major_formatter().set_powerlimits((-4, 4))
    ax.set_ylabel(r'$E^{(s)}-E^{(s)}_{\mathrm{min}}$')
    bx.set_ylabel(r'$\langle E^{(i)} \rangle$')

    cx.plot(trange, sizes_main, color='black', ls='-', linewidth=2.0, marker=None, label='$S_o$ ($E=1$)')
    cx.plot(trange, sizes_winner, color='red', ls='--', linewidth=2.0, marker=None, label='$S_f$ ($E=-1$)')
    cx.set_ylim(0.0, 1.0)
    cx.set_ylabel(r'$\langle S_g/N \rangle$')
    cx.set_xlabel(r'$t$')
    cx.legend(loc=7, frameon=False)

    ax.get_yaxis().set_label_coords(-0.1,0.5)
    bx.get_yaxis().set_label_coords(-0.1,0.5)
    cx.get_yaxis().set_label_coords(-0.1,0.5)

    plt.tight_layout()
    f.subplots_adjust(hspace=0.1)
    plt.savefig(prefix + '-1col.pdf' ,dpi=1200),#bbox_inches='tight')
    plt.close()

def plot_sizes(prefix, model_parameters, size_file=None):
    """
    """

    # Time
    tmod =  model_parameters['niterations'] / model_parameters['npoints']
    trange = np.array([ tmod * i for i in range(model_parameters['npoints']) ]) / float(model_parameters['N'])
    
    # Get data
    if size_file == None:
        sets = read_sizes_output(model_parameters['sizesfile'], model_parameters['npoints'])
    else:
        sets = read_sizes_output(size_file, model_parameters['npoints'])
    
    # Make all sizes together
    for i, set in enumerate(sets):
        
        largest = [ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in set ]
        plt.plot(trange, largest, color=cm.winter(float(i)/ len(sets))) 
        try: # problem going out of index range when no second largest
            second = [ sorted(point, reverse=True)[1][0] / float(model_parameters['N']) for point in set ] 
            plt.plot(trange, second, color=cm.autumn(float(i)/ len(sets)))
        except:
            pass
    
    plt.ylim(0.0,1.0)
    plt.savefig(prefix + '-sizes.png', dpi=150,bbox_inches='tight')
    plt.clf()    
    plt.close()

def concatenate_files(outfile_name, filenames):
    """
    """
    
    # Concatenate files
    with open(outfile_name, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    
def delete_files(filenames):
    """
    """
    
    for file in filenames:
        os.remove(file)

def run_parallel(nprocs, model_parameters):
    """
    """
    
    # Assumes graph file is already written
    
    def worker(out_q, id, parameters):
        """
        """
        
        # Write process parameter file
        write_cparam_file(parameters['paramfile'], parameters)
        # Reset any data files
        clear_data_files(parameters)
        # Run program
        run_cpp(parameters['cpp_exe'], parameters['paramfile'])
        
        out_q.put({id:True})
        
    chunksize = int(math.ceil(model_parameters['num_sets'] / float(nprocs)))
    
    # get chunks of info for each process
    if model_parameters['random'] == False:
        belief_chunks = [ model_parameters['belief_types'][chunksize * i:chunksize * (i + 1)] for i in xrange(nprocs) ]
    else:
        belief_chunks = [ [] for i in xrange(nprocs) ]
    if model_parameters['zealots'] == True:
        zealot_chunks = [ model_parameters['zealot_types'][chunksize * i:chunksize * (i + 1)] for i in xrange(nprocs) ]
    else:
        zealot_chunks = [ [] for i in xrange(nprocs) ]
        
    # generate model parameters for each process
    model_param_list = []
    for i in xrange(nprocs):
        chunk_params = copy.deepcopy(model_parameters)
        chunk_params['belief_types'] = belief_chunks[i]
        chunk_params['zealot_types'] = zealot_chunks[i]
        chunk_params['seed'] = rnd.randint(1,50000)
        chunk_params['paramfile'] = str(i) + '_' + model_parameters['paramfile']
        chunk_params['outfile'] = str(i) + '_' + model_parameters['outfile']
        chunk_params['sizesfile'] = str(i) + '_' + model_parameters['sizesfile']
        chunk_params['out_energy_file'] = str(i) + '_' + model_parameters['out_energy_file']
        if model_parameters['random'] == True and model_parameters['zealots'] == False:
            chunk_params['num_sets'] = chunksize
        elif model_parameters['random'] == False:
            chunk_params['num_sets'] = len(belief_chunks[i])
            if len(belief_chunks[i]) == 0:
                continue
        elif model_parameters['zealots'] == True:
            chunk_params['num_sets'] = len(zealot_chunks[i])
            if len(zealot_chunks[i]) == 0:
                continue
        
        model_param_list.append(chunk_params)
    
    # Start process
    procs = []
    out_q = Queue()
    for i, params in enumerate(model_param_list):
        p = Process(target=worker, args=(out_q, i, params))
        procs.append(p)
        p.start()

    # Collect successful id's
    resultdict = {}
    for i in range(len(procs)):
        resultdict.update(out_q.get())

    for i in xrange(len(procs)):
        procs[i].join()

    # Post processing files

    # Process energy files
    energy_files = [ params['out_energy_file'] for params in model_param_list]
    concatenate_files(model_parameters['out_energy_file'], energy_files)
    delete_files(energy_files)
    
    # process size files
    size_files = [ params['sizesfile'] for params in model_param_list ]
    concatenate_files(model_parameters['sizesfile'], size_files)
    delete_files(size_files)

    # process out files
    out_files = [ params['outfile'] for params in model_param_list ]
    concatenate_files(model_parameters['outfile'], out_files)
    delete_files(out_files)
    
    # clear redundant parameter files
    delete_files([model['paramfile'] for model in model_param_list])

def run_zealot_density_analysis(basic_model_parameters, densities, nprocesses):
    """
    Assumes graph file already created
    
    The model parameters for the simulation. THese will be modified
    accordingly.
    
    densities = list of densities to run the simulation on
    
    nprocesses = number of processes for parallelization
    
    """

    # Run simulation for each density   
    for density in densities:
        model_parameters = copy.deepcopy(basic_model_parameters)
        # Reset ratios
        model_parameters['belief_ratios'] = [1.0 - density]
        model_parameters['zealot_ratios'] = [density]
        # Reset filenames
        model_parameters['paramfile'] = str(density)+ '_' + basic_model_parameters['paramfile']
        model_parameters['outfile'] = str(density)+ '_' + basic_model_parameters['outfile']
        model_parameters['sizesfile'] = str(density)+'_' + basic_model_parameters['sizesfile']
        model_parameters['out_energy_file'] = str(density)+'_' + basic_model_parameters['out_energy_file']
        
        # Run
        run_parallel(nprocesses, model_parameters)

def run_multigraph_zealot_density_analysis(prefix, basic_model_parameters, densities, nprocesses, ngraphs, point_range=(195,200)):
    """
    """

    # Save model parameter object
    save_object(basic_model_parameters, basic_model_parameters['paramfile'] + '.pyobj')

    # Run for each graph
    list_of_tsmeans_for_each_density_for_each_graph = []
    for i in xrange(ngraphs):

        # Create graph
        graph = get_connected_subgraph(nx.fast_gnp_random_graph(basic_model_parameters['N'], basic_model_parameters['multigraph_p']))
        write_graph_file(basic_model_parameters['graph'], graph)
        del graph

        # Run simulation for each density   
        for density in densities:
            model_parameters = copy.deepcopy(basic_model_parameters)
            # Reset ratios
            model_parameters['belief_ratios'] = [1.0 - density]
            model_parameters['zealot_ratios'] = [density]
            # Reset filenames
            model_parameters['paramfile'] = str(density)+ '_' + basic_model_parameters['paramfile']
            model_parameters['outfile'] = str(density)+ '_' + basic_model_parameters['outfile']
            model_parameters['sizesfile'] = str(density)+'_' + basic_model_parameters['sizesfile']
            model_parameters['out_energy_file'] = str(density)+'_' + basic_model_parameters['out_energy_file']
            
            # Run
            run_parallel(nprocesses, model_parameters)

        # run density plot and save needed graph data
        list_of_density_files = [ str(density) + '_' + basic_model_parameters['sizesfile'] for density in basic_model_parameters['densities'] ]
        
        # read in sets
        list_of_tsmeans_for_each_density = []
        for density_size_file in list_of_density_files:
            sets = read_sizes_output(density_size_file, basic_model_parameters['npoints'])
            
            # Calculate average over all sets
            zealot_groups = []
            for j, set in enumerate(sets):
                # Get zealot time-series
                zealot = [ find_zealot_size(point, basic_model_parameters['zealot_types'][j][0]) / float(basic_model_parameters['N']) for point in set ]
                
                zealot_groups.append(zealot)
                        
            # Take mean over all sets
            zealot_groups_tsmean = [ np.mean(set_ts[point_range[0]: point_range[1]]) for set_ts in zealot_groups ]
            list_of_tsmeans_for_each_density.append(zealot_groups_tsmean)

        list_of_tsmeans_for_each_density_for_each_graph.append(list_of_tsmeans_for_each_density)

    list_of_tsmeans_for_each_graph_for_each_density = zip(*list_of_tsmeans_for_each_density_for_each_graph)
    list_of_tsmeans_for_each_density_graphs_folded = [ [ item for subsublist in sublist for item in subsublist] for \
        sublist in list_of_tsmeans_for_each_graph_for_each_density  ]
    save_object((basic_model_parameters['densities'], list_of_tsmeans_for_each_density_graphs_folded), prefix + '_multigraph.pltdat')

def plot_density_pltdat(prefix, pltdat_files, labels):
    # Marker sequence
    markers = ['o', '<','s','D']
    colors = ['darkblue', 'royalblue','skyblue']
    # create zealot density plotter using saved graph data
    for i, pltdat_file in enumerate(pltdat_files):
        # load data for file
        densities, means_per_point = load_object(pltdat_file)
        # calculate errors
        # std_of_largest = np.std(means_per_point, 1)
        mean_err = stats.sem(means_per_point, 1)
        top_err = mean_err
        bot_err = mean_err
        # calculate mean
        mean_of_largest = np.mean(means_per_point, 1)
        # plot w/ some set of shapes and labels
        plt.errorbar(densities, mean_of_largest, yerr=[bot_err, top_err],\
                 marker=markers[i], linewidth=1, color=colors[i], \
                 markerfacecolor='none', markersize=8, label=labels[i],
                 markeredgecolor=colors[i], markeredgewidth=1)

    # Write plots out to file
    plt.xlabel(r'$\rho_o$')
    plt.ylabel(r'$\langle S_z/N \rangle$')
    plt.legend(frameon=False, loc='lower right')
    plt.ylim(0.0, 1.02)
    plt.xlim(0.0, 0.15)
    plt.tight_layout()
    plt.savefig(prefix + '-zealdensity-combo.pdf',dpi=1200,bbox_inches='tight')
    plt.clf()
    plt.close()

def run_mu_analysis(main_param_dict, nprocesses):
    """
    """

    # Get graph file locations from folder
    folder_list = os.listdir(main_param_dict['graphs_folder'])
    model_parameter_file_list = []
    for i, mu in enumerate(main_param_dict['mu_list']):
        os.chdir(main_param_dict['graphs_folder'] + '/' + str(mu) + '-mu_lfr_graph')
        
        model_parameters = copy.deepcopy(main_param_dict)
        model_parameters['graph'] = str(mu) + '-mu_' + model_parameters['graph']
        model_parameters['paramfile'] = str(mu) + '-mu_' + model_parameters['paramfile']
        model_parameters['outfile'] = str(mu) + '-mu_' + model_parameters['outfile']
        model_parameters['sizesfile'] = str(mu) + '-mu_' + model_parameters['sizesfile']
        model_parameters['out_energy_file'] = str(mu) + '-mu_' + model_parameters['out_energy_file']
        model_parameters['mu'] = mu
        
        # Write graph to file
        graph = read_LFR_Benchmark('network.dat', 'community.dat')
        # REset working directory
        os.chdir('../..')
        
        if model_parameters['assignment'] == True:
            graph = relabel_graph(graph, model_parameters['assignment_dict'])
            model_parameters['assignments'] = community_assignment(graph, model_parameters['assignment_dict'])
        write_graph_file(model_parameters['graph'], graph)
        
        del graph
        
        # Write process parameter file
        write_cparam_file(model_parameters['paramfile'], model_parameters)
        save_object(model_parameters, model_parameters['paramfile'] + '.pyobj')
        model_parameter_file_list.append(model_parameters['paramfile'] + '.pyobj')
        # Reset any data files
        clear_data_files(model_parameters)
        # Run program
        run_parallel(nprocesses, model_parameters)

    save_object(model_parameter_file_list, main_param_dict['paramfile'] + '_filelist.pyobj')

def two_variable_analysis(main_param_dict, nprocesses, target_vars, tar_var1, tar_var2):
    """
    VAR must not be MU

    target_vars is a tuple that contains the two variables you want to vary
    two variables should be lists (made by arange or linspace or something), the other two should just be scalars
    the two non scalars should be 
    """

    model_parameter_file_list = []
    # Loop through target variable 1
    for i, var1 in enumerate(tar_var1):
        # Loop through target variable 2
        for j, var2 in enumerate(tar_var2):

            # Copy main parameters
            model_parameters = copy.deepcopy(main_param_dict)
            model_parameters['paramfile'] = str(var1) + '-var1_' + str(var2) + '-var2_' + model_parameters['paramfile']
            model_parameters['outfile'] = str(var1) + '-var1_' + str(var2) + '-var2_' + model_parameters['outfile']
            model_parameters['sizesfile'] = str(var1) + '-var1_' + str(var2) + '-var2_'+ model_parameters['sizesfile']
            model_parameters['out_energy_file'] = str(var1) + '-var1_' + str(var2) + '-var2_' + model_parameters['out_energy_file']

            # Set variable parameters
            if target_vars[0] == 'I':
                model_parameters['I'] = var1 
            elif target_vars[0] == 'J':
                model_parameters['J'] = var1 
            elif target_vars[0] == 'T':
                model_parameters['T'] = var1 

            if target_vars[1] == 'I':
                model_parameters['I'] = var2
            elif target_vars[1] == 'J':
                model_parameters['J'] = var2 
            elif target_vars[1] == 'T':
                model_parameters['T'] = var2
            
            # Write process parameter file
            write_cparam_file(model_parameters['paramfile'], model_parameters)
            save_object(model_parameters, model_parameters['paramfile'] + '.pyobj')
            model_parameter_file_list.append(model_parameters['paramfile'] + '.pyobj')
            # Reset any data files
            clear_data_files(model_parameters)
            # Run program
            run_parallel(nprocesses, model_parameters)

    save_object(model_parameter_file_list, main_param_dict['paramfile'] + '_filelist.pyobj')

def two_variable_mu_analysis(main_param_dict, nprocesses, target_vars, tar_var1, tar_var2):
    """
    VAR1 MUST BE MU

    target_vars is a tuple that contains the two variables you want to vary
    two variables should be lists (made by arange or linspace or something), the other two should just be scalars
    the two non scalars should be 
    """

    # Get graph file locations from folder
    folder_list = os.listdir(main_param_dict['graphs_folder'])
    model_parameter_file_list = []
    # Loop through target variable 1
    for i, var1 in enumerate(tar_var1):
        # Loop through target variable 2
        for j, var2 in enumerate(tar_var2):
            # Change directory to community graph folder
            os.chdir(main_param_dict['graphs_folder'] + '/' + str(var1) + '-mu_lfr_graph')

            # Copy main parameters
            model_parameters = copy.deepcopy(main_param_dict)
            model_parameters['graph'] = str(var1) + '-var1_' + str(var2) + '-var2_' + model_parameters['graph']
            model_parameters['paramfile'] = str(var1) + '-var1_' + str(var2) + '-var2_' + model_parameters['paramfile']
            model_parameters['outfile'] = str(var1) + '-var1_' + str(var2) + '-var2_' + model_parameters['outfile']
            model_parameters['sizesfile'] = str(var1) + '-var1_' + str(var2) + '-var2_'+ model_parameters['sizesfile']
            model_parameters['out_energy_file'] = str(var1) + '-var1_' + str(var2) + '-var2_' + model_parameters['out_energy_file']

            # Set variable parameters
            if target_vars[0] == 'mu':
                model_parameters['mu'] = var1 
            elif target_vars[0] == 'I':
                model_parameters['I'] = var1 
            elif target_vars[0] == 'J':
                model_parameters['J'] = var1 
            elif target_vars[0] == 'T':
                model_parameters['T'] = var1 

            if target_vars[1] == 'mu':
                model_parameters['mu'] = var2
            elif target_vars[1] == 'I':
                model_parameters['I'] = var2
            elif target_vars[1] == 'J':
                model_parameters['J'] = var2 
            elif target_vars[1] == 'T':
                model_parameters['T'] = var2

            # Write graph to file
            graph = read_LFR_Benchmark('network.dat', 'community.dat')
            # REset working directory
            os.chdir('../..')
            
            if model_parameters['assignment'] == True:
                graph = relabel_graph(graph, model_parameters['assignment_dict'])
                model_parameters['assignments'] = community_assignment(graph, model_parameters['assignment_dict'])
            write_graph_file(model_parameters['graph'], graph)
            
            del graph
            
            # Write process parameter file
            write_cparam_file(model_parameters['paramfile'], model_parameters)
            save_object(model_parameters, model_parameters['paramfile'] + '.pyobj')
            model_parameter_file_list.append(model_parameters['paramfile'] + '.pyobj')
            # Reset any data files
            clear_data_files(model_parameters)
            # Run program
            run_parallel(nprocesses, model_parameters)

    save_object(model_parameter_file_list, main_param_dict['paramfile'] + '_filelist.pyobj')

def plot_two_var_contour(prefix, param_file_list_file_name, point_range, var1_key, var2_key, logvar1=True, logvar2=True, switch=False):
    """
    Keywords must be in order according to how the plot was created, use switch to switch the axes
    """

    # Load file list
    param_file_list = load_object(param_file_list_file_name)

    # Loop through each parameter file and get parameter values
    mean_largest_sizes = []
    mean_second_sizes = []
    var1_list = []
    var2_list = []
    for param_file in param_file_list:
        # Load parameters and load size data
        model_parameters = load_object(param_file)
        var1_list.append(model_parameters[var1_key])
        var2_list.append(model_parameters[var2_key])

        # Calculate average of all sets
        sets = read_sizes_output(model_parameters['sizesfile'], model_parameters['npoints'])
        largest_groups = []
        second_groups = []
        for j, output_set in enumerate(sets):
            # Get zealot time-series
            largest = [ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in output_set ]
            try:
                second = [ sorted(point, reverse=True)[1][0] / float(model_parameters['N']) for point in output_set ]
            except:
                second = [ 0.0 for point in output_set ]
            
            largest_groups.append(largest)
            second_groups.append(second)
        
        largest_groups = zip(*largest_groups)
        second_groups = zip(*second_groups)
        
        # Take mean over all sets
        largest_means = [ np.mean(largest_groups[i]) for i in xrange(point_range[0], point_range[1])]
        second_means = [ np.mean(second_groups[i]) for i in xrange(point_range[0], point_range[1]) ]
        
        # Take mean over range (assumes stationary average)
        largest_mean = np.mean(largest_means)
        second_mean = np.mean(second_means)
        
        mean_largest_sizes.append(largest_mean)
        mean_second_sizes.append(second_mean)

    # Combine into gridspace
    var1_range = sorted(list(set(var1_list)))
    var2_range = sorted(list(set(var2_list)))
    vv1, vv2 = np.meshgrid(var1_range, var2_range, indexing='ij')
    zz_large = np.zeros(vv1.shape)
    zz_second = np.zeros(vv1.shape)
    for iii in xrange(len(var1_range)):
        for jjj in xrange(len(var2_range)):
            if (var1_list[iii * len(var2_range) + jjj] == var1_range[iii]) and (var2_list[iii * len(var2_range) + jjj] == var2_range[jjj]):
                zz_large[iii, jjj] = mean_largest_sizes[ iii * len(var2_range) + jjj ]
                zz_second[iii, jjj] = mean_second_sizes[ iii * len(var2_range) + jjj ]
    
    # Plot on two contours, one for largest and one for second largest
    levels = np.linspace(0.0, 1.00001, 11)
    if switch:
        CS = plt.contourf(vv2, vv1, zz_large, levels, cmap=cm.PuBu)        
        plt.ylabel('$' + var1_key +'$')
        plt.xlabel('$' + var2_key +'$')
    else:
        CS = plt.contourf(vv1, vv2, zz_large, levels, cmap=cm.PuBu)
        plt.xlabel('$' + var1_key +'$')
        plt.ylabel('$' + var2_key +'$')

    bar = plt.colorbar(CS)
    bar.ax.set_ylabel(r'$\left \langle S_g/N \right \rangle$')
    if logvar1 == True:
        plt.xscale("log")
    if logvar2 == True:
        plt.yscale("log")
    plt.tight_layout()
    plt.savefig(prefix + "_2var_contour.pdf")
    plt.clf()
    plt.close()

def plot_size_vs_mu(prefix, param_file_list_file_name, point_range, logged=False):
    """
    """
    
    # Load file list
    param_file_list = load_object(param_file_list_file_name)
    
    # Loop through each parameter file
    mean_largest_sizes = []
    mean_second_sizes = []
    mus = []
    for param_file in param_file_list:
        # Load parameters and load size data
        model_parameters = load_object(param_file)
        sets = read_sizes_output(model_parameters['sizesfile'], model_parameters['npoints'])
        mus.append(model_parameters['mu'])
        
        # Calculate average over all sets
        largest_groups = []
        second_groups = []
        for j, set in enumerate(sets):
            # Get zealot time-series
            largest = [ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in set ]
            try:
                second = [ sorted(point, reverse=True)[1][0] / float(model_parameters['N']) for point in set ]
            except:
                second = [ 0.0 for point in set ]
            
            largest_groups.append(largest)
            second_groups.append(second)
        
        largest_groups = zip(*largest_groups)
        second_groups = zip(*second_groups)
        
        # Take mean over all sets
        largest_means = [ np.mean(largest_groups[i]) for i in xrange(point_range[0], point_range[1])]
        second_means = [ np.mean(second_groups[i]) for i in xrange(point_range[0], point_range[1]) ]
        
        # Take mean over range (assumes stationary average)
        largest_mean = np.mean(largest_means)
        second_mean = np.mean(second_means)
        
        mean_largest_sizes.append(largest_mean)
        mean_second_sizes.append(second_mean)
    
    plt.plot(mus, mean_largest_sizes, color='black', marker='o', \
             markerfacecolor='none', markersize=8, markeredgewidth=1, \
             linewidth=1, markeredgecolor='black')
    plt.plot(mus, mean_second_sizes, color='red', marker='^', \
             markerfacecolor='none', markersize=8, markeredgewidth=1, \
             linewidth=1, markeredgecolor='red')


    if logged == True:
        plt.xscale("log")
    plt.xlim(min(mus), max(mus))
    plt.ylim(0.0, 1.02)
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'$\langle S_g/N \rangle$')
    plt.savefig(prefix + '-community_rigidity.pdf',dpi=1200,bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_second_size_vs_mu(prefix, param_file_list_file_names, label_list, point_range, logged=False):
    """
    """

    # Marker sequence
    markers = ['o', '<','s','D']
    colors = ['black', 'blue','purple','teal']
    
    for m, param_file_list_file_name in enumerate(param_file_list_file_names):
        # Load file list
        param_file_list = load_object(param_file_list_file_name)
        
        # Loop through each parameter file
        mean_largest_sizes = []
        mean_second_sizes = []
        largest_errs = []
        mus = []
        for param_file in param_file_list:
            # Load parameters and load size data
            model_parameters = load_object(param_file)
            sets = read_sizes_output(model_parameters['sizesfile'], model_parameters['npoints'])
            mus.append(model_parameters['mu'])
            
            # Calculate average over all sets
            largest_groups = []
            second_groups = []
            for j, set in enumerate(sets):
                # Get zealot time-series
                largest = Normalize([ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in set ], 0.5, 1.0)
                try:
                    second = Normalize([ sorted(point, reverse=True)[1][0] / float(model_parameters['N']) for point in set ], 0.5, 1.0)
                except:
                    second = [ 0.0 for point in set ]
                
                largest_groups.append(largest)
                second_groups.append(second)
            
            largest_groups = zip(*largest_groups)
            second_groups = zip(*second_groups)
            
            # Take mean over all sets
            largest_means = [ np.mean(largest_groups[i]) for i in xrange(point_range[0], point_range[1]) ]
            second_means = [ np.mean(second_groups[i]) for i in xrange(point_range[0], point_range[1]) ]
            
            # Take mean over range (assumes stationary average)
            largest_mean = np.mean(largest_means)
            second_mean = np.mean(second_means)
            
            largest_err = np.std(largest_means)
            largest_errs.append(largest_err)
            mean_largest_sizes.append(largest_mean)
            mean_second_sizes.append(second_mean)

        plt.plot(mus, mean_largest_sizes, color=colors[m], marker=markers[m], \
                 markerfacecolor='none', markersize=8, markeredgewidth=1, \
                 linewidth=1, markeredgecolor=colors[m], label=label_list[m])


    if logged == True:
        plt.xscale("log")
    plt.xlim(min(mus), max(mus))
    plt.ylim(0.0, 1.02)
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'$\textup{Fraction converted}$')
    plt.legend(loc=4, frameon=False)
    plt.tight_layout()
    plt.savefig(prefix + '-community_rigidity.pdf',dpi=1200,bbox_inches='tight')
    plt.clf()
    plt.close()

def find_belief_size(data_point, belief):
    # Convert belief to string
    string_belief = map(str, belief)
    
    for size_belief_pair in data_point:
        if string_belief == size_belief_pair[1]:
            return size_belief_pair[0]
    
    return 0    

def find_zealot_size(data_point, zealot_belief):
    """
    """
    
    # Convert belief to string
    zealot_belief = map(str, zealot_belief)
    
    for size_belief_pair in data_point:
        if zealot_belief == size_belief_pair[1]:
            return size_belief_pair[0]
    
    print "error: could not match: ", zealot_belief
    print "No match in: ", data_point
    sys.exit()

def zealot_densities_plot(prefix, model_parameters, point_range):
    """
    point_range = range in units of data points over which to do average
    """

    # density files lists
    list_of_density_files = [ str(density) +'_'+ model_parameters['sizesfile'] for density in model_parameters['densities'] ]
    
    # read in sets
    density_means_zealots = []
    for density_size_file in list_of_density_files:
        sets = read_sizes_output(density_size_file, model_parameters['npoints'])
        
        # Calculate average over all sets
        zealot_groups = []
        for j, set in enumerate(sets):
            # Get zealot time-series
            zealot = [ find_zealot_size(point, model_parameters['zealot_types'][j][0]) / float(model_parameters['N']) for point in set ]
            
            zealot_groups.append(zealot)
        
        zealot_groups = zip(*zealot_groups)
        
        # Take mean over all sets
        zealot_means = [ np.mean(zealot_groups[i]) for i in xrange(point_range[0], point_range[1])]
        
        # Take mean over range (assumes stationary average)
        zealot_mean = np.mean(zealot_means)
        
        density_means_zealots.append(zealot_mean)
        
    # Plot the means
    plt.plot(model_parameters['densities'], density_means_zealots)
    plt.savefig(prefix + '-zealotDensity.pdf',dpi=1200,bbox_inches='tight')
    plt.clf()
    plt.close()

def zealot_densities_plot_pubversion(prefix, point_range, parameter_object_file_list, labels):
    """
    """
    
    # Marker sequence
    markers = ['o', '<','s','D']
    colors = ['black', 'blue','purple','teal']
    
    # Loop through and do analysis for each file
    for m, parameter_object_file in enumerate(parameter_object_file_list):
        
        # Load object
        model_parameters = load_object(parameter_object_file)
        
        # density files lists (add + '_')
        list_of_density_files = [ str(density) + '_' + model_parameters['sizesfile'] for density in model_parameters['densities'] ]
        
        # read in sets
        density_means_zealots = []
        for density_size_file in list_of_density_files:
            sets = read_sizes_output(density_size_file, model_parameters['npoints'])
            
            # Calculate average over all sets
            zealot_groups = []
            for j, set in enumerate(sets):
                # Get zealot time-series
                zealot = [ find_zealot_size(point, model_parameters['zealot_types'][j][0]) / float(model_parameters['N']) for point in set ]
                
                zealot_groups.append(zealot)
            
            zealot_groups = zip(*zealot_groups)
            
            # Take mean over all sets
            zealot_means = [ np.mean(zealot_groups[i]) for i in xrange(point_range[0], point_range[1])]
            
            # Take mean over range (assumes stationary average)
            zealot_mean = np.mean(zealot_means)
            
            density_means_zealots.append(zealot_mean)
        
        plt.plot(model_parameters['densities'], density_means_zealots, \
                 marker=markers[m], linewidth=1, color=colors[m], \
                 markerfacecolor='none', markersize=8, label=labels[m],
                 markeredgecolor=colors[m], markeredgewidth=1)

    # Write plots out to file
    plt.xlabel(r'$\rho_o$')
    plt.ylabel(r'$\langle S_g/N \rangle$')
    plt.legend(frameon=False, loc='lower right')
    plt.ylim(0.0, 1.02)
    plt.xlim(0.0, 0.25)
    plt.savefig(prefix + '-zealotDensity.pdf',dpi=1200,bbox_inches='tight')
    plt.clf()
    plt.close()

def run_param_analysis(basic_model_parameters, param_values, param_key, nprocesses):
    """
    Assumes graph file already created
    
    The model parameters for the simulation. THese will be modified
    accordingly.
    
    densities = list of densities to run the simulation on
    
    nprocesses = number of processes for parallelization
    
    """

    # Run simulation for each density   
    for val in param_values:
        model_parameters = copy.deepcopy(basic_model_parameters)
        # Reset ratios
        model_parameters[param_key] = val
        # Reset filenames
        model_parameters['paramfile'] = str(val) + '_' + basic_model_parameters['paramfile']
        model_parameters['outfile'] = str(val) + '_' + basic_model_parameters['outfile']
        model_parameters['sizesfile'] = str(val) + '_' + basic_model_parameters['sizesfile']
        model_parameters['out_energy_file'] = str(val) + '_' + basic_model_parameters['out_energy_file']
        
        # Run
        run_parallel(nprocesses, model_parameters)
        
def plot_param_analysis(prefix, model_parameters, param_key, point_range, xlog=True):
    """
    point_range = range in units of data points over which to do average
    """

    # density files lists
    list_of_density_files = [ str(I) +'_'+ model_parameters['sizesfile'] for I in model_parameters['param_values'] ]
    
    # read in sets
    I_means_largest = []
    for I_size_file in list_of_density_files:
        sets = read_sizes_output(I_size_file, model_parameters['npoints'])
        
        # Calculate average over all sets
        largest_groups = []
        for j, set in enumerate(sets):
            # Get zealot time-series
            largest = [ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in set ]
            
            largest_groups.append(largest)
        
        largest_groups = zip(*largest_groups)
        
        # Take mean over all sets
        largest_means = [ np.mean(largest_groups[i]) for i in xrange(point_range[0], point_range[1])]
        
        # Take mean over range (assumes stationary average)
        largest_mean = np.mean(largest_means)
        
        I_means_largest.append(largest_mean)
        
    # Plot the means
    if xlog==True:
        plt.xscale('log')
    plt.plot(model_parameters['param_values'], I_means_largest)
    plt.savefig(prefix + '-param_plot.pdf',dpi=1200,bbox_inches='tight')
    plt.clf()
    plt.close()

def param_plot_pubversion(prefix, param_key, point_range, parameter_object_file_list, xlog=True, plot_second=True):
    """
    point_range = range in units of data points over which to do average
    """
    
    # Marker sequence
    markers = ['o', '^']
    
    # Loop through and do analysis for each file
    for m, parameter_object_file in enumerate(parameter_object_file_list):
        
        # Load object
        model_parameters = load_object(parameter_object_file)
        list_of_peerinf_files = [ str(I) + "_" + model_parameters['sizesfile'] for I in model_parameters['param_values'] ]
        
        # read in sets
        I_means_largest = []
        Largest_bot_err = []
        Largest_top_err = []
        I_means_second = []
        Second_bot_err = []
        Second_top_err = []
        for I_size_file in list_of_peerinf_files:
            sets = read_sizes_output(I_size_file, model_parameters['npoints'])
            
            # Calculate average over all sets
            largest_groups = []
            second_groups = []
            for j, set in enumerate(sets):
                # Get zealot time-series
                largest = [ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in set ]
                try:
                    second = [ sorted(point, reverse=True)[1][0] / float(model_parameters['N']) for point in set ]
                except:
                    second = [ 0.0 for point in set ]
                
                largest_groups.append(largest)
                second_groups.append(second)

            largest_groups_tsmean = [ np.mean(set_ts[point_range[0]: point_range[1]]) for set_ts in largest_groups ]
            second_groups_tsmean = [ np.mean(set_ts[point_range[0]: point_range[1]]) for set_ts in second_groups ]
            
            # Take mean over range (assumes stationary average)
            largest_mean = np.median(largest_groups_tsmean)
            second_mean = np.median(second_groups_tsmean)

            # Find std
            largest_1Q = np.percentile(largest_groups_tsmean, 25)
            largest_3Q = np.percentile(largest_groups_tsmean, 75)
            second_1Q = np.percentile(second_groups_tsmean, 25)
            second_3Q = np.percentile(second_groups_tsmean, 75)

            Largest_bot_err.append( largest_mean - largest_1Q )
            Largest_top_err.append( largest_3Q - largest_mean )
            Second_bot_err.append( second_mean - second_1Q )
            Second_top_err.append( second_3Q - second_mean ) 
            
            I_means_largest.append(largest_mean)
            I_means_second.append(second_mean)
        
        # Plot the means
        if xlog==True:
            plt.xscale('log')
        plt.errorbar(model_parameters['param_values'], I_means_largest, yerr=[Largest_bot_err,Largest_top_err], \
                 marker='o', linewidth=1, color='black',\
                 markerfacecolor='none', markersize=12, label='Largest',\
                 markeredgecolor='black', markeredgewidth=1)
        if plot_second == True:
            plt.errorbar(model_parameters['param_values'], I_means_second, yerr=[Second_bot_err,Second_top_err],\
                     marker='^', linewidth=1, color='red',\
                     markerfacecolor='none', markersize=12, label='Second largest',\
                     markeredgecolor='red', markeredgewidth=1)
    
    plt.xlabel('$'+param_key+'$')
    plt.ylim(0.0, 1.02)
    plt.ylabel(r'$\mathbf{M}[S_g/N]$') #(r'$\langle S_g/N \rangle$')
    plt.savefig(prefix + '-param_values.pdf',dpi=1200,bbox_inches='tight')
    plt.clf()
    plt.close()

def heatmap_param_pubversion(prefix, param_key, parameter_object_file, point_range=(195,200)):
    """
    """

    # Load object
    model_parameters = load_object(parameter_object_file)
    list_of_peerinf_files = [ str(I) + "_" + model_parameters['sizesfile'] for I in model_parameters['param_values'] ]
    
    # read in sets
    scatterx = []
    scattery = []
    for i, I_size_file in enumerate(list_of_peerinf_files):
        sets = read_sizes_output(I_size_file, model_parameters['npoints'])
        
        # Calculate average over all sets
        largest_groups = []
        second_groups = []
        for j, set in enumerate(sets):
            # Get zealot time-series
            largest = [ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in set ]                
            largest_groups.append(largest)

        largest_groups_tsmean = [ np.mean(set_ts[point_range[0]: point_range[1]]) for set_ts in largest_groups ]
        
        # Take mean over range (assumes stationary average)
        scattery += largest_groups_tsmean
        scatterx += [ model_parameters['param_values'][i] for k in xrange(len(largest_groups_tsmean)) ]

    # Determine histogram2D
    xedges = log_linspace(model_parameters['param_values'][0], model_parameters['param_values'][-1], 20) #30
    yedges = np.linspace(0, 1.0 + (1.0 / 20 / 2.), 20) #30
    # print xedges.shape, yedges.shape, len(scatterx), len(scattery)
    H, xbins, ybins = np.histogram2d(scatterx, scattery, bins=(xedges, yedges), normed=False)

    # Interpolate data
    xcenters = np.log(np.exp(xbins[:-1]) + 0.5 * (np.exp(xbins[1:]) - np.exp(xbins[:-1])))
    ycenters = ybins[:-1] + 0.5 * (ybins[1:] - ybins[:-1])
    Xcenters, Ycenters = np.meshgrid(xcenters, ycenters, indexing='ij')
    # Xcenters, Ycenters = np.meshgrid(xbins, ybins, indexing='ij')
    interX, interY = np.meshgrid(log_linspace(model_parameters['param_values'][0], model_parameters['param_values'][-1], 100), 
        np.linspace(0, 1, 100), indexing='ij')
    grid_H = griddata((np.exp(Xcenters.ravel()), Ycenters.ravel()), H.ravel(), (np.exp(interX), interY), \
        method='linear', fill_value=0.0, rescale=False)

    # Plot
    CA = plt.pcolormesh(interX, interY, grid_H, cmap=cm.gist_heat_r, linewidth=0) #rasterized=True vmin=0.0, vmax = 1.0, 
    plt.xlim(np.min(xcenters), np.max(xcenters))
    plt.ylim(0.0, 1.0)
    plt.xscale('log')
    plt.ylabel('$S_g/N$')
    plt.xlabel('$' + param_key + "$")
    bar = plt.colorbar(CA)
    for label in bar.ax.yaxis.get_ticklabels()[::]:
        label.set_visible(False)
    bar.ax.set_ylabel(r'$\textup{density}$')
    plt.savefig(prefix + '_param_heatmap.png', dpi=600,bbox_inches='tight') #1200
    plt.clf()
    plt.close()

def global_energy_v_I(prefix, parameter_object_file_list, xlog=True):
    """
    point_range = range in units of data points over which to do average
    """
    
    # Marker sequence
    markers = ['o', '^']
    
    # Loop through and do analysis for each file
    for m, parameter_object_file in enumerate(parameter_object_file_list):
        
        # Load object
        model_parameters = load_object(parameter_object_file)
        # density files lists (add + '_')
        list_of_peerinf_files = [ str(I) + model_parameters['out_energy_file'] for I in model_parameters['peer_influences'] ]
        
        # read in sets
        global_energy = []
        for energy_file in list_of_peerinf_files:
            global_sets, indiv_sets, total_sets = read_energies(energy_file, model_parameters['npoints'], model_parameters['num_sets'])
            
            # Take mean over last point of all the sets
            global_energy.append(np.mean(zip(*global_sets)[-1]))
        
        # Plot the means
        if xlog==True:
            plt.xscale('log')
        plt.plot(model_parameters['peer_influences'], global_energy, \
                 marker='o', linewidth=1, color='black',\
                 markerfacecolor='none', markersize=12, label='Largest',\
                 markeredgecolor='black', markeredgewidth=1)
    
    plt.xlabel(r'$I$')
    plt.yscale('symlog')
    plt.ylabel(r'$E^{(g)}$')
    plt.savefig(prefix + '-globalE-vs-peerInfluences.pdf',dpi=1200,bbox_inches='tight')
    plt.clf()
    plt.close()

def make_2com_graphs(graph_folder, nprocs, mu_list, net_size, community_size, degree):
    """
    """
    
    ngraphs = len(mu_list)
    # First create # directories equal to number of graphs (not processes)
    folders = [ graph_folder + '/' + str(mu) + '-mu_lfr_graph'  for mu in mu_list ]
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
    
    # next copy benchmark file and parameter file into each folder
    for i, folder in enumerate(folders):
        # Modify benchmark file
        lfr_parameter_file = open('parameters.dat', 'w')
        lines = str(net_size) + '\n'
        lines += str(degree) + '\n'
        lines += str(degree) + '\n'
        lines += '1\n' + '1\n'
        lines += str(mu_list[i]) + '\n'
        lines += str(community_size) + '\n' + str(community_size) + '\n'
        lfr_parameter_file.write(lines)
        lfr_parameter_file.close()
        
        # Copy modified benchmark file
        shutil.copy2('benchmark', folder + '/benchmark')
        shutil.copy2('parameters.dat', folder + '/parameters.dat')
    
    # copy standard parallel code to divide graph making into processes
    def worker(out_q, id, chunksize, folders):
        # worker function gets a list of directors which it must execute
        for folder in folders:
            sub.call('./benchmark', cwd=folder)
            
        outdict = {}
        outdict[id] = id
        out_q.put(outdict)
        
    # No output, just return all processes properly
    chunksize = int(math.ceil(ngraphs / float(nprocs)))
    
    # Get chunks for each process
    folder_chunks = [ folders[chunksize * i:chunksize * (i + 1)] for i in xrange(nprocs) ]
    
    # Start process
    procs = []
    out_q = Queue()
    for i, folder_chunk in enumerate(folder_chunks):
        if len(folder_chunk) > 0:
            p = Process(target=worker, args=(out_q, i, chunksize, folder_chunk))
            procs.append(p)
            p.start()
    
    # Collect successful id's
    resultdict = {}
    for i in range(len(procs)):
        resultdict.update(out_q.get())

    for i in xrange(len(procs)):
        procs[i].join()
    
def make_lfr_graphs(graph_folder, ngraphs, nprocs):
    """
    """
    
    # First create # directories equal to number of graphs (not processes)
    folders = [ graph_folder + '/lfr_graph_' + str(i) for i in xrange(ngraphs) ]
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
    
    # next copy benchmark file and parameter file into each folder
    for folder in folders:
        shutil.copy2('benchmark', folder + '/benchmark')
        shutil.copy2('parameters.dat', folder + '/parameters.dat')
    
    # copy standard parallel code to divide graph making into processes
    def worker(out_q, id, chunksize, folders):
        # worker function gets a list of directors which it must execute
        for folder in folders:
            sub.call('./benchmark', cwd=folder)
            
        outdict = {}
        outdict[id] = id
        out_q.put(outdict)
        
    # No output, just return all processes properly
    chunksize = int(math.ceil(ngraphs / float(nprocs)))
    
    # Get chunks for each process
    folder_chunks = [ folders[chunksize * i:chunksize * (i + 1)] for i in xrange(nprocs) ]
    
    # Start process
    procs = []
    out_q = Queue()
    for i, folder_chunk in enumerate(folder_chunks):
        if len(folder_chunk) > 0:
            p = Process(target=worker, args=(out_q, i, chunksize, folder_chunk))
            procs.append(p)
            p.start()
    
    # Collect successful id's
    resultdict = {}
    for i in range(len(procs)):
        resultdict.update(out_q.get())

    for i in xrange(len(procs)):
        procs[i].join()

def run_multigraph(main_param_dict, nprocs, sets_per_graph):
    """
    """
    
    # Get graph file locations from folder
    folder_list = os.listdir(main_param_dict['graphs_folder'])
    
    main_param_dict['num_sets'] = sets_per_graph * len(folder_list)
    write_cparam_file(main_param_dict['paramfile'], main_param_dict)
    save_object(main_param_dict, main_param_dict['paramfile'] + '.pyobj')
    
    parameter_list = []
    for i, folder in enumerate(folder_list):
        os.chdir(main_param_dict['graphs_folder'] + '/' + folder)
        
        # Setup paramfile
        sub_param_dict = copy.deepcopy(main_param_dict)
        sub_param_dict['graph'] = str(i) + '_' + sub_param_dict['graph']
        sub_param_dict['paramfile'] = str(i) + '_' + sub_param_dict['paramfile']
        sub_param_dict['outfile'] = str(i) + '_' + sub_param_dict['outfile']
        sub_param_dict['sizesfile'] = str(i) + '_' + sub_param_dict['sizesfile'] 
        sub_param_dict['out_energy_file'] = str(i) + '_' + sub_param_dict['out_energy_file']
        sub_param_dict['num_sets'] =  sets_per_graph
        
        # Write graph to file
        graph = read_LFR_Benchmark('network.dat', 'community.dat')
        # REset working directory
        os.chdir('../..')
        
        write_graph_file(sub_param_dict['graph'], graph)
        del graph
        
        parameter_list.append(sub_param_dict)
        
        # Write process parameter file
        write_cparam_file(sub_param_dict['paramfile'], sub_param_dict)
        # Reset any data files
        clear_data_files(sub_param_dict)
        # Run program
        run_influence_analysis(sub_param_dict, sub_param_dict['peer_influences'], nprocs)        
    
    # function that combines the sets of the main param dict from all the files
    
    for i, I in enumerate(main_param_dict['peer_influences']):
        data_files = [ str(I) +'_'+ parameter_dict['sizesfile'] for parameter_dict in parameter_list]

        # Concatenate data files
        concatenate_files(str(I) + '_' + main_param_dict['sizesfile'], data_files)

def run_multigraph_mu_analysis(prefix, basic_model_parameters, nprocesses, ngraphs, point_range=(195,200), avg_degree=5):
    """
    """

    # Save model parameter object
    save_object(basic_model_parameters, basic_model_parameters['paramfile'] + '.pyobj')

    # Generate graphs (temporarily)
    graph_dir_list = []
    for g in xrange(ngraphs):
        graph_dir = prefix + '_lfr_2com_graphs_temp_' + str(g)
        make_2com_graphs(graph_dir, nprocesses, basic_model_parameters['mu_list'], basic_model_parameters['N'], \
            basic_model_parameters['N']/2, avg_degree)
        graph_dir_list.append(graph_dir)

    # Run for each graph
    list_of_tsmeans_for_each_mu_for_each_graph = []
    for g in xrange(ngraphs):
        # Get graph file locations from folder
        folder_list = os.listdir(graph_dir_list[g])
        model_parameter_file_list = []
        for i, mu in enumerate(basic_model_parameters['mu_list']):
            os.chdir(graph_dir_list[g] + '/' + str(mu) + '-mu_lfr_graph')
            
            model_parameters = copy.deepcopy(basic_model_parameters)
            model_parameters['graph'] = str(mu) + '-mu_' + model_parameters['graph']
            model_parameters['paramfile'] = str(mu) + '-mu_' + model_parameters['paramfile']
            model_parameters['outfile'] = str(mu) + '-mu_' + model_parameters['outfile']
            model_parameters['sizesfile'] = str(mu) + '-mu_' + model_parameters['sizesfile']
            model_parameters['out_energy_file'] = str(mu) + '-mu_' + model_parameters['out_energy_file']
            model_parameters['mu'] = mu
            
            # Write graph to file
            graph = read_LFR_Benchmark('network.dat', 'community.dat')
            # REset working directory
            os.chdir('../..')
            
            if model_parameters['assignment'] == True:
                graph = relabel_graph(graph, model_parameters['assignment_dict'])
                model_parameters['assignments'] = community_assignment(graph, model_parameters['assignment_dict'])
            write_graph_file(model_parameters['graph'], graph)
            
            del graph
            
            # Write process parameter file
            write_cparam_file(model_parameters['paramfile'], model_parameters)
            save_object(model_parameters, model_parameters['paramfile'] + '.pyobj')
            model_parameter_file_list.append(model_parameters['paramfile'] + '.pyobj')
            # Reset any data files
            clear_data_files(model_parameters)
            # Run program
            run_parallel(nprocesses, model_parameters)

        list_of_mu_files = [ str(mu) + '-mu_' + basic_model_parameters['sizesfile'] for mu in basic_model_parameters['mu_list'] ]
        # Loop through each parameter file
        list_of_tsmeans_for_each_mu = []
        for mu_file in list_of_mu_files:
            sets = read_sizes_output(mu_file, basic_model_parameters['npoints'])
            
            # Calculate average over all sets
            largest_groups = []
            for j, set in enumerate(sets):
                # Get zealot time-series
                largest = Normalize([ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in set ], 0.5, 1.0)                
                largest_groups.append(largest)
            
            largest_groups_tsmean = [ np.mean(set_ts[point_range[0]: point_range[1]]) for set_ts in largest_groups ]
            list_of_tsmeans_for_each_mu.append(largest_groups_tsmean)

        list_of_tsmeans_for_each_mu_for_each_graph.append(list_of_tsmeans_for_each_mu)

    list_of_tsmeans_for_each_graph_for_each_mu = zip(*list_of_tsmeans_for_each_mu_for_each_graph)
    list_of_tsmeans_for_each_mu_graphs_folded = [ [ item for subsublist in sublist for item in subsublist] for sublist in list_of_tsmeans_for_each_graph_for_each_mu  ]
    save_object((basic_model_parameters['mu_list'], list_of_tsmeans_for_each_mu_graphs_folded), prefix + '_mu_multigraph.pltdat')

def plot_mu_pltdat(prefix, pltdat_files, labels):
    # Marker sequence
    markers = ['o', '<','s','D']
    colors = ['darkblue', 'royalblue','skyblue']
    # create zealot density plotter using saved graph data
    for i, pltdat_file in enumerate(pltdat_files):
        # load data for file
        mus, list_of_tsmeans_for_each_graph_for_each_mu = load_object(pltdat_file)

        # calculate mean
        mean_of_largest = np.mean(list_of_tsmeans_for_each_graph_for_each_mu, 1)

        # calculate errors
        std_of_largest = np.std(list_of_tsmeans_for_each_graph_for_each_mu, 1)
        mean_err = stats.sem(list_of_tsmeans_for_each_graph_for_each_mu, 1)
        top_err = mean_err
        bot_err = mean_err
        top_err = [ 1.0-mean_of_largest[j] if (val+mean_of_largest[j]) > 1.0 else val for j, val in enumerate(std_of_largest) ]# Trim top error bar since you can't have more than 1.0 
        bot_err = [ mean_of_largest[j] if (mean_of_largest[j]-val) < 0.0 else val for j, val in enumerate(std_of_largest) ] # Trim bottom, since nothing is below 0

        # plot w/ some set of shapes and labels
        plt.errorbar(mus, mean_of_largest, yerr=[bot_err, top_err],\
                 marker=markers[i], linewidth=1, color=colors[i], \
                 markerfacecolor='none', markersize=8, label=labels[i],
                 markeredgecolor=colors[i], markeredgewidth=1)

    # Write plots out to file
    plt.xlim(min(mus), max(mus))
    plt.ylim(0.0, 1.02)
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'$\textup{Fraction converted}$')
    plt.legend(loc=4, frameon=False)
    plt.tight_layout()
    plt.savefig(prefix + '-community_rigidity.pdf',dpi=1200,bbox_inches='tight')
    plt.clf()
    plt.close()

def find_critical_I_pll(prefix, nprocs, model_parameters, Iest, num_graphs=20, num_IC=20, haulting_precision=0.001):
    """
    Iest is a tuple with lower bound on est and upper bound on est
    """

    def worker(out_q, id, parameters, Iest, graph_chunk, num_IC=20, haulting_precision=0.001, rndgraph=True):
        """
        """

        # Problem with different sized graphs each time, won't be able to calculate ratio

        # Write process parameter file
        write_cparam_file(parameters['paramfile'], parameters)
        # Reset any data files
        clear_data_files(parameters)
        
        # For different graphs (connected graphs)
        graph_Ic_list = []
        
        # Loop length
        if rndgraph == True:
            num_loops = graph_chunk
        else:
            num_loops = len(graph_chunk)
        
        for i in xrange(num_loops):
            
            if rndgraph == True:
                # Generate new graph -> take connected components only
                graph = relabel_nodes(get_connected_subgraph(nx.fast_gnp_random_graph(parameters['rndGenPar'][0], parameters['rndGenPar'][1])))

                parameters['N'] = len(graph)
                # Save graph with name from model parameters
                write_graph_file(parameters['graph'], graph)
            else:
                # Read in graph from file
                graph = read_LFR_Benchmark(parameters['graphs_folder'] + '/' + graph_chunk[i] +'/network.dat', parameters['graphs_folder'] + '/' + graph_chunk[i] +'/community.dat')
                write_graph_file(parameters['graph'], graph)
            
            # Loop through each randomly generated initial condition (different seed)
            initCond_Ic_list = [] 
            
            for j in xrange(num_IC):
                # Right now seed is just j + 1. It can be changed if you want, as long as it changes each time
                seed = j + 1

                # Update seed
                parameters['seed'] = seed
    
                # Reset Ic variables
                low_bound = Iest[0]
                up_bound = Iest[1]
                Ic = (low_bound + up_bound) / 2.
                
                # Repeat till desired presicion
                percolates = False
                while  (up_bound - low_bound) > haulting_precision:
                    
                    # Update parameter file
                    parameters['I'] = Ic
                    write_cparam_file(parameters['paramfile'], parameters)
                    # Clear old datafiles
                    clear_data_files(parameters)
                    # Run for X timesteps
                    run_cpp(parameters['cpp_exe'], parameters['paramfile'])
                    
                    # Read sizes file
                    sets = read_sizes_output(parameters['sizesfile'], parameters['npoints'])
                    # Take largest in set from last data point and average over all sets
                    largest = np.mean([ max(set[-1])[0] / float(parameters['N']) for set in sets])
                    
                    # If netsized component
                    if largest >= 1.0:
                        percolates = True
                        up_bound = Ic
                    # If no netsized component
                    else:
                        low_bound = Ic
                        
                    Ic = (low_bound + up_bound) / 2.
    
                if percolates == True:
                    initCond_Ic_list.append(Ic)
            
            graph_Ic_list.append(initCond_Ic_list)
    
        out_q.put({id: graph_Ic_list})

    # Create chunks for each process to handle
    chunksize = int(math.ceil(num_graphs / float(nprocs)))
    
    # Either use pre-set graphs or randomly generated ones
    if model_parameters['rndGenPar'] != None:
        rndgraph = True
        graph_chunks = [ chunksize for i in xrange(nprocs) ]
    else:
        rndgraph = False
        # Get the graph folders from the graph folder directory
        folder_list = os.listdir(model_parameters['graphs_folder'])
        # Use ony num_graphs of those available
        folder_list = folder_list[:num_graphs]
        # Assign these graphs to each process
        graph_chunks = [ folder_list[chunksize * i:chunksize * (i + 1)] for i in xrange(nprocs) ]
         
    # generate model parameters for each process
    model_param_list = []
    for i in xrange(nprocs):
        chunk_params = copy.deepcopy(model_parameters)
        chunk_params['paramfile'] = str(i) + '_' + model_parameters['paramfile']
        chunk_params['outfile'] = str(i) + '_' + model_parameters['outfile']
        chunk_params['sizesfile'] = str(i) + '_' + model_parameters['sizesfile']
        chunk_params['out_energy_file'] = str(i) + '_' + model_parameters['out_energy_file']
        chunk_params['graph'] = str(i) + '_' + model_parameters['graph']
        
        model_param_list.append(chunk_params)

    # Start process
    procs = []
    out_q = Queue()
    for i, params in enumerate(model_param_list):
        p = Process(target=worker, args=(out_q, i, params, Iest, graph_chunks[i], num_IC, haulting_precision, rndgraph))
        procs.append(p)
        p.start()

    # Collect successful id's
    resultdict = {}
    for i in range(len(procs)):
        resultdict.update(out_q.get())

    for i in xrange(len(procs)):
        procs[i].join()

    # Take data from processes and combine all the lists into one flat list
    combined_Ic_list = list(itertools.chain(*(resultdict.values())))

    # Save data to file
    save_object(combined_Ic_list, prefix + '_Ic_lists.pyobj')
     
    # Take averages (by flattening lists and taking averages)
    IcList = [ ICIc for graphIc_list in combined_Ic_list \
              for ICIc in graphIc_list]
    best_Ic_mean, dump, dump = stats.bayes_mvs(IcList, 0.9)

    return best_Ic_mean, IcList

def load_outfile_to_graph(outfile, graph):
    """
    works only for nodes that start at 0 and end at N-1
    """

    with open(outfile, 'r') as f:
        node_counter = 0
        for line in f:
            if node_counter == len(graph):
                node_counter = 0

            split_line = line.split(',')
            beliefs = map(int, split_line[:-1])
            energy = float(split_line[-1])
            try:
                graph.node[node_counter]['beliefs'].append(beliefs)
                graph.node[node_counter]['energys'].append(energy)
            except KeyError:
                graph.node[node_counter]['beliefs'] = []
                graph.node[node_counter]['energys'] = []
                graph.node[node_counter]['beliefs'].append(beliefs)
                graph.node[node_counter]['energys'].append(energy)
            
            node_counter += 1

def find_s_clusters(graph):
    """
    """

    # s_clusters_dict = {}
    s_cluster_dist = []
    for i in xrange(len(graph.node[0]['beliefs'])):
        # Find unique belief types:
        unique_beliefs = {}
        for j in graph.nodes_iter():
            try:
                unique_beliefs[tuple(graph.node[j]['beliefs'][i])].append(j)
            except:
                unique_beliefs[tuple(graph.node[j]['beliefs'][i])] = []
                unique_beliefs[tuple(graph.node[j]['beliefs'][i])].append(j)

        # Brake graph into subgraphs based on beliefs
        for nodes in unique_beliefs.values():
            belief_subgraph = graph.subgraph(nodes)
            belief_groups = nx.connected_component_subgraphs(belief_subgraph)
            for sub_group in belief_groups:
                s_cluster_dist.append(len(sub_group) / float(len(graph)))

    return s_cluster_dist

def critical_exp_analysis(outfile, graph_file, paramfile):
    """
    """

    # load graph
    graph = load_object(graph_file)

    # Load outfile
    load_outfile_to_graph(outfile, graph)

    # Load paramfile
    params = load_object(paramfile)

    # Keys are integers of the size of the group, values are the number of such clusters
    # s_clusters_dict = {}

    # Get s_clusters_dist
    s_clusters_dist = find_s_clusters(graph)

    return s_clusters_dist

def plot_size_distribution(prefix, s_clusters_dist):
    """
    """

    s_clusters_dict = Counter(s_clusters_dist)
    inv_cum_freq = []
    sizes = []
    s_clusters_list = sorted(s_clusters_dict.items())
    total_freq = np.sum(zip(s_clusters_list))
    for i, item in enumerate(s_clusters_list):
        cum_freq = 0.0
        for size, freq in s_clusters_list[i:]:
            cum_freq += freq
        inv_cum_freq.append(cum_freq / float(total_freq))
        sizes.append(item[0])
    ax = plt.subplot(111)
    plt.loglog(sizes, inv_cum_freq, marker='o', color='blue', linestyle='none')

    slope, inter, r_val, p_val, stderr = stats.linregress(np.log(sizes), np.log(inv_cum_freq))
    linreg_x = np.array([min(sizes), max(sizes)])
    linreg_y = np.exp(inter) * linreg_x**slope
   
    plt.loglog(linreg_x, linreg_y, 'r--')
    alpha_label = r"$\alpha=" +str(round(slope, 3)) + "$"
    plt.text(0.1,0.15, alpha_label,fontsize=20,transform=ax.transAxes)
    r_squared_label = r"$r^2="+str(round(r_val**2,3))+"$"
    plt.text(0.1,0.07, r_squared_label,fontsize=20,transform=ax.transAxes)
    plt.ylabel(r'$P(S\geq s)$')
    plt.xlabel(r'$s/N$')
   
    plt.savefig(prefix + '_size_distribution.pdf')
    plt.clf()
    plt.close()

def get_connected_subgraph(graph):
    """
    """
    
    subgraph = nx.connected_component_subgraphs(graph)[0]
    new_labels = range(len(subgraph))
    rnd.shuffle(new_labels)
    label_dict = {}
    for i, node in enumerate(subgraph.nodes()):
        label_dict[node] = new_labels[i]
        
    return nx.relabel_nodes(subgraph, label_dict, copy=True)

def figure5(prefix, contour_left_data_file, contour_right_data_file):
    """
    """

    point_range=(195, 200) 
    # Load file list
    param_file_list = load_object(contour_left_data_file)

    # Loop through each parameter file and get parameter values
    mean_largest_sizes = []
    var1_list = []
    var2_list = []
    for param_file in param_file_list:
        # Load parameters and load size data
        model_parameters = load_object(param_file)
        var1_list.append(model_parameters['mu'])
        var2_list.append(model_parameters['J'])

        # Calculate average of all sets
        sets = read_sizes_output(model_parameters['sizesfile'], model_parameters['npoints'])
        largest_groups = []
        for j, output_set in enumerate(sets):
            # Get zealot time-series
            largest = Normalize([ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in output_set ], 0.5, 1.0)
            largest_groups.append(largest)
        
        largest_groups = zip(*largest_groups)
        # Take mean over all sets
        largest_means = [ np.mean(largest_groups[i]) for i in xrange(point_range[0], point_range[1])]
        # Take mean over range (assumes stationary average)
        largest_mean = np.mean(largest_means)
        mean_largest_sizes.append(largest_mean)

    # Combine into gridspace
    var1_rangeL = sorted(list(set(var1_list)))
    var2_rangeL = sorted(list(set(var2_list)))
    vv1L, vv2L = np.meshgrid(var1_rangeL, var2_rangeL, indexing='ij')
    zz_largeL = np.zeros(vv1L.shape)
    for iii in xrange(len(var1_rangeL)):
        for jjj in xrange(len(var2_rangeL)):
            if (var1_list[iii * len(var2_rangeL) + jjj] == var1_rangeL[iii]) and (var2_list[iii * len(var2_rangeL) + jjj] == var2_rangeL[jjj]):
                zz_largeL[iii, jjj] = mean_largest_sizes[ iii * len(var2_rangeL) + jjj ]



    # Load file list
    param_file_list = load_object(contour_right_data_file)

    # Loop through each parameter file and get parameter values
    mean_largest_sizes = []
    var1_list = []
    var2_list = []
    for param_file in param_file_list:
        # Load parameters and load size data
        model_parameters = load_object(param_file)
        var1_list.append(model_parameters['mu'])
        var2_list.append(model_parameters['J'])

        # Calculate average of all sets
        sets = read_sizes_output(model_parameters['sizesfile'], model_parameters['npoints'])
        largest_groups = []
        for j, output_set in enumerate(sets):
            # Get zealot time-series
            largest = Normalize([ sorted(point, reverse=True)[0][0] / float(model_parameters['N']) for point in output_set ], 0.5, 1.0)
            largest_groups.append(largest)
        
        largest_groups = zip(*largest_groups)
        # Take mean over all sets
        largest_means = [ np.mean(largest_groups[i]) for i in xrange(point_range[0], point_range[1])]
        # Take mean over range (assumes stationary average)
        largest_mean = np.mean(largest_means)
        mean_largest_sizes.append(largest_mean)

    # Combine into gridspace
    var1_rangeR = sorted(list(set(var1_list)))
    var2_rangeR = sorted(list(set(var2_list)))
    vv1R, vv2R = np.meshgrid(var1_rangeR, var2_rangeR, indexing='ij')
    zz_largeR = np.zeros(vv1R.shape)
    for iii in xrange(len(var1_rangeR)):
        for jjj in xrange(len(var2_rangeR)):
            if (var1_list[iii * len(var2_rangeR) + jjj] == var1_rangeR[iii]) and (var2_list[iii * len(var2_rangeR) + jjj] == var2_rangeR[jjj]):
                zz_largeR[iii, jjj] = mean_largest_sizes[ iii * len(var2_rangeR) + jjj ]


    # PLOT
    # f, axarr = plt.subplots(1, 2, sharey=True)
    f = plt.figure()
    levels = np.linspace(0.0, 1.0001, 6)
    plt.contourf(vv1L, vv2L, zz_largeL, levels, cmap=cm.PuBu)
    plt.ylabel(r'$J$')
    plt.xlabel(r'$\mu$')
    plt.text(0.05,9,'$I=0.001$', size=28) #adjust as needed

    plt.tight_layout()
    plt.savefig(prefix + "_combo-mu-contour1.pdf")
    plt.clf()
    plt.close()

    plt.figure(figsize=(8.5,6))
    CS = plt.contourf(vv1R, vv2R, zz_largeR, levels, cmap=cm.PuBu)
    plt.xlabel(r'$\mu$')
    bar = plt.colorbar(CS)
    bar.ax.set_ylabel(r'$\textup{Fraction converted}$')
    plt.text(0.05,9,'$I=10.0$', size=28) # adjust as needed

    plt.tight_layout()
    plt.savefig(prefix + "_combo-mu-contour2.pdf")
    plt.clf()
    plt.close()


##------------------------------------------------------------------------------
##
##    Main
##
##------------------------------------------------------------------------------
if __name__ == '__main__':
    """
    If a node has limited # inputs then depending on the parameter values
    chosen, there will be a point at which a node CAN never accept a new edge
    if it has the contradicting edge, even if all his neighbors are sending him
    the contradiction, because the amount subtracted by the contradiction will
    be strictly greater than any neighborly emissions. This will dependon the
    size of the memory. Larger memory will eliminate this issue.
    
    Need term for initial edge ratios w/o setting them immutable.
    """
    
    #===========================================================================
    # Plot settings
    #===========================================================================
    params = {'backend': 'ps', \
              'axes.labelsize': 28, \
              'text.fontsize': 20, \
              'legend.fontsize': 20, \
              'xtick.labelsize': 24, \
              'ytick.labelsize': 24, \
              'text.usetex': True, \
              'xtick.major.size': 10, \
              'ytick.major.size': 10, \
              'xtick.minor.size': 8, \
              'ytick.minor.size': 8, \
              'xtick.major.width': 1.0, \
              'ytick.major.width': 1.0, \
              'xtick.minor.width': 1.0, \
              'ytick.minor.width': 1.0}
    plt.rcParams.update(params)

    #===========================================================================
    # Create cached file
    #===========================================================================
#     filelist = ['coh_file3.dat', 'coh_file4.dat', 'coh_file5.dat', 
#                 'coh_file6.dat', 'coh_file7.dat' ]
#     write_csv_cohfiles(filelist)
#     sys.exit()

    # Create voter model cohfile
#     write_voter_energyfile(5)

   
    #===========================================================================
    # Model setup
    #===========================================================================
 
#     # Pick a graph
    # graph = nx.convert_node_labels_to_integers(nx.grid_2d_graph(20,20,periodic=True))
    # graph = read_LFR_Benchmark('network.dat', 'community.dat')
    graph = get_connected_subgraph(nx.fast_gnp_random_graph(1000, 0.005)) #200, 0.025# 100,000 - 0.00005 (avg 5) # 10,000 - .0005 (avg5)

    # Parameter file
    graph_size = len(graph)
    # graph_size = 1000 #10000
    graphs_folder = 'lfr_2com_graphs_test' #'lfr_2com_graphs_nonlog' #
    zealots = False #True
    random = False #False
    assignment = False #False
    assignment_dict = {1: 'b0', 2: 'b1'}#{1: 'b0', 2: 'z0'}
    num_sets = 160#*160
    param_values = log_linspace(0.001, 10.0, 90)#log_linspace(0.001, 10.0, 60) #log_linspace(0.001, 10.0, 90)
    densities = np.arange(0.003, 0.15, 0.003) # for zealot figure
    mus = np.linspace(0.01, 0.5, 30) #np.linspace(0.01, 0.5, 30) #log_linspace(0.001, 0.5, 30)
    
    # if (random == True) and (zealots == False):
    #     num_sets = 0

    if not ((random == True) and (zealots == False)):
        belief_set_parameters = {'coh_file':'coh_file4.dat',
                                 'zealots': zealots,
                                 'zealot_energy_ranges':[(-1,-1)],
                                 'belief_energy_ranges':[(1,1), (-1,-1)],
                                 'num_sets': num_sets
                                 } # (-0.8,-0.3) (-0.29, -0.1) (-0.08, 0.0) (-1,-1)
         
        belief_system_sets, zealot_sets, num_sets = generate_beliefsys_sets(**belief_set_parameters)

    else:
        belief_system_sets = []
        zealot_sets = []
     
    # Generate assignments if flag is set
    if assignment == True and random == False:
        assignment_list = None
#         graph = relabel_graph(graph, assignment_dict)
#         assignment_list = community_assignment(graph, assignment_dict)
    elif assignment == True and random == True:
        print 'not valid parameter combination'
        sys.exit()
    else:
        assignment_list = []

 
    # Very important:
    # The sum of the ratios of belief_types and zealot_types should be 1.
    # If they exceed 1., then the last belief_type will get crowded out,
    # resulting in a smaller (perhaps 0) ratio of that type.
    # If they are less than 1., then the last belief_type will get enhanced,
    # resulting in a larger ratio of that type.
    prefix = 'unstable_collapse'#"two_var_I_T_J1.0_N1000_Coh4"#"two_var_J_I_T0.001_N1000_Coh4" #'param_analysis_N1000_J_T0.001_I0.1_Coh4'
    model_parameters = {'N': graph_size,
                        'multigraph_p': 0.005,
                        'paramfile': prefix + '_params.dat',
                        'seed': 10, #
                        'outfile': prefix + '_out.dat',
                        'sizesfile': prefix + '_size.dat',
                        'out_energy_file': prefix + '_energy.dat',
                        'T': 2.0, # 0.001
                        'I': 0.001, # 0.1
                        'J': 2.0, # 1.0
                        'cpp_exe': "BeliefNet",
                        'mu': None,
                        'runs': 1,
                        'energy_file': 'coh6_energy.cdat',
                        'update': 'edge',
                        'npoints': 200,#200,
                        'niterations': 200000,#1000000000,#150000000, #50,000,000-150,000,000
                        'random': random,
                        'zealots': zealots,
                        'num_sets': num_sets,
                        'num_belief_types': 2,
                        'belief_types': belief_system_sets,
                        'belief_ratios': [0.99, 0.01], # This must add with zealot ratios to 1.0
                        'num_zealot_types': 1,
                        'zealot_types': zealot_sets,
                        'zealot_ratios': [0.5], # This must add with belief ratios to 1.0
                        'graph':prefix + '-graph.graph',
                        'outfile_prefix': prefix,
                        'densities': densities,
                        'param_values': param_values,
                        'num_nodes': graph_size,
                        'graphs_folder': graphs_folder,
                        'assignment': assignment, #Boolean
                        'assignments': assignment_list,
                        'assignment_dict': assignment_dict,
                        'mu_list': mus,
                        'rndGenPar': (10000,0.0005),
                        'timeseries': True,
                        'ts_sizes': True, #True
                        'finalstate': False, #False,
                        'AppendOut': False
                        }  #(200, 0.05) }#None} 
    # write_cparam_file(model_parameters['paramfile'], model_parameters)
    # save_object(model_parameters, model_parameters['paramfile'] + '.pyobj')
    # save_object(graph, model_parameters['graph'] + '.pyobj')
    # write_graph_file(model_parameters['graph'], graph)
    # nx.write_gexf(graph, 'test_com_graph.gexf')
#     del graph

    # Clear existing data files
    # clear_data_files(model_parameters)
     
#     sys.exit()
    #===========================================================================
    # Run Model
    #===========================================================================

    # function to run c++ program given the name of the program
    # run_cpp(model_parameters['cpp_exe'], model_parameters['paramfile'])
    # Run sets in parallel (assumes graph file already constructed)
    # run_parallel(20, model_parameters)

    # Run zealot analysis over many densities
    # run_zealot_density_analysis(model_parameters, densities, 3)

    # run_multigraph_zealot_density_analysis(prefix, model_parameters, densities, nprocesses=20, ngraphs=20)

    # Run peer influence analysis
    # run_param_analysis(model_parameters, param_values, 'J', 20)

    # Run mu analysis
    # run_mu_analysis(model_parameters, 3)

    # run_multigraph_mu_analysis(prefix, model_parameters, nprocesses=20, ngraphs=20)

    # Run two parameter simulation (MU)
    # mu_range = log_linspace(0.001, 0.5, 30)
    # T_range = log_linspace(0.001, 10.0, 5)
    # two_variable_mu_analysis(model_parameters, 2, ("mu", "J"), mu_range, T_range)

    # Run two parameter simulation
    # v1_range = log_linspace(0.001, 10.0, 30)
    # v2_range = log_linspace(0.001, 10.0, 30)
    # two_variable_analysis(model_parameters, 20, ("I", "T"), v1_range, v2_range)

    #===========================================================================
    # Make pretty graphs
    #===========================================================================

    # plot part of figure 5
    # figure5("figure5-contourpart", 'zealinvade_N1000_coh4_contour_e-1_I0.001_T1.0_params.dat_filelist.pyobj', \
        # 'zealinvade2_N1000_coh4_contour_e-1_I10.0_T1.0_params.dat_filelist.pyobj')

    # Plot multigraph mu
    # plot_mu_pltdat("combo", ['zealot-invade_N1000_e-1_T2.0_I1.0_J2.0_Coh4_mu_multigraph.pltdat',\
    #     'zealot-invade_N1000_e0_T2.0_I1.0_J2.0_Coh4_mu_multigraph.pltdat', \
    #     'zealot-invade_N1000_e1_T2.0_I1.0_J2.0_Coh4_mu_multigraph.pltdat'], ['$E_z=-1.0$','$E_z=0.0$','$E_z=1.0$'])

    # plot_mu_pltdat("combo", ['zealot-invade_N1000_e-1_T2.0_I0.001_J2.0_Coh4_mu_multigraph.pltdat', \
    #     'zealot-invade_N1000_e0_T2.0_I0.001_J2.0_Coh4_mu_multigraph.pltdat',\
    #     'zealot-invade_N1000_e1_T2.0_I1.0_J2.0_Coh4_mu_multigraph.pltdat'], ['$E_z=-1.0$','$E_z=0.0$','$E_z=1.0$'])

    # Plot multigraph
    # plot_density_pltdat("combo", ['zealot-density_N1000_e-1_T2.0_I1.0_J2.0_Coh4' + '_multigraph.pltdat' ,\
    #     'zealot-density_N1000_e0_T2.0_I1.0_J2.0_Coh4' + '_multigraph.pltdat', \
    #     'zealot-density_N1000_e1_T2.0_I1.0_J2.0_Coh4' + '_multigraph.pltdat'], \
    #     ['$E_z=-1.0$ (coherent)', '$E_z=0.0$', '$E_z=1.0$ (incoherent)'])

    # Make 2-variable mu contour
    # prefix = "mu-J_I0.001_T1.0"
    # plot_two_var_contour(prefix, prefix + '_params.dat' + '_filelist.pyobj', (195,200), 'I', 'T', logvar1=True, logvar2=True)
    # plot_two_var_contour("two_var_I_T_J1.0_N1000_Coh4", 'two_var_I_T_J1.0_N1000_Coh4_params.dat' + '_filelist.pyobj', (195,200), 'I', 'T', logvar1=True, logvar2=True, switch=True)
    # plot_two_var_contour("two_var_J_T_I0.001_N1000_Coh4", 'two_var_J_T_I0.001_N1000_Coh4_params.dat' + '_filelist.pyobj', (195,200), 'J', 'T', logvar1=True, logvar2=True, switch=True)

    # Make 3rd figure of paper
    # energies_plot_pubversion('unstable_collapse', 'unstable_collapse_N1000_T2.0_I0.001_J2.0_coh4_params.dat.pyobj')
    # sizes_plot_pubversion('unstable_collapse', 'unstable_collapse_N1000_T2.0_I0.001_J2.0_coh4_params.dat.pyobj')
    unstable_collapes_1col('unstable_collapse', 'unstable_collapse_N1000_T2.0_I0.001_J2.0_coh4_params.dat.pyobj')
    # unstable_collapes_1col('unstable_collapse', 'unstable_collapse_params.dat.pyobj')
    
    # Make sizes plot
    # plot_sizes(prefix, model_parameters)
    # plot_energies(prefix, prefix + '_params.dat.pyobj')

    # Make energy plot
    # prefix = "unstable_collapse"
    # energies_plot_pubversion(prefix, prefix + '_params.dat.pyobj')

    # Make density of figure
    # zealot_densities_plot(prefix, model_parameters, (195,200))
    
    # Make figure 4 of paper
    # zealot_densities_plot_pubversion('zealot_density', (195,200), ['cohzeal-i2_params.dat.pyobj', 'cohzeal-i1.0_params.dat.pyobj', \
    #     'cohzeal-i0.5_params.dat.pyobj', 'hizeal-i2_params.dat.pyobj'], \
    #  [r'Coherent Zealots $I=2.0$', r'Coherent Zealots $I=1.0$', \
     # r'Coherent Zealots $I=0.5$', r'Incoherent zealots $I=2.0$' ])
    # zealot_densities_plot_pubversion('zealot_density', (198,200), ['cohzeal-i2_params.dat.pyobj', 'cohzeal-i0.5_params.dat.pyobj', 'cohzeal-i0.2_params.dat.pyobj', 'hizeal-i2_params.dat.pyobj'], [r'Coherent Zealots $I=2.0$', r'Coherent Zealots $I=0.5$', r'Coherent Zealots $I=0.2$', r'Incoherent zealots $I=2.0$' ])
    # zealot_densities_plot_pubversion('zealot_density', (175,200), ['lowzealot_params.dat.pyobj', 'cohzeal-ilow_params.dat.pyobj' ,'highzealot_params.dat.pyobj'], ['Coherent zealots','Low I Coherence zealots','Incoherent zealots'])

    # Make peer influence figure
    # plot_param_analysis(prefix, model_parameters, 'J', (195,200))

    # Make figure 2
    # param_plot_pubversion(prefix, 'J', (195,200), [prefix + '_params.dat.pyobj'], plot_second=True)
    # heatmap_param_pubversion(prefix, 'T', prefix + '_params.dat.pyobj')


    # Make global energy plot
    # global_energy_v_I('rnd', ['Irnd_params.dat.pyobj'])

    # Make mu figure
    # plot_size_vs_mu(prefix, prefix + '_params.dat_filelist.pyobj', (195,200), True)
    # plot_size_vs_mu('mu_stableZ_stableB', 'mu_stableZeal_stableBelief_params.dat_filelist.pyobj', (170,200), True)
    # plot_second_size_vs_mu('combo', ['zealinvade_N1000_e-1_T2.0_I0.001_J2.0_Coh4_params.dat_filelist.pyobj', \
    #     'zealinvade_N1000_e0_T2.0_I0.001_J2.0_Coh4_params.dat_filelist.pyobj'\
    #     ,'zealinvade_N1000_e1_T2.0_I0.001_J2.0_Coh4_params.dat_filelist.pyobj'], ['$E_z=-1$','$E_z=0$','$E_z=1$'], (195,200), False)

    # Make critical exponent analysis
    # plot_size_distribution('rnd_crit2', critical_exp_analysis('rnd_crit_out.dat', 'graph-rnd_crit.graph.pyobj', 'rnd_crit_params.dat.pyobj'))

    # Investigate peer influence analysis. >> plot all 20 sets for an I that has max 1.0 for LFR
#     model_parameters = load_object('Irnd_params.dat.pyobj')
#     file_list = [ datafile for datafile in os.listdir(os.getcwd()) if 'Irnd_size.dat' in datafile ]
#     for file in file_list:
#         print file
#     
#     for i, file in enumerate(file_list):
#         plot_sizes(file, model_parameters, file)
    
    # make graphs
#     make_lfr_graphs('lfr_graphs', 20, 20)

    # multi-graph
#     run_multigraph(model_parameters, 20, 1)
#     peer_influences_plot('test', load_object('test_params.dat.pyobj'), (150,200))

    # Load Ic results
#     ic_results = load_object('ic_rnd_Ic_lists.pyobj')
#     IcList = [ ICIc for graphIc_list in ic_results \
#               for ICIc in graphIc_list]
#     bestIc, dump, dump = stats.bayes_mvs(IcList, 0.9)
#     print 'Ic:',bestIc[0], '\tL:',bestIc[0]-bestIc[1][0], 'H:',bestIc[1][1]-bestIc[0]
    
#     sys.exit()