import os
import sys
import pandas as pd
import numpy as np
from scipy.special import gamma,psi
from sklearn.neighbors import NearestNeighbors
from scipy import special, spatial
from scipy.integrate import odeint
from scipy import stats, linalg
from sklearn.metrics import mutual_info_score
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nitime
from nitime.timeseries import TimeSeries
from nitime.analysis import CorrelationAnalyzer, CoherenceAnalyzer
from nitime.utils import percent_change
from nitime.viz import drawmatrix_channels, drawgraph_channels, plot_xcorr
from rbig import RBIG, RBIGMI
sys.path.insert(0, './bio_corex')
import corex as ce
import matplotlib.pyplot as plt
import community
import networkx as nx
from community import community_louvain
from nxpd import nxpdParams, draw
#plt.style.use('seaborn-white')


#code modified from https://github.com/cryptofan/FunctionalConnectivity_AoNBrainhackWarsaw/blob/master/codes/kraskov_mutual_info.ipynb

def nearest_distances(X, k=1):
    '''
    X = array(N,M)
    N = number of points
    M = number of dimensions

    returns the distance to the kth nearest neighbor for every point in X
    '''
    knn = NearestNeighbors(n_neighbors=k)
    knn.fit(X)
    d, _ = knn.kneighbors(X) # the first nearest neighbor is itself
    return d[:, -1] # returns the distance to the kth nearest neighbor


def entropy(X, k=1):
    ''' Returns the entropy of the X.

    Parameters
    ===========

    X : array-like, shape (n_samples, n_features)
        The data the entropy of which is computed

    k : int, optional
        number of nearest neighbors for density estimation

    Notes
    ======

    Kozachenko, L. F. & Leonenko, N. N. 1987 Sample estimate of entropy
    of a random vector. Probl. Inf. Transm. 23, 95-101.
    See also: Evans, D. 2008 A computationally efficient estimator for
    mutual information, Proc. R. Soc. A 464 (2093), 1203-1215.
    and:
    Kraskov A, Stogbauer H, Grassberger P. (2004). Estimating mutual
    information. Phys Rev E 69(6 Pt 2):066138.
    '''

    # Distance to kth nearest neighbor
    r = nearest_distances(X, k) # squared distances
    n, d = X.shape
    volume_unit_ball = (np.pi**(.5*d)) / gamma(.5*d + 1)
    '''
    F. Perez-Cruz, (2008). Estimation of Information Theoretic Measures
    for Continuous Random Variables. Advances in Neural Information
    Processing Systems 21 (NIPS). Vancouver (Canada), December.

    return d*mean(log(r))+log(volume_unit_ball)+log(n-1)-log(k)
    '''
    return (d*np.mean(np.log(r + np.finfo(X.dtype).eps))
            + np.log(volume_unit_ball) + psi(n) - psi(k))

def mutual_information(variables, k=1):
    '''
    Returns the mutual information between any number of variables.
    Each variable is a matrix X = array(n_samples, n_features)
    where
      n = number of samples
      dx,dy = number of dimensions
    Optionally, the following keyword argument can be specified:
      k = number of nearest neighbors for density estimation
    Example: mutual_information((X, Y)), mutual_information((X, Y, Z), k=5)
    '''
    if len(variables) < 2:
        raise AttributeError(
                "Mutual information must involve at least 2 variables")
    all_vars = np.hstack(variables)
    
    MI = sum([entropy(X, k=k) for X in variables]) - entropy(all_vars, k=k)
    MI = np.log(2**MI)
    
    return MI


# compute FUNCTIONAL CONNECTIVITY on the data with KRASKOV MUTUAL INFORMATION:        
# you need to specify the number of neighbours

def KMI(data, k):
    matMI = np.zeros((Nvars,Nvars))

    for ix in np.arange(Nvars):
        for jx in np.arange(ix+1,Nvars):
            matMI[ix,jx] = mutual_information((data[:,ix].reshape(-1, 1), data[:,jx].reshape(-1, 1)), k)
            matMI[jx,ix] = matMI[ix,jx]
    return matMI

def pearson_corr(matrix):
    """
    Function returns a matrix with Pearson correlation coefficients
    between pairs of variables (columns of the matrix) along with p-values.
    
    Parameters
    ----------
    matrix : array-like, shape (n, m)
             Array with each column being a different variable.
             
    Returns
    -------
    corr : array-like, shape (m, m)
           corr[i, j] contains the partial correlation of matrix[:, i] and matrix[:, j].
    prob : array-like, shape (m, m)
           prob[i, j] contains the p-value of a coefficient in corr[i, j].
    """
    (n, m) = matrix.shape

    DO = matrix - (np.sum(matrix, 0) / np.double(n))
    # note that mean row will be applyed row-wise to original matrices
    
    cov = np.einsum("nt,nm->tm", DO, DO)

    varO = np.sum(DO ** 2, 0)
    tmp = np.outer(varO, varO)
    
    corr = cov / np.sqrt(tmp)
    
    df = n-2
    
    diag = (np.diag(np.diag(corr)))
    corr -= diag
    t_squared = corr*corr*(df / ((1.0 - corr) * (1.0 + corr)))
    prob = special.betainc(0.5*df, 0.5, df / (df + t_squared))
    np.fill_diagonal(corr, 1)
    np.fill_diagonal(prob, 0)
    
    return corr, prob


def partial_corr(matrix):
    """
    Returns the sample linear partial correlation coefficients between pairs of variables in a matrix,
    controlling for the remaining variables in that matrix.
    
    Parameters
    ----------
    matrix : array-like, shape (n, m)
             Array with the different variables. Each column of the matrix is taken as a variable.
    
    Returns
    -------
    partial : array-like, shape (m, m)
              partial[i, j] contains the partial correlation of matrix[:, i] and matrix[:, j],
              controlling for the remaining variables in the matrix.
    prob : array-like, shape (m, m)
           prob[i, j] contains the p-value of a coefficient in corr[i, j].
    """

    n = matrix.shape[0]
    m = matrix.shape[1]
    ic = -np.linalg.pinv(np.cov(matrix, rowvar=0))
    diag1 = np.tile(np.sqrt(np.abs(np.diag(ic))),[m,1]).T
    diag2 = np.tile(np.sqrt(np.abs(np.diag(ic))),[m,1])
    partial = ((ic/diag1)/diag2)+2*np.eye(m)
    
    if n > m:
        df = n-m
    
        diag = (np.diag(np.diag(partial)))
        partial -= diag
        t_squared = partial*partial*(df / ((1.0 - partial) * (1.0 + partial)))
        prob = special.betainc(0.5*df, 0.5, df / (df + t_squared))
        np.fill_diagonal(partial, 1)
        np.fill_diagonal(prob, 0)
        return partial, prob
    else:
        return partial
    

def partial_partial_corr(matrix, alpha):
    """
    Returns the sample linear partial correlation coefficients between pairs of variables in a matrix,
    controlling for the remaining variables in that matrix, with additional option to partially partial out the variables.
    
    Parameters
    ----------
    matrix : array-like, shape (n, m)
             Array with the different variables. Each column of the matrix is taken as a variable.
    alpha :  parameter controlling the extent to which remaining variables are partialled out (0 - not at all, 1 - fully)
    
    Returns
    -------
    partial : array-like, shape (m, m)
              partial[i, j] contains the partial correlation of matrix[:, i] and matrix[:, j],
              controlling for the remaining variables in the matrix.
    prob : array-like, shape (m, m)
           prob[i, j] contains the p-value of a coefficient in corr[i, j].
    """
 
    n = matrix.shape[0]
    m = matrix.shape[1]
    partial = np.zeros((m, m), dtype=np.float)
    for i in range(m):
        partial[i, i] = 1
        for j in range(i+1, m):
            idx = np.ones(m, dtype=np.bool)
            idx[i] = False
            idx[j] = False
            beta_i = linalg.lstsq(matrix[:, idx], matrix[:, j])[0]
            beta_j = linalg.lstsq(matrix[:, idx], matrix[:, i])[0]

            res_j = matrix[:, j] - alpha*matrix[:, idx].dot(beta_i)
            res_i = matrix[:, i] - alpha*matrix[:, idx].dot(beta_j)

            corr = stats.pearsonr(res_i, res_j)[0]
            partial[i, j] = corr
            partial[j, i] = corr

    if n > m:
        df = n-m
    
        diag = (np.diag(np.diag(partial)))
        partial -= diag
        t_squared = partial*partial*(df / ((1.0 - partial) * (1.0 + partial)))
        prob = special.betainc(0.5*df, 0.5, df / (df + t_squared))
        np.fill_diagonal(partial, 1)
        np.fill_diagonal(prob, 0)
        return partial, prob
    else:
        return partial
     

def distance(matrix, metric="euclidean"):
    
    if metric=="euclidean":
        mat = spatial.distance.pdist(matrix.T)
    elif metric=="manhattan":
        mat = spatial.distance.pdist(matrix.T, metric='cityblock')
    
    return spatial.distance.squareform(mat)
    
    
def euclidean_distance(x,y):
    """Returns euclidean distance between two lists or numpy arrays"""
    x = np.array(x)
    y = np.array(y)
    return np.sqrt(sum((x-y)**2))


def manhattan_distance(x,y):
    """Returns manhattan distance between two lists or numpy arrays"""
    x = np.array(x)
    y = np.array(y)
    return sum(abs(x - y))


def calc_MI(X,Y,bins):
    """Returns Shannon's mutual information between two array-like objects"""
    c_XY = np.histogram2d(X,Y,bins)[0]
    c_X = np.histogram(X,bins)[0]
    c_Y = np.histogram(Y,bins)[0]

    H_X = shan_entropy(c_X)
    H_Y = shan_entropy(c_Y)
    H_XY = shan_entropy(c_XY)

    MI = H_X + H_Y - H_XY
    MI = np.log(2**MI)
    return MI

def shan_entropy(c):
    """Retuurns Shannon's information (entropy) for arrray-like object"""
    c_normalized = c / float(np.sum(c))
    c_normalized = c_normalized[np.nonzero(c_normalized)]
    H = -sum(c_normalized* np.log2(c_normalized))
    return H

def detrend(matrix):
    # remove global mean from the data:
    detrended = matrix - np.matlib.repmat(np.mean(matrix,axis=0), matrix.shape[0],1)
    return detrended


def calc_conn(data, func, bins=10):
    
    n_vars = data.shape[1]
    
    conn = np.zeros((n_vars, n_vars))
    conn_p = np.zeros((n_vars, n_vars))
    
    if str(func)[10:14] == 'pear':
        print('Pearson')
        for i in range(n_vars):
            for j in range(n_vars):
                    (conn[i, j], 
                     conn_p[i, j]) = func(data.iloc[:,i], data.iloc[:,j])

    elif str(func)[10:14] == 'calc':
        print('MI')    
        for i in range(n_vars):
            for j in range(i+1, n_vars):
                conn[i, j] = func(data.iloc[:,i], data.iloc[:,j], bins)
                conn[j, i] = conn[i, j]
                
    return conn


def visualize(what, title, figsize=(5,5)):
    f = plt.figure(figsize=figsize)
    ax = plt.gca()
    im = plt.imshow(what, clim=[-0.5,0.5],cmap='coolwarm')
    im = plt.imshow(what, clim=[-1,1],cmap='coolwarm')
    plt.title(title, fontsize=20)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    plt.show()
    
 
def RBIG_MI(data):
    """ Estimate Mutual Information with RBIG """
    n_vars = data.shape[1]
    conn = np.zeros((n_vars, n_vars))
    n_layers = 10
    rotation_type = 'PCA'
    random_state = 0
    zero_tolerance = 60
    tolerance = None
    print('RBIGMI')
    # Initialize RBIG class
    rbig_model = RBIGMI(n_layers=n_layers, rotation_type=rotation_type, random_state=random_state, 
        zero_tolerance=zero_tolerance, tolerance=tolerance)

    for i in range(n_vars):
        for j in range(i+1, n_vars):
            rbig_model.fit(data[:,i].reshape(-1, 1), data[:,j].reshape(-1, 1))
            mi_rbig = rbig_model.mutual_information() * np.log(2)
            conn[i, j] = mi_rbig
            conn[j, i] = conn[i, j]
    return conn


def RBIG_TC(data):
    """ Estimate Total Correlation with RBIG """
    n_vars = data.shape[1]
    conn = np.zeros((n_vars, n_vars))
    n_layers = 10
    rotation_type = 'PCA'
    random_state = 0
    zero_tolerance = 60
    tolerance = None
    print('RBIGTC')
    # Initialize RBIG class
    rbig_model = RBIG(n_layers=n_layers, rotation_type=rotation_type, random_state=random_state, 
        zero_tolerance=zero_tolerance, tolerance=tolerance)

    for i in range(n_vars):
        for j in range(i+1, n_vars):
            rbig_model.fit(np.hstack((data[:,i].reshape(-1, 1), data[:,j].reshape(-1, 1))))
            tc_rbig = rbig_model.mutual_information * np.log(2)
            conn[i, j] = tc_rbig
            conn[j, i] = conn[i, j]
    return conn

    
def Corex_TC(data):
    """ Estimate Total Correlation with Corex """
    n_vars = data.shape[1]
    conn = np.zeros((n_vars, n_vars))
    layer = ce.Corex(n_hidden=1)
    print('CorExTC')
    for i in range(n_vars):
        for j in range(i+1, n_vars):
            layer = layer.fit(np.hstack((data[:,i].reshape(-1, 1), data[:,j].reshape(-1, 1))))
            mi_corex= layer.tc * 10
            conn[i, j] = mi_corex
            conn[j, i] = conn[i, j]          
    return conn

def graph_info(infos, labels, colors=None):
    '''help docs: https://github.com/gabeschamberg/whats_that_smell/blob/master/notebooks/olfactory_bulb_te.ipynb'''
    G = nx.DiGraph()
    G.graph['rankdir'] = 'LR'
    G.graph['dpi'] = 120
    # Create nodes
    for label in labels:
        if colors is not None:
            G.add_node(label, style='filled', fillcolor=colors[label])
        else:
            G.add_node(label, style='filled')
    
    for info in infos:
        G.add_edge(labels[info[0]], labels[info[1]], label=str(round(info[2], 4)), penwidth=str(1+10*info[2]**2))
    draw(G)
    return G


def create_corr_network_1(G):
    '''help docs:https://towardsdatascience.com/visualising-stocks-correlations-with-networkx-88f2ee25362e'''
    #crates a list for edges and for the weights
    edges, weights = zip(*nx.get_edge_attributes(G,'weight').items())

    positions=nx.circular_layout(G)
    
    #Figure size
    plt.figure(figsize=(15,15))

    nx.draw_networkx_nodes(G,positions,node_color='#DA70D6',
                           node_size=500,alpha=0.8)
    
    nx.draw_networkx_labels(G, positions, font_size=8, 
                            font_family='sans-serif')
        
    nx.draw_networkx_edges(G, positions, edgelist=edges, style='solid')
    
    plt.axis('off')
    #plt.savefig("part1.png", format="PNG")
    plt.show() 


def make_graph(weights, node_weights, column_label, max_edges=100):
    '''The function adative from https://github.com/gregversteeg/bio_corex/blob/master/vis_corex.py'''
    all_edges = np.hstack(list(map(np.ravel, weights)))
    max_edges = min(max_edges, len(all_edges))
    w_thresh = np.sort(all_edges)[-max_edges]
    print('weight threshold is %f for graph with max of %f edges ' % (w_thresh, max_edges))
    g = nx.DiGraph()
    max_node_weight = max([max(w) for w in node_weights])
    for layer, weight in enumerate(weights):
        m, n = weight.shape
        for j in range(m):
            g.add_node((layer + 1, j))
            g._node[(layer + 1, j)]['weight'] = 0.3 * node_weights[layer][j] / max_node_weight
            for i in range(n):
                if weight[j, i] > w_thresh:
                    if weight[j, i] > w_thresh / 2:
                        g.add_weighted_edges_from([( (layer, i), (layer + 1, j), 10 * weight[j, i])])
                    else:
                        g.add_weighted_edges_from([( (layer, i), (layer + 1, j), 0)])

    # Label layer 0
    for i, lab in enumerate(column_label):
        g.add_node((0, i))
        g._node[(0, i)]['label'] = lab
        g._node[(0, i)]['name'] = lab  # JSON uses this field
        g._node[(0, i)]['weight'] = 1
    return g

def make_graph_ql(weights, node_weights, column_label, max_edges=100):
    '''The function adative from https://github.com/gregversteeg/bio_corex/blob/master/vis_corex.py'''
    all_edges = np.hstack(list(map(np.ravel, weights)))
    max_edges = min(max_edges, len(all_edges))
    w_thresh = np.sort(all_edges)[-max_edges]
    print('weight threshold is %f for graph with max of %f edges ' % (w_thresh, max_edges))
    g = nx.DiGraph()
    max_node_weight = max([w for w in node_weights])
    for weight in enumerate(weights):
        m, n = weight.shape
        for j in range(m):
            g.add_node((layer + 1, j))
            g._node[(layer + 1, j)]['weight'] = 0.3 * node_weights[layer][j] / max_node_weight
            for i in range(n):
                if weight[j, i] > w_thresh:
                    if weight[j, i] > w_thresh / 2:
                        g.add_weighted_edges_from([( (layer, i), (layer + 1, j), 10 * weight[j, i])])
                    else:
                        g.add_weighted_edges_from([( (layer, i), (layer + 1, j), 0)])

    # Label layer 0
    for i, lab in enumerate(column_label):
        g.add_node((0, i))
        g._node[(0, i)]['label'] = lab
        g._node[(0, i)]['name'] = lab  # JSON uses this field
        g._node[(0, i)]['weight'] = 1
    return g

def trim(g, max_parents=False, max_children=False):
    '''The function adative from https://github.com/gregversteeg/bio_corex/blob/master/vis_corex.py'''
    if float(nx.__version__) < 2:
        edgedict = g.edge
    else:
        edgedict = g.adj
    for node in g:
        if max_parents:
            parents = list(g.successors(node))
            weights = [edgedict[node][parent]['weight'] for parent in parents]
            for weak_parent in np.argsort(weights)[:-max_parents]:
                g.remove_edge(node, parents[weak_parent])
        if max_children:
            children = g.predecessors(node)
            weights = [edgedict[child][node]['weight'] for child in children]
            for weak_child in np.argsort(weights)[:-max_children]:
                g.remove_edge(children[weak_child], node)
    return g


def safe_open(filename, mode):
    '''The function adative from https://github.com/gregversteeg/bio_corex/blob/master/vis_corex.py'''
    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    return open(filename, mode)


# Visualization utilities

def neato(fname, position=None, directed=False):
    '''The function adative from https://github.com/gregversteeg/bio_corex/blob/master/vis_corex.py'''
    if directed:
        os.system(
            "sfdp " + fname + ".dot -Tpdf -Earrowhead=none -Nfontsize=12  -GK=2 -Gmaxiter=1000 -Goverlap=False -Gpack=True -Gpackmode=clust -Gsep=0.01 -Gsplines=False -o " + fname + "_sfdp.pdf")
        os.system(
            "sfdp " + fname + ".dot -Tpdf -Earrowhead=none -Nfontsize=12  -GK=2 -Gmaxiter=1000 -Goverlap=False -Gpack=True -Gpackmode=clust -Gsep=0.01 -Gsplines=True -o " + fname + "_sfdp_w_splines.pdf")
        return True
    if position is None:
        os.system("neato " + fname + ".dot -Tpdf -o " + fname + ".pdf")
        os.system("fdp " + fname + ".dot -Tpdf -o " + fname + "fdp.pdf")
    else:
        os.system("neato " + fname + ".dot -Tpdf -n -o " + fname + ".pdf")
    return True


def extract_color(label):
    '''The function adative from https://github.com/gregversteeg/bio_corex/blob/master/vis_corex.py'''
    colors = 'indigo,gold,hotpink,firebrick,indianred,yellow,mistyrose,darkolivegreen,darkseagreen,pink,tomato,lightcoral,orangered,navajowhite,palegreen,darkslategrey,greenyellow,burlywood,seashell,mediumspringgreen,papayawhip,blanchedalmond,chartreuse,dimgray,black,peachpuff,springgreen,aquamarine,white,orange,lightsalmon,darkslategray,brown,ivory,dodgerblue,peru,lawngreen,chocolate,crimson,forestgreen,slateblue,lightseagreen,cyan,mintcream,antiquewhite,mediumorchid,skyblue,gray,darkturquoise,goldenrod,darkgreen,floralwhite,darkviolet,moccasin,saddlebrown,grey,darkslateblue,lightskyblue,lightpink,mediumvioletred,slategrey,red,deeppink,limegreen,palegoldenrod,plum,turquoise,lightgrey,lightgoldenrodyellow,darkgoldenrod,lavender,maroon,yellowgreen,sandybrown,thistle,violet,navy,magenta,dimgrey,tan,rosybrown,blue,lightblue,ghostwhite,honeydew,cornflowerblue,linen,powderblue,seagreen,darkkhaki,snow,sienna,mediumblue,royalblue,lightcyan,green,mediumpurple,midnightblue,cornsilk,paleturquoise,bisque,slategray,khaki,wheat,darkorchid,deepskyblue,salmon,steelblue,palevioletred,lightslategray,aliceblue,lightslategrey,orchid,gainsboro,mediumseagreen,lightgray,mediumturquoise,lemonchiffon,cadetblue,lightyellow,lavenderblush,coral,purple,whitesmoke,mediumslateblue,darkorange,mediumaquamarine,darksalmon,beige,blueviolet,azure,lightsteelblue,oldlace'.split(',')

    parts = label.split('_')
    for part in parts:
        if part in colors:
            parts.remove(part)
            return '_'.join(parts), part
    return label, 'black'


def edge2pdf(g, filename, threshold=0, position=None, labels=None, connected=True, directed=False, makepdf=True):
    #This function will takes list of edges and a filename
    #and write a file in .dot format. Readable, eg. by omnigraffle
    # OR use "neato file.dot -Tpng -n -o file.png"
    # The -n option says whether to use included node positions or to generate new ones
    # for a grid, positions = [(i%28,i/28) for i in range(784)]
    '''The function adative from https://github.com/gregversteeg/bio_corex/blob/master/vis_corex.py'''
    def cnn(node):
        #change node names for dot format
        if type(node) is tuple or type(node) is list:
            return u'n' + u'_'.join(map(str, node))
        else:
            return str(node)

    if connected:
        touching = list(set(sum([[a, b] for a, b in g.edges()], [])))
        g = nx.subgraph(g, touching)
        print('non-isolated nodes,edges', len(list(g.nodes())), len(list(g.edges())))
    f = safe_open(filename + '.dot', 'wb+')
    if directed:
        f.write("strict digraph {\n".encode('utf-8'))
    else:
        f.write("strict graph {\n".encode('utf-8'))
    #f.write("\tgraph [overlap=scale];\n".encode('utf-8'))
    f.write("\tnode [shape=point];\n".encode('utf-8'))
    for a, b, d in g.edges(data=True):
        if 'weight' in d:
            if directed:
                f.write(("\t" + cnn(a) + ' -> ' + cnn(b) + ' [penwidth=%.2f' % float(
                    np.clip(d['weight'], 0, 9)) + '];\n').encode('utf-8'))
            else:
                if d['weight'] > threshold:
                    f.write(("\t" + cnn(a) + ' -- ' + cnn(b) + ' [penwidth=' + str(3 * d['weight']) + '];\n').encode(
                        'utf-8'))
        else:
            if directed:
                f.write(("\t" + cnn(a) + ' -> ' + cnn(b) + ';\n').encode('utf-8'))
            else:
                f.write(("\t" + cnn(a) + ' -- ' + cnn(b) + ';\n').encode('utf-8'))
    for n in g.nodes():
        if labels is not None:
            if type(labels) == dict or type(labels) == list:
                thislabel = labels[n].replace(u'"', u'\\"')
                lstring = u'label="' + thislabel + u'",shape=none'
            elif type(labels) == str:
                if 'label' in g._node[n]:
                    thislabel = g._node[n][labels].replace(u'"', u'\\"')
                    # combine dupes
                    #llist = thislabel.split(',')
                    #thislabel = ','.join([l for l in set(llist)])
                    thislabel, thiscolor = extract_color(thislabel)
                    lstring = u'label="%s",shape=none,fontcolor="%s"' % (thislabel, thiscolor)
                else:
                    weight = g._node[n].get('weight', 0.1)
                    if n[0] == 1:
                        lstring = u'shape=circle,margin="0,0",style=filled,fillcolor=black,fontcolor=white,height=%0.2f,label="%d"' % (
                            2 * weight, n[1])
                    else:
                        lstring = u'shape=point,height=%0.2f' % weight
            else:
                lstring = 'label="' + str(n) + '",shape=none'
            #lstring = unicode(lstring)
        else:
            lstring = False
        if position is not None:
            if position == 'grid':
                position = [(i % 28, 28 - i / 28) for i in range(784)]
            posstring = 'pos="' + str(position[n][0]) + ',' + str(position[n][1]) + '"'
        else:
            posstring = False
        finalstring = u' [' + u','.join([ts for ts in [posstring, lstring] if ts]) + u']\n'
        #finalstring = u' ['+lstring+u']\n'
        f.write((u'\t' + cnn(n) + finalstring).encode('utf-8'))
    f.write("}".encode('utf-8'))
    f.close()
    if makepdf:
        neato(filename, position=position, directed=directed)
    return True


def create_corr_network(G, corr_direction,threshold):
    '''Graph visualization functional connectivity'''
    H = G.copy()
    for sk1, sk2, weight in G.edges(data=True):
        if corr_direction == "positive":
            if weight["weight"] <threshold:
                H.remove_edge(sk1, sk2)
        else:
            if weight["weight"] >=threshold:
                H.remove_edge(sk1, sk2)
                
    #crates a list for edges and for the weights
    edges,weights = zip(*nx.get_edge_attributes(H,'weight').items())
    
    ### increases the value of weights, so that they are more visible in the graph
    weights = tuple([(1+abs(x))**3.5 for x in weights])
        
    #####calculates the degree of each node
    #d = nx.degree(H) #only works on network 1.x
    d =dict(H.degree)
    #####creates list of nodes and a list their degrees that will be used later for their sizes
    
    nodelist, node_sizes = zip(*d.items())
    
    #positions
    positions=nx.circular_layout(H)
    
    #g = make_graph(weights, node_sizes, nodelist, max_edges=max_edges)
    # Display pruned version
    #h = g.copy()  # trim(g.copy(), max_parents=max_parents, max_children=max_children)
    #edge2pdf(h, prefix + '/graphs/graph_prune_' + str(max_edges), labels='label', directed=True, makepdf=True)
    #Display tree version
    #tree = g.copy()
    #tree = trim(tree, max_parents=1, max_children=False)
    #edge2pdf(tree, prefix + '/graphs/tree', labels='label', directed=True, makepdf=True)

    #Figure size
    plt.figure(figsize=(10,10))

    #draws nodes (nodes color: #DA70D6)
    nx.draw_networkx_nodes(H,positions, node_color='#000000', nodelist=nodelist,       
                           #####the node size will be now based on its degree
                           node_size=tuple([x**4.3 for x in node_sizes]),alpha=0)    
    #Styling for labels
    nx.draw_networkx_labels(H, positions, font_size=12, font_family='bold', font_color='b')
    ###edge colors based on weight direction
    if corr_direction == "positive":  
        edge_colour = plt.cm.inferno#plt.cm.spring
    else:
        edge_colour = plt.cm.GnBu
        
    #draws the edges
    nx.draw_networkx_edges(H, positions, edgelist=edges, style='solid',
                          ###adds width=weights and edge_color = weights 
                          ###so that edges are based on the weight parameter 
                          ###edge_cmap is for the color scale based on the weight
                          ### edge_vmin and edge_vmax assign the min and max weights for the width
                          width=weights, edge_color = weights, edge_cmap = edge_colour,
                          edge_vmin = min(weights), edge_vmax=max(weights))

    # displays the graph without axis
    plt.axis('off')
    #saves image
    #plt.savefig("part5" + corr_direction + ".png", format="PNG")
    #plt.show() 
    return G

##############################################################################################################
# PDF  
##############################################################################################################
TR = 1.89
df1 = pd.read_csv(r"fmri_timeseries_noHeader.csv", delimiter=',', header=None)
data1  = df1.values
Nvars  = data1.shape[1]
data2=data1[:, 1:]
data2=np.array(data2)
print(data2.shape)
data3 = percent_change(data2)
Ta = TimeSeries(data3, sampling_interval=TR)
data = Ta.data
data_pc = pd.DataFrame(data)
labels= pd.read_csv('fmri_timeseries.csv', nrows=0).columns.tolist()

# Pearson Correlation
pearson = calc_conn(data_pc, stats.pearsonr)
cor_matrix = np.asmatrix(pearson)
df = pd.DataFrame(cor_matrix)
df.columns = labels
df.index = labels
df.to_csv('fmri_timeseries_corr.csv', index=True)

G = nx.from_numpy_matrix(cor_matrix)
G = nx.relabel_nodes(G,lambda x: labels[x-1])
G.edges(data=True)
#create_corr_network_1(G)
G = create_corr_network(G, corr_direction='positive', threshold=0.14)


tree = G.copy()
tree = trim(G, max_parents=1, max_children=False)
prefix = 'cc_graph'  
edge2pdf(tree, prefix + '/graphs/tree', labels='label', directed=True, makepdf=True)


'''
corr_direction="positive"
H= G.copy()
for sk1, sk2, weight in G.edges(data=True):
    if corr_direction == "positive":
        if weight["weight"] <0.001:
            H.remove_edge(sk1, sk2)
    else:
        if weight["weight"] >=0:
            H.remove_edge(sk1, sk2)
            
#crates a list for edges and for the weights
edges,weights = zip(*nx.get_edge_attributes(H,'weight').items())

### increases the value of weights, so that they are more visible in the graph
weights = tuple([(1+abs(x))**2 for x in weights])
d =dict(H.degree)
nodelist, node_sizes = zip(*d.items())
node_weights=tuple([x for x in node_sizes]) 
'''
# Display pruned version
h = G.copy()  # trim(g.copy(), max_parents=max_parents, max_children=max_children)
edge2pdf(h, 'test/graphs/tree'  + '/graphs/graph_prune_' + str(35), labels=str(labels), directed=True, makepdf=True)

# Display tree version
tree = G.copy()
tree = trim(tree, max_parents=1, max_children=False)
edge2pdf(tree, 'test/graphs/tree' + '/graphs/tree', labels=str(labels), directed=True, makepdf=True)



# Shannon Mutual Information
#mi = calc_conn(data_pc, calc_MI, bins=int(np.sqrt(data_pc.shape[0])))
#mi_matrix = np.asmatrix(mi)
#df = pd.DataFrame(mi_matrix)
#df.columns = labels
#df.index = labels
#df.to_csv('fmri_timeseries_mi.csv', index=True)
#G = nx.from_numpy_matrix(mi_matrix)
#G = nx.relabel_nodes(G,lambda x: labels[x-1])
#G.edges(data=True)
#create_corr_network(G, corr_direction='positive',threshold=0.32)

'''

# Kraskov Mutual Information
kmi16 = KMI(data, 25)
kmi_matrix = np.asmatrix(kmi16)
G = nx.from_numpy_matrix(kmi_matrix)
G = nx.relabel_nodes(G,lambda x: labels[x-1])
G.edges(data=True)
create_corr_network(trim(G), corr_direction='positive')

# RBIG Mutual Information
rbigmi=RBIG_MI(data)
rbigmi_matrix = np.asmatrix(rbigmi)
G = nx.from_numpy_matrix(rbigmi_matrix)
G = nx.relabel_nodes(G,lambda x: labels[x-1])
G.edges(data=True)
create_corr_network(trim(G), corr_direction='positive')

#RBIG Total Correlation=RBIG Mutual Information
rbigtc=RBIG_TC(data)
rbigtc_matrix = np.asmatrix(rbigtc)
G = nx.from_numpy_matrix(rbigtc_matrix)
G = nx.relabel_nodes(G,lambda x: labels[x-1])
G.edges(data=True)
create_corr_network(trim(G), corr_direction='positive')

#Corex Total Information = CorEx Mutual Information
corextc=Corex_TC(data)
corextc_matrix = np.asmatrix(corextc)
G = nx.from_numpy_matrix(corextc_matrix)
G = nx.relabel_nodes(G,lambda x: labels[x-1])
G.edges(data=True)
create_corr_network(trim(G), corr_direction='positive')
'''
