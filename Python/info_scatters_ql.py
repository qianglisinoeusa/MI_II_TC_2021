"""
Info vs correlation scatter plots for intro to MI slides

"""
import numpy as np
import numpy.random
import matplotlib.pyplot as plt
import sys
sys.path.append('pyentropy-master/')  
from pyentropy.utils import quantise
from pyentropy.systems import DiscreteSystem
from pyitlib import discrete_random_variable as drv
import infotheory
pi = np.pi

def plot_calc(x,y,ax=None, xlim=[-4, 4], ylim=[-4,4]):
    if ax is None:
        f = plt.figure()
        ax = f.add_subplot(111)
    # correlation
    cor = np.corrcoef(x,y)[0,1]

    # information
    #m = 8 
    #qx = quantise(x, m, uniform='sampling')[0]
    #qy = quantise(y, m, uniform='sampling')[0]
    #s = DiscreteSystem(qx, (1,m), qy, (1,m))
    #s.calculate_entropies(method='plugin', calc=['HX','HXY'])
    #I = s.I()
    #mi = drv.information_mutual(x,y)
    #I = mi
    it = infotheory.InfoTools(2, 3)
    it.set_equal_interval_binning([10]*2, [np.min(x), np.min(y)], [np.max(x), np.max(y)])
    it.add_data(np.vstack([x, y]).T)
    I = it.mutual_info([0, 1]) / np.log2(8)

    Nplot = 600
    plotidx = np.random.permutation(len(x))[:Nplot]
    ax.scatter(x[plotidx],y[plotidx], s=5, edgecolors='none', c="blue")
    ax.set_title("Corr=%.1f  MI=%.2f" % (cor, I), fontsize=10, fontweight='bold')
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    for sp in ['left','right','bottom','top']:
        ax.spines[sp].set_color('none')
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if np.abs(cor)<1e-3:
        cor = 0.0
    if np.isnan(cor):
        cor = '-'
    else:
        cor = "%.2f"%cor
    #ax.text(0.2, 1, cor, horizontalalignment='center', verticalalignment='bottom', 
    #        transform = ax.transAxes, color='#663300', fontsize=12, fontweight='bold')
    #ax.text(0.8, 1, "%.2f"%I, horizontalalignment='center', verticalalignment='bottom', 
    #        transform = ax.transAxes, color='#b20000', fontsize=12, fontweight='bold')
    return ax

def gen_mv_normal(n, cor, ax=None):
    sd = np.array([ [1, cor], [cor, 1] ])
    mean = np.array([0, 0])
    xy = np.random.multivariate_normal(mean, sd, size=n)
    x = xy[:,0]
    y = xy[:,1]
    plot_calc(x,y,ax)
    
def rotate(xy, t):
    R = np.array([ [np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)] ]).T
    return np.dot(xy, R)

def gen_rot_normal(n, t, ax=None):
    sd = np.array([ [1,1],[1,1] ])
    mean = np.array([0,0])
    xy = np.random.multivariate_normal(mean, sd, size=n)
    xyr = rotate(xy, t)
    x = xyr[:,0]
    y = xyr[:,1]
    # round to 0 for flat one
    y[np.abs(y)<1e-13] = 0
    plot_calc(x,y,ax)

def gen_others(n, idx, ax=None):
    x = np.linspace(-1,1,n)

    # idx = index of plot to generate
    if idx==1:
        # curvy
        r = (np.random.random(n)*3) - 1
        y = 4.0* (x**2 - 0.5)**2 + (r/3)
        plot_calc(x,y,ax,xlim=[-1,1],ylim=[-1/3.0, 1+(1/3.0)])
    if idx==2:
        # rotated uniform
        y = np.random.random(n)*2 - 1
        xy = rotate(np.c_[x,y], -np.pi/8.0)
        lim = np.sqrt(2 + np.sqrt(2)) / np.sqrt(2)
        plot_calc(xy[:,0],xy[:,1],ax,xlim=[-lim,lim], ylim=[-lim,lim])
    if idx==3:
        # circle
        r = np.random.normal(0, 1/8.0, n)
        y = np.cos(x*pi) + r
        r = np.random.normal(0, 1/8.0, n)
        x = np.sin(x*pi) + r
        plot_calc(x,y,ax,xlim=[-1.5,1.5],ylim=[-1.5,1.5])
    if idx==4:
        # smile
        r = (np.random.random(n)*2) - 1
        y = 2*(x**2) + r
        plot_calc(x,y,ax,xlim=[-1,1],ylim=[-1,3])
    if idx==5:
        # rotated uniform
        y = np.random.random(n)*2 - 1
        xy = rotate(np.c_[x,y], -np.pi/4.0)
        lim = np.sqrt(2)
        plot_calc(xy[:,0],xy[:,1],ax,xlim=[-lim,lim], ylim=[-lim,lim])
    
    
def build_figure():
    n = 100000
    f = plt.figure(figsize=(12,4))

    # MV NORMAL
    mv_ax = []
    for i in range(4):
        mv_ax.append( f.add_subplot(2,4,i+1) )
    mv_cors = [0.8, 0.4, -0.4, -0.8]
    for i in range(4):
        gen_mv_normal(n, mv_cors[i], mv_ax[i])

    # ROT NORMAL
    #rot_ax = []
    #for i in range(7):
        #rot_ax.append( f.add_subplot(3,7,i+8) )
    #rot_vals = [0, pi/12, pi/6, pi/4, pi/2 - pi/6, pi/2 - pi/12, pi/2]
    #for i in range(7):
        #gen_rot_normal(n, rot_vals[i], rot_ax[i])

    # OTHERS
    other_ax = []
    for i in range(4):
        other_ax.append( f.add_subplot(2,4,i+5) )
    other_vals = range(1,5)
    for i in range(4):
        gen_others(n, other_vals[i], other_ax[i])

    f.subplots_adjust(left=0.01, top=0.95, right=0.99, bottom=0.01)

    return f
    
if __name__ == '__main__':
    f=build_figure()
    plt.show()