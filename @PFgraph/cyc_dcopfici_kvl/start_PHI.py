"""
The corresponding file that calls the optimize() method for the model
is bnc_dcopf.py. There the MIPFOCUS parameter for the solution strategy 
is set. 
"""

import os 
import sys 
sys.path.append(os.path.dirname(os.path.dirname( __file__ )))
from grb_read_graph import read_file, unidirect 
from gurobipy  import setParam
from pnm_cflow import Pnm_Cflow
from make_mheur import Make_Mheur

# Input parameters for optimization
MIPGap = 1e-2
START = False
HEUR = dict()
HEUR['enabled'] = False
ici_pth = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'ici'))
       
# Initial solution (mostly to debug Gurobi by fixing feasible x(i,k))
x_ini, _ = read_file(os.path.join(ici_pth, 'x_ini.csv'))
#z_ini, _ = read_file(os.path.join(ici_pth, 'z_ini.csv'))
y_ini, _ = read_file(os.path.join(ici_pth, 'y_ini.csv'))

# Graph parameters
bdcpf, nodes = read_file(os.path.join(ici_pth, 'graph.csv'))
flow0, _ = read_file(os.path.join(ici_pth, 'flow0.csv'))
brlim, _ = read_file(os.path.join(ici_pth, 'brlim.csv'))

# Sets
buses = read_file(os.path.join(ici_pth, 'nodes.csv'))
islid = read_file(os.path.join(ici_pth, 'islid.csv'))
islsw = read_file(os.path.join(ici_pth, 'islsw.csv'))
loads = read_file(os.path.join(ici_pth, 'loads.csv'))
genss = read_file(os.path.join(ici_pth, 'genss.csv'))
mu    = read_file(os.path.join(ici_pth, 'mu.csv'))
nu    = read_file(os.path.join(ici_pth, 'nu.csv'))
mu    = float(mu.pop())
nu    = float(nu.pop())

# Study case name to load minimal cycle basis
f = open(os.path.join(ici_pth, 'caseid.txt'), 'r')
caseid = f.read(); f.close()

# Nodal parameters  
pL0   = read_file(os.path.join(ici_pth,   'pL0.csv'))
pG0   = read_file(os.path.join(ici_pth,   'pG0.csv'))
trm2k = read_file(os.path.join(ici_pth, 'cohgr.csv'))
pLmax = read_file(os.path.join(ici_pth, 'pLmax.csv'))
pLmin = read_file(os.path.join(ici_pth, 'pLmin.csv'))
pGmax = read_file(os.path.join(ici_pth, 'pGmax.csv'))
pGmin = read_file(os.path.join(ici_pth, 'pGmin.csv'))
costL = read_file(os.path.join(ici_pth, 'costL.csv'))

# Unidirect input graphs
bdcpf = unidirect(bdcpf)
flow0 = unidirect(flow0)
brlim = unidirect(brlim)

# Create optimization model
method = 'pnm_cflow'    #'bnc_xarcs' 'bnc_nodes'
setParam('DisplayInterval', 10)
assert set(nodes)==buses, "Set of nodes implied from admittance graph should be equal to set of nodes read from nodes.csv"
args = {'nodes':nodes,
        'islid':islid,
        'islsw':islsw,
        'loads':loads,
        'genss':genss,
        'bdcpf':bdcpf,
        'flow0':flow0,
        'brlim':brlim,
        'x_ini':x_ini,
        'y_ini':y_ini,        
        'z_ini':None,        
        'pL0'  :pL0,
        'pG0'  :pG0,
        'costL':costL,
        'trm2k':trm2k,
        'pLmax':pLmax,
        'pLmin':pLmin,
        'pGmax':pGmax,
        'pGmin':pGmin,
        'icase':caseid,
        'mu'   :mu,
        'nu'   :nu}
if len(buses)<500:
    TimeLimit = 480
else:
    TimeLimit = 720
HEUR['Theur'] = TimeLimit*0.06
HEUR['T_lim'] = TimeLimit*0.02
if method=='pnm_cflow':
    m = Pnm_Cflow(args, RESLIM=TimeLimit, OPTCR=MIPGap, MIPSTART=START, CONNTYP='flow')   #'tree' 'flow'
    m.model.Params.PreCrush = 1
    # Create the auxiliary MIP heuristics model here to avoid recursive module imports in pnm_cflow.py
    if HEUR['enabled']:
        m.model = Make_Mheur(m.model, args, HEUR)
    else:
        m.model._MipHeur = HEUR
    ##m.model.Params.CutPasses = 10000000
    #m.model.Params.ImproveStartTime = 0.02*TimeLimit
    m.model.Params.LazyConstraints = 1 
    m.solve()
else:
    raise RuntimeError('The value of method argument can be either "bnc_arbor" or "bnc_nodes"')







