
from pnm_cflow import Pnm_Cflow


def Make_Mheur(model, args, MIPHEUR):
    """
    Create MIP model for heuristics based on the Pnm_Cflow template
    (more complex alternatives are also possible).
    """
    if MIPHEUR is not None and MIPHEUR['enabled']:
        HEUR = dict()
        HEUR['enabled'] = False
        HEUR['Theur'] = -1  
        HEUR['T_lim'] =  0 
        o_heur = Pnm_Cflow(args, RESLIM=MIPHEUR['T_lim'], OPTCR=0.999999, MIPSTART=False, CONNTYP='flow') 
        m_heur = o_heur.model        
        m_heur.Params.TimeLimit = MIPHEUR['T_lim']   # max time per heuristic run
        m_heur._heurst = False
        m_heur._Theur  = -1
        m_heur.update()
        model.update()
        #---------------------------------------------------------------------------------
        m_heur.Params.PreCrush = model.Params.PreCrush
        m_heur.Params.LazyConstraints = model.Params.LazyConstraints
        #---------------------------------------------------------------------------------
        x_lbnd = dict()
        for i,k in m_heur._x:
            x_lbnd[i,k] = m_heur._x[i,k].lb
        x_ubnd = dict()
        for i,k in m_heur._x:
            x_ubnd[i,k] = m_heur._x[i,k].ub            
        y_lbnd = dict()
        for i,j in m_heur._y:
            y_lbnd[i,j] = m_heur._y[i,j].lb
        y_ubnd = dict()
        for i,j in m_heur._y:
            y_ubnd[i,j] = m_heur._y[i,j].ub
        # To reset the original bounds after each run of the heuristic 
        m_heur._x_lbnd = x_lbnd
        m_heur._x_ubnd = x_ubnd        
        m_heur._y_lbnd = y_lbnd
        m_heur._y_ubnd = y_ubnd    
        m_heur._thrsh  = 0.9    # INITIAL rounding threshold for LP relaxation    
        m_heur._thrshacc = 1    # number of significant digits in _thrsh
        m_heur._roudg  = 0.71   # INITIAL percentage of fixed nodes from LP relaxation        
        m_heur.Params.LogToConsole = 0     
        # Copy the objective from the original model
        for k in model._IMBLOD:
            m_heur._IMBLOD[k].Obj = model._IMBLOD[k].Obj        
        for l in model._pLS:
            m_heur._pLS[l].Obj = model._pLS[l].Obj        
        for g in model._pGS:
            m_heur._pGS[g].Obj = model._pGS[g].Obj        
        for i,j in model._y:
            m_heur._y[i,j].Obj = model._y[i,j].Obj
        # Set model attributes
        m_heur.update()
        model._mheur   = m_heur
        model._Theur   = MIPHEUR['Theur']    # remaining time to be 
        model._heurst  = MIPHEUR['enabled']  # used on MIP heuristics
        model._MipHeur = MIPHEUR
        model._infeascount = 0
        model._infeascycle = 0
    else:
        model._heurst = False
        model._MipHeur = None
    return model


