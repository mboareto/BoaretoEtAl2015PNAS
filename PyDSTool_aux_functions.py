#-------------------------------------------------------------------------------------------------------------------#
#  Auxilary functions to use with PyTDSTool.
#  Writen by Marcelo Boareto: marceloboareto@gmail.com
#  Last update: 11/2014
#-------------------------------------------------------------------------------------------------------------------#

from PyDSTool import *
from PyDSTool.Toolbox import phaseplane as pp
import random as rand
import copy


#-------------------------------------------------------------------------------------------------------------------#
def functions():
### User defined functions. 
    return {'Hp' : (['X','X0','nX'],'((X/X0)**nX)/(1. + (X/X0)**nX)'     )  # Positive Hill function
           ,'Hn' : (['X','X0','nX'],       '(1.0)/(1. + (X/X0)**nX)'     )  # Negative Hill function
           ,'HS' : (['X','X0','nX','lamb'],'(1-lamb)*Hn(X,X0,nX) + lamb' )  # Shifted Hill function
           }
#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#
def eliminate_redundants(fp, eps=10):
### Eliminate redundants fixed points. eps=10 
    for i in range(len(fp)):
        for k, v in fp[i].items():
            v = round(v,eps) 
            fp[i][k] = v
    seen = set()
    new_l = []
    for d in fp:
        t = tuple(d.items())
        if t not in seen:
            seen.add(t)
            new_l.append(d)
    return new_l
#-------------------------------------------------------------------------------------------------------------------#
 
#-------------------------------------------------------------------------------------------------------------------#
def PyCont_args(name, freepar, maxnumpoints, maxstep=1e+1, minstep=1e-1, 
                step=1e-0, LocBifPoints='all', saveeigen=False, Type='EP-C'):
### Define PyCont arguments
    PCargs = PyDSTool.args(name=name, type=Type)      # 'EP-C' stands for Equilibrium Point Curve.
    PCargs.freepars     = [freepar]                   # control parameter 
    PCargs.MaxNumPoints = maxnumpoints                # The following 3 parameters are set after trial-and-error
    PCargs.MaxStepSize  = maxstep
    PCargs.MinStepSize  = minstep
    PCargs.StepSize     = step
    PCargs.StopAtPoints = ['B']                       
    if LocBifPoints != None:
        PCargs.LocBifPoints = LocBifPoints            # detect limit points / saddle-node bifurcations
    PCargs.SaveEigen    = saveeigen                   # to tell unstable from stable branches
    return PCargs
#-------------------------------------------------------------------------------------------------------------------#
   
#-------------------------------------------------------------------------------------------------------------------#
def plot_continuation(ODE, PCargs, nmodel, keys, freepar, names, fs=[6,5], off_points=False, nrow=None, ncol=None, 
                      fontsize=18, c='k', dpi=500, save_fig=False, xmin=None, xlim=None, ylim=None, LimitPoints='',
                      showcurve=True, FPs=None, linewidth=3, normal_form_coef=False):
### Plot continuation curves
    if showcurve:
        if ncol == None:
            ncol = len(keys)
        if nrow == None:
            nrow = 1
        figure(figsize=(fs[0]*ncol,fs[1]*nrow), dpi=dpi)

    if FPs==None:
        FPs = [ODE.initialconditions]
        
    for j in range(len(FPs)):
        ODE.set(ics  = FPs[j])

        PyCont = PyDSTool.ContClass(ODE)     
        PyCont.newCurve(PCargs)
        PyCont[nmodel].forward()
        PyCont[nmodel].backward()
        if showcurve:        
            for i in range(len(keys)):
                PyCont.display((freepar,keys[i]), stability=True, axes=(nrow,ncol,i+1), color=c, linewidth=linewidth)
                if off_points:
                    PyCont.plot.toggleLabels('off')
                plt.xlabel(freepar, fontsize=fontsize)
                plt.ylabel(names[keys[i]], fontsize=fontsize)
                if xmin!= None:
                    plt.xlim(xmin=xmin)
                plt.title('')
                if xlim != None:
                    plt.xlim([xlim[0],xlim[1]])
                if ylim != None:
                    plt.ylim([ylim[0],ylim[1]])
    if save_fig:
        plt.savefig(save_fig, format='pdf', dpi=dpi)

    if normal_form_coef:
        i = 1
        while PyCont[nmodel].getSpecialPoint('LP'+str(i)):
            print "LP"+str(i), PyCont[nmodel].getSpecialPoint('LP'+str(i)).labels['LP']['data']
            i += 1
    
    if LimitPoints!='':
        P = []
        i=1
        while PyCont[nmodel].getSpecialPoint(LimitPoints+str(i)):
            P += [PyCont[nmodel].getSpecialPoint(LimitPoints+str(i))[freepar]]
            i +=1
        return P
#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#
def stability(FPs, ODE, eps=0.1):
### Determine the stability of a fixed point. eps=0.1
    out = []
    for i in range(len(FPs)):
        X = {}
        stable = True
        for k in FPs[0].keys():
            X[k] = FPs[i][k]*(1 + eps*rand.sample(list([-1,1]),1)[0])
        ODE.set(ics  = X)  
        traj = ODE.compute('traj')
        X = traj.sample()[-1]
        for k in FPs[0].keys():
            if np.abs(X[k]-FPs[i][k]) > eps*FPs[i][k]:
                stable = False
        out += ['S'] if stable else ['I']
    return out
#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#
def nullclines(axis, DSargs, stab, fp, names=None, nfp=0, vlim=None, c = ['b','g'], maxpoints=[1000,1000], step=5e+1, 
               minstep=1e-1, maxstep=1e+3, fs=[6,5], save_fig=False, plotaxis=[0,1], loc=0, fontsize=16,
               pcontour=None):
### Plot nullclines   
    figure(figsize=(fs[0],fs[1]), dpi=500)
    DSnc = copy.deepcopy(DSargs) 
    for i in plotaxis:
        keys = DSargs.varspecs.keys()
        keys.remove(axis[i])
        DSnc.pars[axis[i]] = fp[nfp][axis[i]]
            
        DSnc.varspecs = {}
        DSnc.ics = {}
        DSnc.xdomain = {}
        DSnc.pdomain = {}
        if vlim != None:
            DSnc.pdomain[axis[i]] = vlim[axis[i]]
        for k in keys:
            DSnc.varspecs[k] = DSargs.varspecs[k]
            DSnc.ics[k]      = fp[nfp][k]
            DSnc.xdomain[k]  = DSargs.xdomain[k]
        ODEnc = Vode_ODEsystem(DSnc)
        PCargs = PyCont_args('nullclines', axis[i], maxpoints[i], maxstep, minstep, step, LocBifPoints=None)
        PCargs.StopAtPoints = ['B']
        PyCont = PyDSTool.ContClass(ODEnc)
        PyCont.newCurve(PCargs)
        PyCont['nullclines'].forward()
        PyCont['nullclines'].backward()
        PyCont.display((axis[0],axis[1]), stability=True, linewidth=3, color=c[i], label='d'+names[axis[i-1]]+'/dt'+'=0'
                       if names!=None else 'd'+axis[i-1]+'/dt'+'=0' )
        PyCont.plot.toggleLabels('off')
        PyCont.plot.togglePoints('off')
        del DSnc.pars[axis[i]]

    for i in range(len(fp)):
        plt.plot(fp[i][axis[0]],fp[i][axis[1]], 'ok', markersize=12, markerfacecolor='r' if stab[i]=='S' else 'w')
        
    if pcontour!=None:  
        H, xedges, yedges = np.histogram2d(asarray(pcontour[axis[0]])[:,0], asarray(pcontour[axis[1]])[:,0], bins=100)
        H = np.rot90(H)
        H = np.flipud(H)
        xbin = 0.5*(xedges[1] - xedges[0])
        ybin = 0.5*(yedges[1] - yedges[0])
        plt.contour(xedges[1:]-xbin, yedges[1:]-ybin, H)
    plt.xlabel(names[axis[0]] if names != None else axis[0], fontsize=fontsize)
    plt.ylabel(names[axis[1]] if names != None else axis[1], fontsize=fontsize)
    plt.title('')
    plt.legend(loc=loc)
    if vlim != None:
        plt.xlim((vlim[axis[0]][0],vlim[axis[0]][1]))
        plt.ylim((vlim[axis[1]][0],vlim[axis[1]][1]))
    if save_fig:
        plt.savefig(save_fig, format='pdf', dpi=500)
#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#
def param_sensitivity_bars(list_pars, ODE, DSargs, var, save_fig=False, fs=[10,5], delta=[0.0, 0.1, -0.1]):
### Plot the changes in one variable as changes of the values of the paramters
    change = {}
    for pars in list_pars:
        if DSargs.pars[pars] != 0:
            a = []
            for d in delta:
                ODE.set(pars = {pars: (1.0 + d)*DSargs.pars[pars]} ) 
                a += [eliminate_redundants(pp.find_fixedpoints(ODE, n=4, maxsearch=1e+4, eps=1e-12),6)[0][var]]
            change[pars] = [100*(a[2] - a[0])/a[0],100*(a[1] - a[0])/a[0]] 
        else:
            change[pars] = [0,0] 
    l = change.keys()
    isort = np.argsort([np.abs(change[i][0])+np.abs(change[i][1]) for i in l])[::-1]

    figure(figsize=(fs[0],fs[1]), dpi=500)
    plt.bar(range(len(change.keys())), [change[l[i]][0] for i in isort], color='r', align='center', alpha=0.8)
    plt.bar(range(len(change.keys())), [change[l[i]][1] for i in isort], color='b', align='center', alpha=0.8)
    plt.xticks(np.arange(len(list_pars)+1), [l[i] for i in isort])
    plt.xlim([-1,len(list_pars)])
    plt.ylabel('Change in the signal (%)', fontsize= 18)
    plt.legend( ('- 10%', '+10%'), loc='upper right')
    if save_fig:
        plt.savefig(save_fig, format='pdf', dpi=900) 
    plt.show()
#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#
def param_sensitivity_bifurcations(DSargs, pars, freepar, var, nmodel, delta=[0.0, 0.1, -0.1], c=['k', 'b', 'r'], 
                                   ylim=None, save_fig=None, maxstep=1e+1, minstep=1e-3, step=5e+0):
### plot bifurcation curves for changes in 10% in the parameters 
    figure(figsize=(6*len(pars),5), dpi=500)
    for i in range(len(pars)):
        ODE = Vode_ODEsystem(DSargs)
        for j in range(len(delta)):
            ODE.set(pars = {pars[i]: (1.0 + delta[j])*DSargs.pars[pars[i]]} ) 
            fp_coord = eliminate_redundants(pp.find_fixedpoints(ODE, n=4, maxsearch=1e+4, eps=1e-10),6)
            ODE.set(ics  = fp_coord[0])  
            PCargs = PyCont_args(nmodel, freepar, 200, saveeigen=True, maxstep=maxstep, minstep=minstep, step=step)
            PyCont = PyDSTool.ContClass(ODE)     
            PyCont.newCurve(PCargs)
            PyCont[nmodel].forward()
            PyCont[nmodel].backward()
            PyCont.display((freepar,var), stability=True, axes=(1,len(pars),i+1), color=c[j], linewidth=3)
            plt.plot(0,0, linewidth=3, color=c[j])
            PyCont.plot.toggleLabels('off')
        plt.title(pars[i])
        if i == 0:
            plt.legend(('0%','+10%','- 10%'))
        if ylim != None:
            plt.ylim(ylim)
    if save_fig:
        plt.savefig(save_fig, format='pdf', dpi=900) 
#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#
def plot_PhaseDiagram(ODE, nmodel, names, freepar, vfreepar, par, rpar, xlim, ylim, figname, keys=['N','D','J','I'],
		      saveeigen=False, maxstep=5e-1, minstep=1e-2, step=1e-1, silence=False, LimitPoints='LP', showBif=True):
    if silence:
        class NullDevice():
            def write(self, s):
                pass
        original_stdout = sys.stdout
        sys.stdout = NullDevice()

    x  = []
    for i in rpar:
        ODE.set(pars = {freepar: vfreepar, par: i} )
        fp = eliminate_redundants(pp.find_fixedpoints(ODE, n=2, maxsearch=1e+3, eps=1e-10),6)
        ODE.set(ics  = fp[0]) 
        PCargs = PyCont_args(nmodel, freepar, 10000, saveeigen=saveeigen, maxstep=maxstep, minstep=minstep, step=step)
        p = plot_continuation(ODE, PCargs, nmodel,keys, freepar, names, FPs=[fp[0]], showcurve=showBif,
			      off_points=True, fontsize=20, xlim=xlim, LimitPoints=LimitPoints)
        for j in range(len(p)):
            x += [[i,p[j],j]]

    figure(figsize=(6,5), dpi=500)
    plt.xlabel(freepar, fontsize= 18)
    plt.ylabel(par, fontsize= 18)

    plt.xlim(xlim)
    plt.ylim(ylim)

    x = asarray(x)
    color = ['k', 'r', 'r', 'b']
    for i in range(int(max(x[:,2]))+1):
        plot(x[x[:,2]==i][:,1],x[x[:,2]==i][:,0], color=color[i])
    plt.savefig(figname, format='pdf', dpi=500)
    if silence:
        sys.stdout = original_stdout
#-------------------------------------------------------------------------------------------------------------------#