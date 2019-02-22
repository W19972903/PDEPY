# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 12:18:52 2019

@author: Michele Boreale
 
Python code for the examples in "Coalgebra, partial differential equations and boundary problems"
by Michele Boreale
February 2019
Experimented under the Anaconda distribution for Windows 10.

Summary: given a set of PDE systems ("stratified system"), representing a PDE boundary problem,
the DC algorithm computes (weakest) preconditions and (strongest) postconditions
as specified below (subject to certain restrictions, see paper).



########## USAGE ###################################################################################

_ = initvar(indepX=['t','x'],dependU=['u'],maxder=4,extrader=None)
Defines independent variables, and dependent variables and all their derivatives up to order maxder; these will be 
accessible in the global variable Xlist; optional argument extrader permits more flexible definitions of derivatives (see examples)

qt,J = DC(pt,Huser,P,Plist,Monotony=False)
INPUT:
pt = polynomial template (possibly defined with genpt, see examples)
Huser = a list of subsystems, defining the coherent, stratified system (boundary problem) to analyse. NB: coherence is not checked!
P = list of polynomials, defining the variety V(P) of initial boundary conditions. The set of ind. variables X (Xlist) will be included automatically.
Plist = list of parameters appearing in pt
Monotony = True assumes the chain of generated ideals is ascending, rather than just eventually ascending, leading to a more efficient algorithm. 
           This is safe only when pt has 0 or 1 parameter.
           
OUTPUT:
qt = result template with s parameters. The set of instantiations of qt, qt[R^s], is a postcondition of V(P). It is the *strongest* (largest) such postcondition in case <P>=I(V(P)) (i.e. <P> is a real radical)
J =  (Groebner basis of the) ideal defining the weakest  precondition of qt[R^s].
If maxiter iterations are exceeded, the execution stops with a warning.


########## SYNTAX ##################################################################################

- A jet-like notation for derivatives is used;
- dependent variables must be ONE-letter, say u,v,..;
- indices denote order of derivatives. E.g. if X={t,x} and U={u}, then u00 denotes u, u12 denotes u_{txx}, and so on;
- however, if u is specified via the flexible 'extrader' option (see example on heat eq. below), 
  indices are separated by an underscore '_': u0_0 denotes u, u1_2 denotes u_{txx}, and so on;
- an equation eq=(d,p) is a pair, where d is a derivative (e.g. u00, u12, ...) and p is a polynomial;
- a system S=[eq1,...,eqk] is a list of equations;
- a subsystem C=(S,Y), is a pair with S a system and Y a list of independent variables
- a stratified system H=[C1,..,Cm] is a list of subsystems. The main subsystem must be C1.


########## EXAMPLES IN THE PAPER ###################################################################

1) INVISCID BURGERS' EQUATION
A nonlinear equation plus a linear boundary condition at t=0. For generic real constants c,b:
u_t(t,x) = u(t,x)*u_x(t,x)
u(0,x)=c*x+b

We want to find all valid polynomial postconditions of total degree <=3. 
Try the following snippet.


_=initvar(['t','x'],['c','b','u'],maxder=1)
eq1b=[(u10,-u00*u01),(c10,0),(c01,0),(b10,0),(b01,0)]
eq2b=[(u01,c00)]
C1b=(eq1b,[t,x])
C2b=(eq2b,[x])
Hb=[C1b,C2b]  
P0b=[u00-b00]
pvb=[c00, b00, u00]
pt,pl=genpt([t,x]+pvb,3)  # generates complete polynomial template of total degree 3 with indeterminates {t,x,c00, b00, u00}
qt,_=DC(pt,Hb,P0b,pl)
Poly(qt/1,[u00])/1        # pretty printing of qt. Setting the only parameter of qt to 1 yields the equation u(t,x)=(c*x+b)/(c*t+1)



2) HEAT EQUATION
Heat equation in 1 spatial dimension, with a generic sinusoidal boundary condition at t=0 (b,c real constants): 
u_t(t,x)=b*u_xx(t,x)
u(0,x)=sin(c*x)

We look for a solution u of the form: (exponential of time) * (sinusoid of space)
Try the following snippet.

_=initvar(['t','x'],['a','b','c','d','f','g','h','j','i','u'],maxder=1,extrader={'u':':2_(:3)'})
C1=([(u1_0, b00*u0_2),(a10, 0),(a01, 0),(b10, 0),(b01, 0),(c10, 0),(c01, 0),(f10, 0),(g10, 0),(h01, 0),(i10, 0),(i01, 0),(j10, 0),(j01, 0),(d01, 0),(d10, 0)],[t, x])
C2=([(u0_1, g01), (f01, -c00*g00), (g01, c00*f00)], [x])
C3=([(h10, -d00*h00)], [t])
Heat= [C1, C2, C3]
a1,e=var('a1,e')  # dummy variables
pt=a1*a00*(u0_0+i00*g00*h00+j00*f00*h00) 
P0=[a00]  
_,J=DC(pt,Heat,P0,[a1],True)
J1=GroebnerBasis(list(J)+[a00-1,f00-1,g00,h00-1,c00*e-1],Xlist+[e],domain='QQ')  # build a set of new equations with specific values for a00,f00,g00,h00 and c00!=0
solve(J1,[i00,j00,d00,e],dict=True)  # replace the values obtained here for i00,j00,d00 in E, then solve for u



3) TRIGONOMETRIC FUNCTIONS

xtri=initvar(['x','y'],['f','g','i','j','h','k'],maxder=1)
eqtr=[(f10,-g00), (f01,-g00), (g10,f00), (g01,f00), (i10,-j00), (i01,0), (j10,i00), (j01,0),(h10,0), (h01,-k00), (k10,0), (k01,h00) ]
Htr=[(eqtr,[x,y])]
pt,pl=genpt([f00, g00, i00, j00, h00, k00],2) # f=cos(x+y), g=sin(x+y), i=cos(x), j=sin(x), h=cos(y), k=sin(y)
P0=[f00-1, g00, i00-1, j00, h00-1, k00]
qt,J=DC(pt,Htr,P0,pl,False)
Poly(qt/1,pl)/1    # pretty printing of result template
sigma={f00:cos(x+y), g00:sin(x+y), i00:cos(x), j00:sin(x), h00:cos(y), k00:sin(y)}
tridlist=[term.subs(sigma) for term in list(Poly(qt/1,pl).as_dict().values())]    # all trigonometric laws of degree <=2 (with the above positions f=cos(x+y) etc.)
print(tridlist)
####################################################################################################
"""



from sympy  import *
from sympy  import itermonomials as itm
import time
import numpy as np
from sympy.polys.orderings import monomial_key

global D
global DI
global NIND
global Xlist
 

maxiter = 60

def genpt(Xlist,deg):
    monlist=list(itm(Xlist,deg))
    l=len(monlist)
    parlist = list(var('a%d' % j) for j in range(l))
    prod=(Matrix(monlist).T *Matrix(parlist))[0,0]
    return prod, parlist

def linind(list,p):     # stub: should check linear independence of p from res
    return not(p==0)


def instantiate0(ptlist, zeropar,xPar):
    global Xlist, NIND
    res = []
    for pt in ptlist:
        if zeropar.keys()==[]:
            if linind(res,pt):
                res.append(pt)
        else:
            for a in zeropar.keys():
                tmp=zeropar.copy()
                tmp.update({a:1})
                p=Poly(pt.subs(tmp)/1,Xlist[:NIND]+xPar)
                if linind(res,p):
                    res.append(p)
    return res # returns list of polynomials in R



def check0(base, newinst,Xpar):
    global Xlist, NIND
    if Xpar==[]:
        return (newinst==0), base
    G = groebner(base,Xlist[:NIND]+Xpar,order='grevlex')
    lg=len(G)
    for q in newinst:
        if q!=0:
            if lg==0:
                return False,G
            #print('q=',q)
            l,r = reduced(q,G)
            if r!=0:
                return False,G
    return True,G
    

def seqsub(sigma,pt):   # sequentially applies a list og substitutions in sigma to template pt
    for s in sigma:
        pt=pt.subs(s)
    return pt

def ext2int(Huser):
    global Xlist
    return [(formeq(C[0]),C[1]) for C in Huser]
    
def checkinvariance(Jm,freepar):  # Stub: should check invariance: S_H(delta(E,x)) in Jm, for each E in Hashtable (visited), x in X
    return None
    
    
def DC(pt,Huser,P,Plist,Monotony=False):   
    start_time = time.time()
    global D, DI, NIND, HT, Xlist
    H=ext2int(Huser) 
    pt=Poly(pt/1,Xlist)
    Xind=Xlist[:NIND]
    zeroind={x:0 for x in Xind}
    Xpar=parvar(H)
    G0=groebner(P,Plist+Xpar,domain='QQ',order='lex')  
    print('Groebner basis for <P0,X> = ',Xlist[:NIND]+[p/1 for p in list(G0)])
    print("Search will proceed by exploring derivatives of increasing order of input template pt (\"levels\")")
    print("")
    if Monotony:
        print("Monotony = True: *monotonic search optimization* will be used.")
        print("Monotonic search keeps track of monomials tau s.t. the derivative pt_tau is in the ideal of previously visited derivatives (\"frontier\"). ")
        print("Stops as soon as there is no monomial left to explore below the frontier.")
        print("Warning: monotonic search oprimization is heuristic. Invariance of {pt_tau} (available as global dictionary variable HT) should be checked ex-post.")
    print("")
    qt=onestepstarHpt(pt,H)       # compute S_H(pt_0)
    if qt==0:
        rt=0
    else:
        _,rt=reduced(Poly(qt.subs(zeroind),Plist+Xpar),G0)  # rt used to extract linear constraints for membership in JR
    border=[]
    print('*** Level m=0 ***')
    if Monotony:
        print('  Frontier=',border)
    print('  tau   =1')
    #print("pt0  =                   ",pt/1)
    #print("qt0  = S_H(pt0) =        ",qt/1)
    #print("rt0  = S_H(pt0) mod J  = ",rt/1)

    coeffs = Poly(rt/1,Xpar).coeffs() # list of linear expressions (constraints) extracted as coefficients of rt    
    sigma=[]    # list of individual substitutions for parameters 
    for lin in coeffs:      # solve the list of linear constraints, one by one
        lin= (seqsub(sigma,lin))  # apply previously accumulated subst. to current constraint lin
        lin=Add(lin,0)
        if lin.free_symbols==set():
            if lin==0:
                C=True
            else:
                sigma=sigma+[{lin:0}]
                print("  Linear constraint for V_0: "+str(sigma))
                print("  Cannot find instance of given template that is an invariant.")
                print('m=0')
                print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                return False
        else:
            C=solve(lin,list(lin.free_symbols),solution_dict=True)      # solve one linear constraint lin
            s = C[0]
            if type(s)==dict:
                sigma=sigma+[s]
            else:
                sigma=sigma+[{list(lin.free_symbols)[0]:s}]                                      # append s to list of substitutions
    print("  New linear constraint for V_0: "+str(sigma))        

    pt=Poly(seqsub(sigma,pt),Xlist)   # apply sigma to pt
    qt=Poly(seqsub(sigma,qt),Xlist)   # apply sigma to qt
    derlist=[qt] # derlist contains {S_H(pt_0),...,S_H(pt_m)}; used to generate B_m after instantiation of the parameters a_j 
    freepar = pt.free_symbols-set(Xlist) 
    zeropar = {a:0 for a in freepar}
    base = []
    updbase = False
    m=(0,)*len(Xind) # m is the last generated monomial
    initHT(len(Xind),pt)
    levelchecked=False
        
    for j in range(maxiter):
        print("")
        oldeg=sum(m)
        m,newpt=nextder(m,H[0],border,Monotony)  # newpt = new total derivative of pt
        if m==None:
            print("Frontier =",border)
            print('No new derivative below frontier: both chains stabilized')
            print('m='+str(oldeg))
            print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
            if Monotony:
                checkinvariance(Jm,freepar)   # Shoud check invariance ex-post; not yest implemented
            return (HT[(0,)*len(Xind)],GroebnerBasis(Xlist[:NIND]+list(Jm), Xlist[:NIND]+Xpar))
        if (not(Monotony)):
            if ((sum(m)>oldeg) & (not(Monotony))):
                if levelchecked:
                    print('Entire level checked: both chains stabilized')
                    print('m='+str(sum(m)))
                    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                    return (HT[(0,)*len(Xind)],GroebnerBasis(Xlist[:NIND]+list(Jm), Xlist[:NIND]+Xpar))
                else:
                    levelchecked=True
        if sum(m)>oldeg:
            print('*** Level m=',str(sum(m))+' ***')
        if Monotony:
            print('  Frontier=',border)        
        print('  tau   =',Poly({m:1},Xind)/1)
        #print("pt"+str(j+1)+"= ",newpt/1) 
        qt=onestepstarHpt(newpt,H)   # qt=S_H(newpt) 
        #print("qt"+str(j+1)+"= ",qt/1) 
        if qt==0:
            rt=0
        else:
            _,rt=reduced(Poly(qt.subs(zeroind),Plist+Xpar),G0)     # rt= S_H(newpt) mod G0
        #print("rt"+str(j+1)+"= ",rt/1)
        coeffs = Poly(rt/1,Xpar).coeffs() # list of linear expressions (constraints) extracted as coefficients of rt    
        sigma=[]                                                
        for lin in coeffs:      # solve the list of linear constraints, one by one
            lin=seqsub(sigma,lin)  # apply previously accumulated subst. to current constraint lin
            lin=Add(lin,0)           # forces lin to symbolic expression
            if lin.free_symbols==set():
                if lin==0:
                    C=True
                else:
                    sigma=sigma+[{lin:0}]
                    print("  New linear constraint for V_"+str(sum(m))+": "+str(sigma))
                    print("Cannot find instance of given template that is an invariant.")
                    print('m='+str(j+1))
                    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                    return False
            else:
                C=solve(lin,list(lin.free_symbols),solution_dict=True)      # solve one linear constraint lin
                s = C[0]
                if type(s)==dict:
                    sigma=sigma+[s]
                else:
                    sigma=sigma+[{list(lin.free_symbols)[0]:s}]                                      # append s to list of substitutions
        print("  New linear constraint for V_"+str(sum(m))+": "+str(sigma))

        if sigma==[]:
            print("  No new linear constraint detected for V_"+str(sum(m)))
            print("  Checking equality with previous ideal...")
            if not(updbase):
                zeropar = { a:0 for a in freepar }
                base = instantiate0(derlist,zeropar,Xpar)
                updbase = True
            newinst=instantiate0([qt],zeropar,Xpar)
            flag,Jm=check0(base,newinst,Xpar)
            #print('base=',[t/1 for t in base])
            if flag:
                if Monotony:
                    print("  Equality holds, updating frontier.")
                else:
                    print("  Equality holds.")
                border.append(m)
            else:
                print("  Equality does not hold, chains not yet stabilized; level not yet cleared")
                levelchecked=False                
            derlist.append(qt)
            base = base + newinst
        else:
            updbase=False
            levelchecked=False
            derlist = [Poly(seqsub(sigma,qtt)/1,Xlist) for qtt in derlist]
            for mon in HT.keys():
                ptm=HT[mon]
                HT[mon]=Poly(seqsub(sigma,ptm)/1,Xlist)
            qt = Poly(seqsub(sigma,qt)/1,Xlist)
            if qt!=0:
                derlist.append(qt)
            freepar = freepar-set([ list(s.keys())[0] for s in sigma]) 
    print('m='+str(sum(m)))
    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))    
    return (HT[(0,)*len(Xind)],base,'WARNING: maximum number of iterations reached, result is likely to be not correct.')

    
    
 

  



##### PDE

# generate dependent var, their derivatives and the corresponding dictionary and complete var list
def initvar(indepX=['t','x'],dependU=['u'],maxder=4,extrader=None):  # indepX,dependU must be disjoint lists of chars; ex: extrader={'u':':11_(:3)'}
    global D, DI, NIND, Xlist
    NIND=len(indepX)
    derlist=[]
    for u in dependU:
        if extrader!=None:
            if u in extrader.keys():
                dimstr=extrader[u]
            else:
                dimstr=('(:'+str(maxder+1)+')')*NIND
        else:
            dimstr=('(:'+str(maxder+1)+')')*NIND
        derlist=derlist+list(var(u+dimstr))
    D={}
    DI={}
    for der in derlist:
        derstr=str(der)
        if extrader!=None:
            if ((derstr[0] in extrader.keys()) & ('_' in derstr)):
                D[NIND+derlist.index(der)]=[int(c) for c in derstr[1:].split('_')]+[derstr[0]]
            else:
                D[NIND+derlist.index(der)]=[int(c) for c in list(derstr[1:])]+[derstr[0]]
        else:
            D[NIND+derlist.index(der)]=[int(c) for c in list(derstr[1:])]+[derstr[0]]
    DI={tuple(D[k]):k for k in D.keys()}
    Xlist=var(indepX)+derlist    
    return Xlist  

def formeq(E):
    global Xlist
    d={}
    for (u,q) in E:
        d[tuple(D[Xlist.index(u)])]=q
    return d
    
def mon2list(alpha,Xind):
    m=Poly(alpha,Xind).monoms()[0]
    taulist=[[Xind[j]]*m[j] for j in range(len(m))]
    tau= [item for sublist in taulist for item in sublist] # flatten
    return tau
    
    
        
def dermon(m,x):
    global D, DI, NIND, Xlist
    res={}
    changexl=False
    j=Xlist.index(x)
    for i in np.nonzero(m)[0]:
        mi=m[i]    # degree of i-th variable
        mnew=list(m)+[0]*(len(Xlist)-len(m))  # create new monomial, with padding
        if (i>=NIND):  # dependent variable
            mnew[i]=mi-1  # decrease degree of variable
            dxiplus=list(D[i])   # copy dictionary's value to variable
            dxiplus[j]=dxiplus[j]+1
            dxiplust=tuple(dxiplus)
            if not(dxiplust in DI.keys()):   # create new variable, update Xlist, D, DI
                l=len(dxiplus)
                xnew=var(''.join(list([dxiplus[l-1]]+['_'+str(c) for c in dxiplus[:l-1] ]))  )
                Xlist.append(xnew)
                lx=len(Xlist)#-NIND
                D[lx-1]=dxiplus
                DI[dxiplust]=lx-1
                mnew=mnew+[0]
                changexl=True
            h=DI[dxiplust]
            mnew[h]=mnew[h]+1
            mnew=tuple(mnew)
            if mnew in res.keys():
                res[mnew]=res[mnew]+mi
            else:
                res[mnew]=mi
        elif (i==j):            # x variable
            mnew=list(m)
            mnew[i]=mi-1
            mnew=tuple(mnew)
            if mnew in res.keys():
                res[mnew]=res[mnew]+mi
            else:
                res[mnew]=mi
    return res, changexl
                
            
def sumdict(LC,changexl):
    global Xlist
    lx=len(Xlist)
    if changexl:   # padding needed
        LD=[]
        for d in LC:
            d1={}          
            for k in d.keys():
                lk=len(k)
                d1[tuple(list(k)+[0]*(lx-lk))]=d[k]
            LD.append(d1)
    else:
        LD=LC
    sumd={}
    for d in LD:
        for k in  d.keys():
            if k in sumd.keys():
                sumd[k]=sumd[k]+d[k]
            else:
                sumd[k]=d[k]
    return sumd
            
        
def totalder(p,x):
    global Xlist
    pdic=p.as_dict()
    LC=[]
    chxl=False
    for m in pdic.keys():
        v=pdic[m]
        dmon,changexl=dermon(m,x)
        LC.append({k:dmon[k]*v for k in dmon.keys()})
        if changexl:
            chxl=True
    sd=sumdict(LC,chxl)
    return Poly(sd,Xlist) 

def totalderstar(p,tau):
    global Xlist
    q=Poly(p/1,Xlist)
    for x in tau:
        q=totalder(q,x)
    return q
    
def searchprinc(i,Comp):
    global D, NIND, Xlist
    if not(i in D.keys()):
        return None
    der=tuple(D[i])
    Yind=[Xlist.index(y) for y in Comp[1]] # independent variables indices
    Eqns=Comp[0]
    for k in Eqns.keys():
        if k[NIND]==der[NIND]:
            if (all(x <= y for x, y in zip(k,der))):
                diff=[der[j]-k[j] for j in range(NIND)]
                if set(np.nonzero(diff)[0]).issubset(Yind):
                    return diff, Eqns[k]
    return None
    
def onestep(p,Comp):  # important, p MUST be a polynomial with Xlist as generators
    global NIND, Xlist
    mons=p.as_dict()
    for m in mons:
        for i in [j for j in np.nonzero(m)[0] if j>=NIND]:
            search=searchprinc(i,Comp)
            if search!=None:
                tau= [item for sublist in [[Xlist[j]]*search[0][j] for j in range(NIND)] for item in sublist]      
                G=totalderstar(search[1],tau)/1
                q=p/1
                return Poly(q.subs({Xlist[i]:G}),Xlist)#q.subs({Xlist[i]:G})#Poly(q.subs({Xlist[i]:G}),Xlist)
    return None

            
def onestepstar(p,C): # important, p MUST be a polynomial with Xlist as generators
    global Xlist    
    while (True):
        q=onestep(p,C)
        if q==None:
            return p
        p=q
                
def onestepH(p,H): # important, p MUST be a polynomial with Xlist as generators
    global NIND, Xlist
    for C in H:
        q=onestepstar(p,C)
        if q!=p:
            return q    
    return None
    
def onestepstarH(p,H):
    global Xlist
    while (True):
        q=onestepH(p,H)
        if q==None:
            return p
        p=q


def onestepstarHpt(pt,H):
    global Xlist
    ptdic=pt.as_dict()
    newpt=0
    for m in ptdic:
        m1=tuple(list(m)+[0]*(len(Xlist)-len(m))) # padding
        e=onestepstarH(Poly({m1:1},Xlist),H)/1
        newpt=newpt+ptdic[m]*e
    return Poly(newpt,Xlist)
    
        
        
def delta(p,x,C):  # important, p MUST be a polynomial with Xlist as generators
    global Xlist
    p=Poly(p/1,Xlist)
    q=onestepstar(p,C)
    q1=totalder(q,x)
    return onestepstar(q1,C)
    
def deltapt(pt,x,C):  # important, pt MUST be a polynomial with Xlist as generators
    global Xlist
    pt=Poly(pt/1,Xlist)    
    ptdic=pt.as_dict()
    newpt=0
    for m in ptdic:
        m1=tuple(list(m)+[0]*(len(Xlist)-len(m))) # padding
        e=delta(Poly({m1:1},Xlist),x,C)/1
        newpt=newpt+ptdic[m]*e
    return Poly(newpt,Xlist)
    
def deltastar(p,tau,C):
    global Xlist
    q=p
    for x in tau:
        q=delta(q,x,C)
    return q
    
def parvar(H):
    global Xlist
    S=[]
    found=False
    for k in D.keys():
        for C in H:
            res=searchprinc(k,C)
            if res!=None:
                found=True
                break
        if found:
            found=False
        else:
            S.append(k)
    return [Xlist[k] for k in S]

def fact(m,xind):
    mt=Poly(m,xind).monoms()[0]
    prod=1
    for i in mt:
        prod=prod*factorial(i)
    return prod

def taylor(p,H,n=3,mord='grlex'):  # important, p MUST be a polynomial with Xlist as generators
    global NIND,Xlist
    xInd=Xlist[:NIND]
    Hi=ext2int(H)
    mons=sorted(itermonomials(xInd, n), key=monomial_key(mord, xInd))
    coefflist=[ onestepstarH(deltastar(p,mon2list(m,xInd),Hi[0]),Hi)/fact(m,xInd) for m in mons ] 
    return sum([c*m for c,m in zip(coefflist,mons) ])
    

def nextord(m):  # next in graded lex order
    n=len(m)
    totaldeg=sum(m)
    if totaldeg==0:
        return tuple([0]*(n-1)+[1])
    nz=np.nonzero(m)[0]
    j=max(nz)
    if j==0:
        return tuple([0]*(n-1)+[totaldeg+1])
    mnew=list(m)
    mnew[j-1]=mnew[j-1]+1
    mnew[j]=mnew[j]-1 
    return tuple(mnew)

def onediff(m): # smaller monomial differing of deg-1
    nz=np.nonzero(m)[0]
    j=min(nz)
    mnew=list(m)
    mnew[j]=mnew[j]-1 
    return tuple(mnew), j

global HT

   
def initHT(n,p):
    global HT
    HT={}
    HT[(0,)*n]=p



def nextder(m,C,border,Monotony=True):
    global HT,Xlist
    if Monotony:
        mnew=nextordborder(m,border)
    else:
        mnew=nextord(m)
    if mnew==None:
        return None, None
    moneminus,j=onediff(mnew)
    x=Xlist[j]
    pt=HT[moneminus]
    ptnew=delta(pt,x,C)#deltapt(pt,x,C)
    HT[mnew]=ptnew
    return mnew, ptnew

def nextordborder(m,border):
    if border==[]:
        return(nextord(m))
    mnew=m
    newlevel=False
    while(true):
        mnew1=nextord(mnew)
        if (sum(mnew1)>sum(mnew)): # mnew1 is one level up mnew's (total degree increased by 1)
            if newlevel:     # if this is the second total degree increase
                return None  # then traversed an entire degree level without finding any derivative extra border
            newlevel=True
        mnew=mnew1
        if ( all(  (any(x <y for x, y in zip(mnew,m1)))   for m1 in border   )  ):   # mnew1 is below the border
            return mnew
    





