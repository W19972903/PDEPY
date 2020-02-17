# -*- coding: utf-8 -*-
"""
Last modified on February 15, 2020, 9:44

@author: Michele Boreale

Python code for the examples in the paper 
"Automatic pre- and postconditions for partial differential equations"
by Michele Boreale 

Summary
Given a set of PDE systems ("stratified system"), the POST algorithm computes (weakest) 
preconditions and (strongest) postconditions as specified below (subject to certain restrictions; 
see paper).

Code tested under the Python Anaconda distribution for Windows 10.

########## USAGE ###################################################################################

From the Python console, import the script with:

from PDE import *

Main procedures.

_ = initvar(indepX=['t','x'],dependU=['u'],maxder=4,extrader=None)
Introduce independent variables, and dependent variables with all their derivatives up to order maxder; 
these will be accessible in the global variable Xlist; optional argument extrader permits more flexible 
definitions of derivatives (see examples).

qt,J = post(pt,H,P,Plist,Monotony=False)
INPUT:
pt = polynomial template (possibly defined with genpt, see examples)
H  = list of subsystems, defining the coherent (cf. paper), stratified system (initial value problem) 
     to analyse (NB: coherence is not checked)
P  = list of polynomials, defining the variety V(P) of initial conditions
Plist = list of parameters appearing in pt
Monotony = True: monotonic search optimization will be used, leading to a more efficient algorithm. 
           This option is *heuristic*.
           
OUTPUT:
qt = result template with s parameters. The set of instantiations of qt, qt[R^s], is a postcondition
     of V(P). It is the *strongest* (largest) such postcondition in case <P>=I(V(P))
J =  (Groebner basis of the) ideal defining the weakest  precondition of qt[R^s].
If maxiter iterations are exceeded, the execution halts with a warning.

OTHER useful functions.
pt,pl=genpt(vars,n)  : generates the complete polynomial template pt of total degree n 
                       with indeterminates in the list vars, and the corresponding list of 
                       parameters pl
taylor(p,H,n=3,mord='grlex') : Taylor expansion of the differential polynomial p up to order n.

NT,T=findCL(H,P0=[],base=[],degree=2) : 
   NT (nontrivial) + T (trivial) is a basis of all polynomial Conservation Laws (density-flux pairs) 
   of deg<=degree, with indeterminates in the given base, for the IVP H with initial data in V(P0).
   Works only for the case of 2 independent variables.


########## SYNTAX ##################################################################################

- Dependent variables must be 1-letter, say u,v,...
- A jet-like notation for derivatives is used: indices denote order of derivatives. 
  E.g. if X={t,x} and U={u}, then u00 denotes u, u12 denotes u_{txx}, and so on.
- However, if u is specified via the flexible 'extrader' option (see example on heat eq. below), 
  indices are separated by an underscore '_': u0_0 denotes u, u1_2 denotes u_{txx}, and so on.
- An equation eq is a pair (d,p), where d is a derivative (e.g. u00, u12, ...) and p is a polynomial.
- A system S is a list of equations [eq1,...,eqk].
- A subsystem C is a pair (S,Y), with S a system and Y a list of independent variables.
- A stratified system H is a list of subsystems [C1,..,Cm]. Here C1 must be the only main subsystem.

########## EXAMPLES IN THE PAPER ###################################################################

1) INVISCID BURGERS' EQUATION
A nonlinear equation plus a linear boundary condition at t=0. For generic real constants c,b:
u_t(t,x) = -u(t,x)*u_x(t,x)
u(0,x)=c*x+b

We want to find all valid polynomial postconditions of total degree <=3. 
Try the following snippet.

_=initvar(['t','x'],['c','b','u'],maxder=1)
eq1b=[(u10,-u00*u01),(b10,0),(c10,0),(c01,0)]
eq2b=[(u01,b00),(b01,0)]
C1b=(eq1b,[t,x])
C2b=(eq2b,[x])
Hb=[C1b,C2b]  
P0b=[u00-c00]
pvb=[t,x,c00, b00, u00]
pt,pl=genpt(pvb,3)  # generates complete polynomial template of total degree 3 with indeterminates {t,x,b,c,u}
qt,J=post(pt,Hb,P0b,pl)
Poly(qt,pl)/1        # pretty printing of qt. Setting the only parameter of qt to 1 yields the equation u(t,x)=(c*x+b)/(c*t+1)
solve(qt,[u00],dict=True)   # finds an explicit formula for the solution u


2) HEAT EQUATION
Heat equation in 1 spatial dimension, with a generic sinusoidal boundary condition at t=0 (b,c real constants): 
u_t(t,x)=b*u_xx(t,x)
u(0,x)=sin(c*x)

We look for a solution u of the form: (exponential of time) * (sinusoid of space)
Try the following snippet.

_=initvar(['t','x'],['a','b','c','d','f','g','h','j','i','u'],maxder=1,extrader={'u':':2_(:3)'})
C1=([(u1_0, b00*u0_2),(a10, 0),(a01, 0),(b10, 0),(b01, 0),(c10, 0),(d01, 0),(f10, 0),(g10, 0),(h01, 0),(i10, 0),(i01, 0),(j10, 0),(j01, 0)],[t, x])
C2=([(u0_0, g00), (f01, -c00*g00), (g01, c00*f00), (c01,0)], [x])
C3=([(h10, -d00*h00),(d10, 0)], [t])
Heat= [C1, C2, C3]
a1,e=var('a1,e')  # dummy variables
pt=a1*a00*(u0_0+i00*g00*h00+j00*f00*h00) 
P0=[a00]  
_,J=post(pt,Heat,P0,[a1])
J1=GroebnerBasis(list(J)+[a00-1,f00-1,g00,h00-1,c00*e-1],Xlist+[e],domain='QQ')  # build a set of new equations with specific values for a00,f00,g00,h00 and c00!=0
solve(J1,[i00,j00,d00,c00,e],dict=True)  # replace the values obtained here for i00,j00,d00 in E, then solve for u


3) A BOUNDARY PROBLEM
(u_x)^2+(u_y)^2=1 (Eikonal equation)
u=0 on the unit circle

A standard transformation by the method of characteristics
makes this equivalent to an IVP in the (r,s)-variables

x_s=2p
y_s=2q
z_s=2p^2+2q^2
p_s=0
q_s=0
with  initial conditions the s=0 line:
x(r,0)=cos(r), y(r,0)=sin(r), z(r,0)=0, p(r,0)=cos(r), q(r,0)=sin(r).

Try the following snippet.

_=initvar(['r','s'],['x','y','z','p','q'],maxder=1)
S1=[(x01,2*p00),(y01,2*q00),(z01,2*p00**2+2*q00**2),(p01,0),(q01,0)]
S2=[(x10,-q00),(y10,p00),(z10,0),(p10,-q00),(q10,p00)]
HC=[(S1,[r,s]),(S2,[r])]
P0=[z00,x00-1,y00,p00-1,q00]
pt,pl=genpt([x00,y00,z00],2)
qt,J=post(pt,HC,P0,pl)
solve(Poly(qt,pl).coeffs(),[z00],dict=True)   # finds an explicit formula for the solution z=u


4) CONSERVATION LAWS FOR A WAVE EQUATION IVP
u_tt=u_xx, u_t(0,x)=0, u(0,x)=A*sin(x)+B*cos(x)  (A,B arbitrary constants)

_=initvar(['t','x'],['u','v','w'],maxder=1,extrader={'u':':3_(:3)'})
Hw=[([(u2_0, u0_2),   (v10, 0),   (w10, 0)], [t, x]),
    ([(u1_0, 0), (u0_1, w00), (v01, -w00), (w01, v00)], [x]) ]
NT,T=findCL(Hw,P0=[],base=[u0_0,u1_0,u0_1,u1_1,u0_2],degree=2) # NT+T is a basis of all polynomial C.L.s (density-flux pairs) of degree <=2 in the given base of variables

Changing the first intial condition to u_t(0,x)=C*exp(-x^2) (C arbitrary constant)

_=initvar(['t','x'],['u','v','w','h'],maxder=1,extrader={'u':':3_(:3)'})
Hw=[([(u2_0, u0_2),   (v10, 0),   (w10, 0), (h10, 0)], [t, x]),
    ([(u1_0, h00), (u0_1, w00), (v01, -w00), (w01, v00), (h01, -h00*2*x)], [x]) ]
NT,T=findCL(Hw,P0=[],base=[u0_0,u1_0,u0_1,u1_1,u0_2],degree=2) 


########## ADDITIONAL EXAMPLES ###################################################################

5) ONE- AND TWO-VARIABLE TRIGONOMETRIC FUNCTIONS
Here f=cos(x+y), g=sin(x+y), i=cos(x), j=sin(x), h=cos(y), k=sin(y)

xtri=initvar(['x','y'],['f','g','i','j','h','k'],maxder=1)
eqtr=[(f10,-g00), (f01,-g00), (g10,f00), (g01,f00), (i10,-j00), (i01,0), (j10,i00), (j01,0),(h10,0), (h01,-k00), (k10,0), (k01,h00) ]
Htr=[(eqtr,[x,y])]
pt,pl=genpt([f00, g00, i00, j00, h00, k00],2) # f=cos(x+y), g=sin(x+y), i=cos(x), j=sin(x), h=cos(y), k=sin(y)
P0=[f00-1, g00, i00-1, j00, h00-1, k00]
qt,J=post(pt,Htr,P0,pl)
Poly(qt,pl)/1    # pretty printing of result template
sigma={f00:cos(x+y), g00:sin(x+y), i00:cos(x), j00:sin(x), h00:cos(y), k00:sin(y)}
tridlist=[term.subs(sigma) for term in list(Poly(qt/1,pl).as_dict().values())]    # all trigonometric laws of degree <=2 (with the above positions f=cos(x+y) etc.)
print(tridlist)


6) WAVE EQUATION
u_tt =b*u_xx 
u_t(0,x)=0       # initial condition 1
u(0,x)=sin(c*x)  # initial condition 2
h''(x)=s*h(x)    # generic exponential, spatial
k''(t)=-r*k(t)   # generic sinusoid, temporal

Determine, if any, conditions on s,r and i s.t. u(t,x)=i*h(x)*k(t)
Try the following snippet.

_=initvar(['t','x'],['a','b','c','f','g','i','h','k','s','r','u'],maxder=1,extrader={'u':':3_(:3)','h':':2_(:3)','k':':3_(:2)'})
Hw=[([(u2_0, b00**2*u0_2),
   (a10, 0),
   (a01, 0),
   (b10, 0),
   (b01, 0),
   (c10, 0),
   (i10, 0),
   (i01, 0),
   (f10, 0),
   (g10, 0),
   (h1_0, 0),
   (k0_1, 0),
   (r01, 0),
   (s10, 0)],
  [t, x]),
 ([(u1_0, 0),
   (u0_1, c00*f00),
   (f01, -c00*g00),
   (g01, c00*f00),
   (h0_2, h0_0*s00), (c01, 0), (s01, 0) ],
  [x]),
 ([(k2_0, -k0_0*r00) ,    (r10, 0) ], [t]), ]
P0 = [a00, f00-1,g00]
a1,e,l=var('a1,e,l')  # dummy variables
pt=a1*a00*(u0_0+i00*h0_0*k0_0)
_,J=post(pt,Hw,P0,[a1])
J1=GroebnerBasis(list(J)+[a00-1,f00-1,g00,c00*e-1,b00*l-1],Xlist+[e,l],domain='QQ',order='grevlex')  # build a set of new equations with specific values for a00,f00,g00,h00 and c00!=0
solve(J1,[r00,s00,i00,e,l],dict=True)  # This gives all possible nontrivial solutions
 
 
7) KdV 
u_t+u*u_x+u_xxx = 0
u(0,x)=c*x+b   # initial condition

_=initvar(['t','x'],['c','b','u'],maxder=1,extrader={'u':':2_(:4)'})
E1= [(u1_0,-(u0_0*u0_1+u0_3)), (c10,0), (b10,0), (b01,0)] 
E2= [(u0_1,c00),(c01,0)] # u(0,x)=c*x+b
E3= [(u0_1,-c00*u0_0),(c01,0)]  # u(0,x)=exp(-c*x)
Hkdv=[(E1,[t,x]), (E2,[x])]
Hkdv2=[(E1,[t,x]), (E3,[x])]
P0=[u0_0-b00]
P1=[u0_0-1]
pt,pl=genpt([u0_0,u0_1,u1_0,u1_1],2)
qt,J=post(pt,Hkdv,P0,pl)
Poly(qt,pl)/1  # pretty printing of qt, gives all valid polynomial equations of degree <= 2 for the considered variables.


8) Boundary problem n. 2
(u+2*y)*u_x+u*u_y=0
u(x,1)=1/x 

_=initvar(['r','s'],['x','y','z'],maxder=1)
S1=[(x01,z00+2*y00),(y01,z00),(z01,0)]
S2=[(x10,1),(y10,0),(z10,-z00**2)]
HC2=[(S1,[r,s]),(S2,[r])]
P0=[x00-1,y00-1,z00-1]
pt,pl=genpt([x00,y00,z00],2)
qt,J=post(pt,HC2,P0,pl)
solve(Poly(qt,pl).coeffs(),[z00],dict=True)   # finds an explicit formula for the solution z=u


9) Boundary problem n. 3
(y+u)*u_x+y*u_y=x-y
u(x,1)=1+x

_=initvar(['t','s'],['x','y','z'],maxder=1)
S1=[(x10,z00+y00),(y10,y00),(z10,x00-y00)]
S2=[(x01,1),(y01,0),(z01,1)]
HC3=[(S1,[t,s]),(S2,[s])]
P0=[x00,y00-1,z00-1]
pt,pl=genpt([x00,y00,z00],2)
qt,J=post(pt,HC3,P0,pl)
solve(Poly(qt,pl).coeffs(),[z00],dict=True)   # finds an explicit formula for the solution z=u

####################################################################################################
"""

#from sage.all import *

from sympy  import *
from sympy  import itermonomials as itm
import time
import numpy as np
from collections import OrderedDict
from sympy.polys.orderings import monomial_key

global D
global DI
global NIND
global Xlist

from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout
            

maxiter = 60

def genpt(Xlist,deg,offset=0):
    monlist=list(itm(Xlist,deg))
    l=len(monlist)
    parlist = list(var('a%d' % j) for j in range(offset,l+offset))
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
    if type(base)==type([]):
        G = groebner(base,Xlist[:NIND]+Xpar,order='grevlex',domain='QQ')#order='grevlex',field=True)
    else:
        G = base # assume base is already a Groebner basis
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

#from copy import *

def check01(base, newinst,Xpar):
    global Xlist, NIND
    if Xpar==[]:
        return (newinst==0), base
    if type(base)==type([]):
        G=groebner(base,Xlist[:NIND]+Xpar,order='lex',domain='QQ')#field=True)
    else:
        G = base # assume base is already a Groebner basis
    lg=len(G)
    for q in newinst:
        if q!=0:
            if lg==0:
                return False,q
            #print('q=',q)
            l,r = reduced(q,G)
            if r!=0:
                #print(G,q)
                return False,q
    return True,None
    
    
def checkinvariance(J):  # ex-post invariance check: S_H(delta(pt_tau,x)[v]) in J, for each tau,x and v in Vm
    global HT, Xlist, H, NIND,zeropar,Xpar    
    #zeropar={a:0 for a in freepar}
    HTkeys=list(HT.keys())
    F=frontier(HTkeys)
    print('Frontier of visited derivation monomials=',F)
    for m in F:
        pt=HT[m]
        for x in Xlist[:NIND]:
            newpt=delta(pt,x,H[0])#nextder(m,H[0],[],Monotony=False)
            qt=onestepstarHpt(newpt,H)
            qtinst=instantiate0([qt], zeropar,Xpar)
            flag, q = check01(J,qtinst,Xpar)
            if flag==False:
                print('Invariance is not true:')
                print('  tau=', Poly({m:1},Xlist[:NIND])/1)
                print('  x=', x)
                print('  pt_{tau}=', pt/1)
                print('  pt_{tau,x}=', newpt/1)
                print('  S_H(pt_{tau,x})=', qt/1)
                print('  S_H(pt_{tau,x})[v]=', q/1)
                return False
    print('Invariance is true.')
    return True
        


    
    
def post(pt,Huser,P,Plist,Monotony=False,extraptlist=None):   
    start_time = time.time()
    global D, DI, NIND, HT, Xlist, H, zeropar,Xpar
    H=ext2int(Huser) 
    pt=Poly(pt,Xlist)
    Xind=Xlist[:NIND]
    zeroind={x:0 for x in Xind}
    Xpar=parvar(H)
    G0=groebner(Xind+P,Plist+Xind+Xpar,order='lex',domain='QQ')#field=True)  
    print('Groebner basis for <P0 U X> = ',Xlist[:NIND]+[p/1 for p in list(G0)])
    print("Search will proceed by exploring derivatives of increasing order of input template pt (\"levels\")")
    print("")
    if Monotony:
        print("Monotony = True: *monotonic search optimization* will be used.")
        print("Monotonic search keeps track of the set of monomials tau s.t. pt_tau[Vm] is included in the ideal of previously visited derivatives (\"frontier\"). ")
        print("Stops as soon as there is no monomial left to explore below the frontier.")
        print("*Warning*: monotonic search optimization is heuristic; actual invariance of obtained ideal Jm (see below) can be checked ex-post by calling checkinvariance(Jm).")
    print("")
    qt=onestepstarHpt(pt,H)       # compute S_H(pt_0)
    if qt==0:
        rt=0
    else:
        _,rt=reduced(Poly(qt.subs(zeroind),Plist+Xind+Xpar,domain='QQ'),G0)  # rt used to extract linear constraints for membership in JR
    border=[]
    print('*** Level m=0 ***')
    if Monotony:
        print('  Frontier=',border)
    print('  tau   =1')

    coeffs = Poly(rt,Xpar).coeffs() # list of linear expressions (constraints) extracted as coefficients of rt    
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
                sigma=sigma+[{list(lin.free_symbols)[0]:s}]
    if extraptlist!=None:
        extraptlist=[Poly(seqsub(sigma,qtt),Xlist) for qtt in extraptlist]                                      # append s to list of substitutions
    print("  New linear constraint for V_0: "+str(sigma))        

    pt=Poly(seqsub(sigma,pt),Xlist)   # apply sigma to pt
    qt=Poly(seqsub(sigma,qt),Xlist)   # apply sigma to qt
    derlist=[qt] # derlist contains {S_H(pt_0),...,S_H(pt_m)}; used to generate J_m after instantiation of the parameters a_j 
    freepar = pt.free_symbols-set(Xlist) 
    zeropar = {a:0 for a in freepar}
    base = []
    updbase = False
    m=(0,)*len(Xind) # m is the last generated monomial
    initHT(len(Xind),pt)
    levelchecked=False
        
    for j in range(maxiter):
        print("")
        oldm=m
        oldeg=sum(m)
        if Monotony:
            m=nextordborder(m,border)
        else:
            m=nextordborder(m,[])
        if m==None:
            print("Frontier =",border)
            print('No new derivative below frontier: both chains stabilized')
            print('You can check invariance of Jm, that is  S_H(pt_tau)[v] in Jm for each tau and v in Vm, by calling checkinvariance(Jm)')
            print('m='+str(oldeg))
            G = groebner(Xind+base,Xind+Xpar,order='lex',domain='QQ')
            print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
            if extraptlist!=None:
                return HT[(0,)*len(Xind)],G,extraptlist 
            else:
                return HT[(0,)*len(Xind)],G
        if (not(Monotony)):
            if ((sum(m)>oldeg) & (not(Monotony))):
                if levelchecked:
                    print('Entire level checked: both chains stabilized')
                    print('m='+str(sum(m)))                    
                    G = groebner(Xind+base,Xind+Xpar,order='grevlex',domain='QQ')
                    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                    if extraptlist!=None:
                        return HT[(0,)*len(Xind)],G,extraptlist 
                    else:
                        return HT[(0,)*len(Xind)],G
                else:
                    levelchecked=True
        if sum(m)>oldeg:
            print('*** Level m=',str(sum(m))+' ***')
        if Monotony:
            print('  Frontier=',border)        
        print('  tau   =',Poly({m:1},Xind)/1)
        m,newpt=nextder(oldm,H[0],border,Monotony)  # newpt = new total derivative of pt
        qt=onestepstarHpt(newpt,H)   # qt=S_H(newpt) 
        if qt==0:
            rt=0
        else:
            _,rt=reduced(Poly(qt.subs(zeroind),Plist+Xlist[:NIND]+Xpar,domain='QQ'),G0)     # rt= S_H(newpt) mod G0
        coeffs = Poly(rt,Xpar).coeffs() # list of linear expressions (constraints) extracted as coefficients of rt      
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
                C=solve(lin,list(lin.free_symbols),solution_dict=True,rational=True)      # solve one linear constraint lin
                s = C[0]
                if type(s)==dict:
                    sigma=sigma+[s]
                else:
                    sigma=sigma+[{list(lin.free_symbols)[0]:s}]                                      # append s to list of substitutions

        if sigma==[]:
            print("  No new linear constraint detected for V_"+str(sum(m)))
            if (levelchecked | Monotony):
                print("  Checking equality with previous ideal...")
                if not(updbase):
                    zeropar = { a:0 for a in freepar }
                    base = instantiate0(derlist,zeropar,Xpar)
                    updbase = True
                newinst=instantiate0([qt],zeropar,Xpar)
                flag,Jm=check0(base,newinst,Xpar)
                if flag:
                    if Monotony:
                        print("  Equality holds, updating frontier.")
                    else:
                        print("  Equality holds.")
                    border.append(m)
                else:
                    print("  Equality does not hold, chains not yet stabilized; level not yet cleared")
                    levelchecked=False
                    base = base + newinst
            else:
                updbase=False
            derlist.append(qt)
        else:
            print("  New linear constraint for V_"+str(sum(m))+": "+str(sigma))
            updbase=False
            levelchecked=False
            derlist = [Poly(seqsub(sigma,qtt),Xlist) for qtt in derlist]
            for mon in HT.keys():
                ptm=HT[mon]
                HT[mon]=Poly(seqsub(sigma,ptm),Xlist)
            qt = Poly(seqsub(sigma,qt),Xlist)
            if qt!=0:
                derlist.append(qt)
            if extraptlist!=None:
                extraptlist=[Poly(seqsub(sigma,qtt),Xlist) for qtt in extraptlist]
            freepar = freepar-set([ list(s.keys())[0] for s in sigma]) 
    print('m='+str(sum(m)))
    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))    
    return (HT[(0,)*len(Xind)],base,extraptlist,'WARNING: maximum number of iterations reached, result is likely to be not correct.')

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
        derlist=derlist+list(var(u+dimstr))#(symbols(u+dimstr))
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
    Xlist=var(indepX)+derlist    #symbols(indepX)+derlist#var(indepX)+derlist    
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
    if type(p)!=type(Poly(0,Xlist)):
        p=Poly(p,Xlist)
    xInd=Xlist[:NIND]
    sub0={x:0 for x in xInd}
    Hi=ext2int(H)
    mons=sorted(itermonomials(xInd, n), key=monomial_key(mord, xInd))
    coefflist=[ ((onestepstarH(deltastar(p,mon2list(m,xInd),Hi[0]),Hi)/1).subs(sub0))/fact(m,xInd) for m in mons ] 
    #print(mons,coefflist)    
    return sum([c*m for c,m in zip(coefflist,mons) ])# Poly(sum([c*m for c,m in zip(coefflist,mons) ]),Xlist[:NIND])/1
    

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
    HT=OrderedDict()
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
    

def frontier(S): # computes frontier of a set of downward closed monomial; S ordered by grlex
    lmax=max([sum(m) for m in S])
    F=[]
    S=list(reversed(S))
    for m in S:
        if sum(m)==lmax:
            F.append(m)
        elif(( any(  (all(x <=y for x, y in zip(m,m1)))   for m1 in F   )  )==False):
            F.append(m)
    return F
            

### Conservation laws


def checknotzerotd(pair,t,x,Hsys,P0):
    global  Xlist
    pt=pair[0]
    px=pair[1]
    dpt=totalder(Poly(pt,Xlist),t)
    dpx=totalder(Poly(px,Xlist),x)
    if ((dpt==0)|(dpx==0)):
        return False
    if Poly(dpt+dpx,Xlist)==0:
        return False
    #print('dpt=',dpt)
    res=post(dpt,Hsys,P0,[])
    if res==False:
        return True
    #print('dpx=',dpx)
    res=post(dpx,Hsys,P0,[])
    if res==False:
        return True
    return False
 
def checknottrivial(cp,t,x,Hsys,P0):
    ntcp=[]
    tcp=[]
    for pair in cp:       
        if checknotzerotd(pair,t,x,Hsys,P0):
            ntcp.append(pair)
        else:
            if (pair[0],pair[1])!=(0,0):
                tcp.append(pair)
    return ntcp,tcp
    


def findCL(H,P0=[],base=None,ext=[],degree=2):
    global Xlist
    if base==None:
        newbase=[u0_0,u1_0,u0_1,u1_1,u2_0,u0_2]
    newbase=base+ext
    ptt,plt=genpt(newbase,degree)
    ptx,plx=genpt(newbase,degree,offset=len(plt))
    t=Xlist[0]
    x=Xlist[1]
    Dt=totalder(Poly(ptt,Xlist),t)/1
    Dx=totalder(Poly(ptx,Xlist),x)/1
    qt,J,lc=post(Dt+Dx,H,P0,plt+plx,extraptlist=[ptt,ptx])
    start_timefreepar = time.time()
    print("   Starting calculation of free parameters left...")
    Pt=lc[0]/1
    Px=lc[1]/1
    freepar= (Pt.free_symbols).union(Px.free_symbols).difference(set(Xlist))
    #sigma0={a:0 for a in freepar}
    print("   Elapsed time(freepar): %s seconds ---" % (time.time() - start_timefreepar))
    start_timecc = time.time()
    print("   Starting construction of conserved currents...")    
    J=[]
    Pt=Poly(Pt,freepar)
    Px=Poly(Px,freepar)
    zeroarg=[0]*len(freepar)
    for j in range(len(freepar)):
        sigmaj=zeroarg.copy()
        sigmaj[j]=1
        J.append((Pt(*sigmaj),Px(*sigmaj))  )
    print("   Elapsed time(cc): %s seconds ---" % (time.time() - start_timecc))
    #print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
    print("   Filtering trivial laws...")
    with suppress_stdout():
        NT,T=checknottrivial(J,t,x,H,P0)
    return NT,T

