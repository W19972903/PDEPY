# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 12:18:52 2019

@author: Miko
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 16:46:44 2019

@author: Michele
"""

# -*- coding: utf-8 -*-
# rli='import os,imp\nos.chdir(\"C:\\Users\\michele\\Dropbox\\PAPERI MIEI\\Coalgebra of PDE\")\nimport PDE\nimp.reload(PDE)\nfrom PDE import *'
# exec(rli)

#from sage.all import *
import sympy
from sympy  import *
from sympy  import itermonomials as itm
import time

 

maxiter = 60

def genpt(Xlist,deg):
    monlist=list(itm(Xlist,deg))
    l=len(monlist)
    parlist = list(var('a%d' % j) for j in range(l))
    prod=(Matrix(monlist).T *Matrix(parlist))[0,0]
    return prod, parlist

def linind(list,p):     # stub: should check linear independence of p from res
    return not(p==0)

def instantiate(ptlist, zeropar,R):
    res = []
    for pt in ptlist:
        if zeropar.keys()==[]:
            if linind(res,pt):
                res.append(pt)
        else:
            for a in zeropar.keys():
                tmp=zeropar.copy()
                tmp.update({a:1})
                p=R(pt.subs(tmp))
                if linind(res,p):
                    res.append(p)
    return res # returns list of polynomials in R


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


def check(base, newinst,R):
    J = R.ideal(base)
    for p in newinst:
        if (not(p in J)):
            return False,J
    return True,J

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
            print('q=',q)
            l,r = reduced(q,G)
            if r!=0:
                return False,G
    return True,G
    
    
step = 10

expandflag = true

def seqsub(sigma,pt):   # sequentially applies a list og substitutions in sigma to template pt
    #print('inside seqsub, sigma=', sigma)
    for s in sigma:
        #print('!!!!s=',s)
        pt=pt.subs(s)
    return pt

def expansion(p):
    if expandflag:
        return expand(p)
    else:
        return p

algebraiclosure = False

def ext2int(Huser):
    global Xlist
    return [(formeq(C[0]),C[1]) for C in Huser]

def post0(pt,Huser,P,Plist,Monotony=True):  # non sage version. Important, pt MUST be a polynomial with Xlist as generators
    start_time = time.time()
    global D, DI, NIND, HT, Xlist
    H=ext2int(Huser)#[(formeq(C[0],Xlist),C[1]) for C in Huser]
    pt=Poly(pt/1,Xlist)
    Xind=Xlist[:NIND]
    zeroind={x:0 for x in Xind}
    Xpar=parvar(H)
    G0=groebner(P,Plist+Xpar,domain='QQ',order='lex')  
    print('G0=',G0)
    qt=onestepstarHpt(pt,H)       # compute S_H(pt_0)
    if qt==0:
        rt=0
    else:
        _,rt=reduced(Poly(qt.subs(zeroind),Plist+Xpar),G0)  # rt used to extract linear constraints for membership in JR
    border=[]
    print('*** iteration j=0 ***')
    if Monotony:
        print('border=',border)
    print('tau   =1')
    print("pt0  =                   ",pt/1)
    print("qt0  = S_H(pt0) =        ",qt/1)
    print("rt0  = S_H(pt0) mod J  = ",rt/1)

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
                print("Linear constraint for V_0: "+str(sigma))
                print("Cannot find instance of given template that is an invariant.")
                print('m=0')
                print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                return False
        else:
            C=solve(lin,list(lin.free_symbols),solution_dict=True)      # solve one linear constraint lin
            #print('C=',C)
            s = C[0]
            if type(s)==dict:
                sigma=sigma+[s]
            else:
                sigma=sigma+[{list(lin.free_symbols)[0]:s}]                                      # append s to list of substitutions
            #print('sigma updated =',sigma,s)
    print("Linear constraint for V_0: "+str(sigma))        

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
        print('*** iteration j=',str(j+1)+' ***')
        if Monotony:
            print('border=',border)
        else:
            oldeg=sum(m)
        m,newpt=nextder(m,H[0],border,Monotony)  # newpt = new total derivative of pt
        if m==None:
            print('No new derivative extra border: both chains stabilized')
            print('m='+str(j+1))
            print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))                
            return (HT[(0,)*len(Xind)],Jm)
        if (not(Monotony)):
            if ((sum(m)>oldeg) & (not(Monotony))):
                if levelchecked:
                    print('Entire level checked: both chains stabilized')
                    print('m='+str(j+1))
                    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                    return (HT[(0,)*len(Xind)],Jm)
                else:
                    levelchecked=True
        print('tau   =',Poly({m:1},Xind)/1)
        print("pt"+str(j+1)+"= ",newpt/1) 
        qt=onestepstarHpt(newpt,H)   # qt=S_H(newpt) 
        print("qt"+str(j+1)+"= ",qt/1) 
        if qt==0:
            rt=0
        else:
            _,rt=reduced(Poly(qt.subs(zeroind),Plist+Xpar),G0)     # rt= S_H(newpt) mod G0
        print("rt"+str(j+1)+"= ",rt/1)
        coeffs = Poly(rt/1,Xpar).coeffs() # list of linear expressions (constraints) extracted as coefficients of rt    
        sigma=[]                                                
        for lin in coeffs:      # solve the list of linear constraints, one by one
            lin=seqsub(sigma,lin)  # apply previously accumulated subst. to current constraint lin
            lin=Add(lin,0)           # forces lin to symbolic expression
            #print('lin=',lin)
            #print('lin free symbols=',lin.free_symbols)
            if lin.free_symbols==set():
                if lin==0:
                    C=True
                else:
                    sigma=sigma+[{lin:0}]
                    print("Linear constraint for V_"+str(j+1)+": "+str(sigma))
                    print("Cannot find instance of given template that is an invariant.")
                    print('m='+str(j+1))
                    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                    return False
            else:
                C=solve(lin,list(lin.free_symbols),solution_dict=True)      # solve one linear constraint lin
                #print('C=',C)
                s = C[0]
                if type(s)==dict:
                    sigma=sigma+[s]
                else:
                    sigma=sigma+[{list(lin.free_symbols)[0]:s}]                                      # append s to list of substitutions
                #print('sigma updated =',sigma,s)
        print("Linear constraint for V_"+str(j+1)+": "+str(sigma))

        if sigma==[]:
            print("Vector space equality detected: V_"+str(j)+" = V_"+str(j+1))
            print("Checking ideal equality J_"+str(j)+" = J_"+str(j+1)+"...")
            if not(updbase):
                zeropar = { a:0 for a in freepar }
                base = instantiate0(derlist,zeropar,Xpar)
                updbase = True
            newinst=instantiate0([qt],zeropar,Xpar)
            flag,Jm=check0(base,newinst,Xpar)
            if flag:
                print("Equality holds, updating border")
                border.append(m)
            else:
                print("Equality does not hold, chains not yet stabilized; level not yet checked")
                derlist.append(qt)
                base = base + newinst
                levelchecked=False
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
            #print('sigma=',sigma)
            freepar = freepar-set([ list(s.keys())[0] for s in sigma]) 
            #print('freepar=',freepar)
    print('m='+str(j))
    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))    
    return (HT[(0,)*len(Xind)],base,'WARNING: maximum number of iterations reached, result is likely to be not correct.')

    
    


    
def stringify(L):
    s = "".join(str(x)+',' for x in L)
    return  s[0:len(s)-1]


MAXITER = 50        



##### PDE

# generate dependent var, their derivatives and the corresponding dictionary and complete var list

#from collections import Counter 
import numpy as np
from sympy.polys.orderings import monomial_key

global D
global DI
global NIND
global Xlist
#import copy

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
        #dimstr=('(:'+str(n)+')')*NIND
        #print(dimstr)
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
    #D={NIND+j: ([int(i) for i in  (str(derlist[j])[1:])])+[str(derlist[j])[0]] for j in range(len(derlist))}
    #DI={tuple([int(i) for i in  (str(derlist[j])[1:])]+[str(derlist[j])[0]]):NIND+j for j in range(len(derlist))}
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
                #print('***dxiplust=',dxiplust)
                l=len(dxiplus)
                xnew=var(''.join(list([dxiplus[l-1]]+['_'+str(c) for c in dxiplus[:l-1] ]))  )
                Xlist.append(xnew)
                #print('***xnew=',xnew)
                lx=len(Xlist)#-NIND
                D[lx-1]=dxiplus
                DI[dxiplust]=lx-1
                mnew=mnew+[0]
                changexl=True
            h=DI[dxiplust]
            #print('mnew,h=',mnew,h)
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
    #if not(type(p)==sympy.polys.polytools.Poly):
    #p=Poly(p/1,Xlist)
    pdic=p.as_dict()
    LC=[]#Counter({})
    chxl=False
    for m in pdic.keys():
        v=pdic[m]
        dmon,changexl=dermon(m,x)
        #print('m,changexl=',m,changexl)
        LC.append({k:dmon[k]*v for k in dmon.keys()})#C+Counter({k:dmon[k]*v for k in dmon.keys()})
        if changexl:
            chxl=True
    sd=sumdict(LC,chxl)
    #print('sd=',sd)
    return Poly(sd,Xlist) 

def totalderstar(p,tau):
    global Xlist
    #if not(type(p)==sympy.polys.polytools.Poly):
    q=Poly(p/1,Xlist)
    #else:
    #    q=p
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
    #p=Poly(p/1,Xlist)
    mons=p.as_dict()
    for m in mons:
        for i in [j for j in np.nonzero(m)[0] if j>=NIND]:
            search=searchprinc(i,Comp)
            if search!=None:
                tau= [item for sublist in [[Xlist[j]]*search[0][j] for j in range(NIND)] for item in sublist]
                #print(tau, search[1])                
                G=totalderstar(search[1],tau)/1
                q=p/1
                return Poly(q.subs({Xlist[i]:G}),Xlist)#q.subs({Xlist[i]:G})#Poly(q.subs({Xlist[i]:G}),Xlist)
    return None

            
def onestepstar(p,C): # important, p MUST be a polynomial with Xlist as generators
    #p=Poly(p/1,Xlist)
    global Xlist    
    while (True):
        q=onestep(p,C)
        if q==None:
            return p
        p=q
                
def onestepH(p,H): # important, p MUST be a polynomial with Xlist as generators
    global NIND, Xlist
    #q=p/1
    #sigma0={ x:0 for x in Xlist[:NIND]}
    #q1=q.subs(sigma0)
    #if q1!=q:
    #    return Poly(q1,Xlist)
    for C in H:
        q=onestepstar(p,C)
        if q!=p:
            return q    
    return None
    
def onestepstarH(p,H):
    global Xlist
    #p=Poly(p/1,Xlist)
    while (True):
        q=onestepH(p,H)
        if q==None:
            return p
        p=q


def onestepstarHpt(pt,H):
    global Xlist
    #pt=Poly(pt/1,Xlist)
    ptdic=pt.as_dict()
    newpt=0
    for m in ptdic:
        #print('m=',m)
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
    #print('x,pt=',x,pt/1)
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
    


"""
EXPERIMENTS  Feb 2, 2019. Machine: Microsoft Surface Pro 4
@author: Michele

Everywhere, ordering of independent variables: t>x, u=u(t,x)

1) HEAT EQUATION
u_t = b*u_xx
Boundary and initial conditions: u(0,x)= sin(c*x),  u(t,0)=exp(-b*c^2*t), for generic b,c.

xheat=initvar(['t','x'],['b','c','f','g','h','u'],maxder=1,extrader={'u':':2_(:3)'})
eq1=[(u1_0,b00*u0_2),(b10,0),(b01,0),(c10,0),(c01,0),(f10,0),(g10,0),(h01,0)]
eq2=[(u0_1,g01),(f01,-c00*g00), (g01,c00*f00)]
eq3=[(h10,-b00*c00**2*h00)]

C1h=(eq1,[t,x])
C2h=(eq2,[x])
C3h=(eq3,[t])
Hheat=[C1h,C2h,C3h]
P0=[u0_0,g00,f00-1,h00-1]
pt=Poly(a0+a1*u0_0+a2*g00*h00+a3*f00*h00+a4*f00**2+a5*g00**2,xheat)

resh=post0(pt,Hheat,P0,[a0,a1,a2,a3,a4,a5])

Output:
(...)
*** iteration j= 10 ***
Entire level checked: both chains stabilized
m=10
--- Elapsed time: 4.9931557178497314 seconds ---

res[0]/1
Out[117]: -a1*g00*h00 + a1*u0_0 + a5*f00**2 + a5*g00**2 - a5
a1=1,a5=0: separation of variable solution of heat equation; 
a1=0,a5=1: trigonometric identity (cos^2+sin^2=1)


2) INVISCID BURGER EQUATION
u_t = u*u_x
Boundary condition: u(0,x)=c*x+b, for generic c,b

xburger=initvar(['t','x'],['c','b','u'],maxder=1)
eq1b=[(u10,-u00*u01),(c10,0),(c01,0),(b10,0),(b01,0)]
eq2b=[(u01,c00)]

C1b=(eq1b,[t,x])
C2b=(eq2b,[x])
Hb=[C1b,C2b]
P0b=[u00-b00]
pvb=[c00, b00, u00]
pt,pl=genpt([t,x]+pvb,3)

resb=post0(pt,Hb,P0b,pl,False)

Output:
(...)
*** iteration j= 15 ***
Entire level checked: both chains stabilized
m=15
--- Elapsed time: 14.150203227996826 seconds ---

Poly(res[0]/1,[u00])/1
Out[314]: a50*b00 + a50*c00*x + u00*(-a50*c00*t - a50)
a50=1: complete integral of B.E., u(t,x)=(c*x+b)/(c*t+1)


3) KdV
u_t+u*u_x+δ*u_xxx = 0
To do.

4) mKdv
u_t+u^2*u_x+δ*u_xxx = 0
To do.


PRECONDITIONS
Introduce new 'control variables' (constant) i,j,k. 

xheat=initvar(['t','x'],['b','c','d','f','g','h','j','i','k','u'],maxder=1,extrader={ 'u':':2_(:3)'})

Express 'generic' initial and boundary conditions: 
u(0,x)= sin(c*x),  u(t,0)=exp(-i*t), for generic c,i.

Heat= [([(u1_0, b00*u0_2), (b10, 0),
   (b01, 0),
   (c10, 0),
   (c01, 0),
   (f10, 0),
   (g10, 0),
   (h01, 0),
   (i10, 0),
   (i01, 0),
   (j10, 0),
   (j01, 0),
   (k10, 0),
   (k01, 0),
   (d01, 0),
   (d10, 0)], 
   [t, x]), ([(u0_1, g01), (f01, -c00*g00), (g01, c00*f00)], [x]),
 ([(h10, -h00*i00)], [t])]

P0=[b00,c00,i00,u0_0,h00-1,j00,k00]  
pt=a0*(u0_0+j00*g00*h00+k00*f00*h00) 

P0 ensures pt is trivially an invariant for any a0. Hence the returned Groebner basis encodes 
the *weakest* precondition (on variables in xheat) such that u0_0+j00*g00*h00+k00*f00*h00 is an invariant.
Indeed:
    
res=post0(pt,Heat,P0,[a0],True)

Outcome:
(...)
*** iteration j= 10 ***
Entire level checked: both chains stabilized
m=10
--- Elapsed time: 7.865466594696045 seconds ---

res[1]
Out[163]: GroebnerBasis([c00**2*h00*j00**2*u0_0 + c00**2*h00*k00**2*u0_0 - c00**2*f00*k00 + c00**2*j00*u0_0, c00*h00*i00*j00**2*u0_0 + c00*h00*i00*k00**2*u0_0 - c00*f00*i00*k00 + c00*i00*j00*u0_0, h00**2*i00*j00**2*u0_0 + h00**2*i00*k00**2*u0_0 + 2*h00*i00*j00*u0_0 + i00*u0_0, h00*i00*j00**2*u0_0**2 + h00*i00*k00**2*u0_0**2 - f00*i00*k00*u0_0 + i00*j00*u0_0**2, b00*c00**3*f00 - c00*f00*i00, b00*c00**2*f00*k00 - h00*i00*j00**2*u0_0 - h00*i00*k00**2*u0_0 - i00*j00*u0_0, c00**2*f00**2*k00 + c00**2*k00*u0_0**2, c00**2*f00*h00*k00 + c00**2*h00*j00*u0_0 + c00**2*u0_0, c00*f00**2*h00*k00 + c00*g00**2*h00*k00 - c00*f00*g00 + c00*f00*u0_0, c00*f00**2*i00*k00 + c00*i00*k00*u0_0**2, c00*f00*h00*i00*k00 + c00*h00*i00*j00*u0_0 + c00*i00*u0_0, c00*g00*h00*i00*k00 - c00*h00*i00*k00*u0_0, f00*h00*i00*j00*u0_0 - h00*i00*k00*u0_0**2 + f00*i00*u0_0, f00**2*i00*k00*u0_0 + i00*k00*u0_0**3, f00*h00*i00*k00*u0_0 + h00*i00*j00*u0_0**2 + i00*u0_0**2, c00*f00*g00*j00 + c00*f00**2*k00 - c00*f00*j00*u0_0 + c00*g00*k00*u0_0, c00*f00*h00*j00 - c00*g00*h00*k00 + c00*f00, c00*f00*g00*i00 - c00*f00*i00*u0_0, b00*c00**2*u0_0 - i00*u0_0, c00**2*g00 - c00**2*u0_0, g00*h00*j00 + f00*h00*k00 + u0_0, g00*i00*u0_0 - i00*u0_0**2], t, x, b00, c00, f00, g00, h00, j00, i00, k00, u0_0, domain='ZZ', order='grevlex')

E=list(res[1])
E0=[p.subs({g00:0,f00:1,h00:1}) for p in E]  # set known initial values for cos, sin, exp

# solve for the control parameters i, j, k and initial condition u(0,0) 
solve(E0+[c00*f00-1], [k00,j00,i00,u0_0],dict=True) # use f00 as a dummy variable to encode c!=0
Out[173]: [{u0_0: 0, k00: 0, i00: b00*c00**2, j00: -1}]
We see: j=-1, k=0, i=b*c^2, u(0,0)=0, as expected.


TRIGONOMETRIC
xtri=initvar(['x','y'],['f','g','i','j','h','k'],maxder=1,extrader={'f':':2_(:2)','g':':2_(:2)','i':':2_(:2)','j':':2_(:2)','h':':2_(:2)','k':':2_(:2)'})
eqtr=[(f1_0,-g0_0), (f0_1,-g0_0), (g1_0,f0_0), (g0_1,f0_0), (i1_0,-j0_0), (i0_1,0), (j1_0,i0_0), (j0_1,0),  (h1_0,0), (h0_1,-k0_0), (k1_0,0), (k0_1,h0_0) ]
Htr=[(eqtr,[x,y])]
pt,pl=genpt([f0_0, g0_0, i0_0, j0_0, h0_0, k0_0],2) # f=cos(x+y), g=sin(x+y), i=cos(x), j=sin(x), h=cos(y), k=sin(y)
P0=[f0_0-1, g0_0, i0_0-1, j0_0, h0_0-1, k0_0]
restr=post0(pt,Htr,P0,pl,False)
""" 


