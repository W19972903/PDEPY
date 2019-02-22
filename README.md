# PDEPY
Python code for the examples in the paper "Coalgebra, partial differential equations and boundary problems"
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

