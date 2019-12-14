#----------------------------------------------------
# Author: Lucas Marques Moreno
# date: Apr/2018 
# Script used in the SAE AeroDesign 2018 Competiton
# This script is responsble for the Performance study
#----------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from CD import Cd
from propeller import propeller



Cruise = True  
takeoff = True 
V_n = True 
CD0 = False # disable the Roskam algorithm
plotCD = True
Pont = True

#============== Input Data ================

#----------- Nature Constants --------------------

g = 9.81#[m/s²] gravity
u = 1.7894*10**(-5)#[N.s/m²] air dynamic viscosity 
p = 1.2 #[kg/m³] air mass specific

#-----------Aircraft Input ------------------
m_a = 2.5#[kg] aircraft empty mass 
cp = 6.#[kg] mass of the load
m = m_a+cp
w = m*g #[N] aircraft gross weight
CDo = 0.022#[]parasite drag coefficient
k = 0.0535#[]constant for induced drag coefficient
CL_max = 1.5# maximum lift coefficient

'''
Propellers Aviable :
    12X4 - 12X5
    13X4 - 13X5
    14X4 - 14X5
'''
hel = '13X4'

#///////////// Wing ////////////

S = 0.192#[m²] wing area
b = 1.2#[m] wing span
c_asa = 0.16#[m] wing chord
h_asa = 0.16929#[m] distance between wing and ground
airfoil_asa = 'e61'#[] name of the airfoil (without extension)

#////////// Fuselage ///////

L_fus = 0.398#[m] fuselage length
D_fus = 0.150#[m] fuselage diameter
A_fus = np.pi*D_fus*L_fus#[m²] superficial area 

#///////// Landing gear ////////

area_lg = (90*10**(-3))*(9*10**(-3))*3 #[m²] frotnal area

#////////////// Empennage ////////////

#_____________ Horizontal Stabilizer _________

airfoil_estH= 'n0009sm'# horizontal stabilizer airfoil name (without extenson)
c_estH = 0.091#[m] horizontal stabilizer chord
b_estH = 0.32#[m] horizontal stabilizer span
S_H = c_estH * b_estH  #[m²] horizontal stabilizer area


#_____________ Vertical Stabilizer _________

airfoil_estV= 'n0009sm'# vertical stabilizer airfoil name (without extenson)
c_estV = 0.091#[m] vertical stabilizer chord
b_estV = 0.191#[m] vertical stabilizerspan
S_V = c_estV * b_estV#[m²] vertical stabilizer area  

#--------------- Velocity Range -------------
Vstall = (w/(0.5*p*S*CL_max))**0.5
Vi = round(Vstall*0.9)#[m/s] initial velocity
Vf =40#[m/s] final velocity
precisao = 95#[%] velocity definition precison


#=========== Dados distancia de takeoff ==============

ua = 0.05#[] friction coefficient
CLto = CL_max/1.44#[] Hoerner Lift
seg = 20.#[%] safe margin in the stall velocity for takeoff 


#===========  Diagram V-n Data ========================

Clmax = CL_max# wing profile maximum lift coefficient 
pvn = 0.1#diagram V-n step calculation
v_Max = 26#[m/s] maximum flight velocity
c_barra = c_asa #[m] aerodynimic chord
n_max = 2.# maximum load factor
v_c= 18.#[m/s] maximum cruise velocity
c1=0.25#Ratio between positives and negatives loads
c2=1.5#safety coefficent
c3 = 1.25#Ratio Vd/Vmax
n_lim = 2.5#Limit for Load factor
Cl_vn = 1.427 #tangent of CL/alfa curve
Alfa_vn = 10#[º] angle of attack



#------------- Determinar Velocidade para o diagrama -----------------

U_g1=8#[m/s] Vertical gust velocity
U_g2=10#[m/s] Vertical gust velocity
U_g3=13#[m/s] Vertical gust velocity

#================== Structure Input ============================
D_asa = 23.4#[mm] wing spar diameter
Esp_asa = 1.#[mm] wing spar thickness 
d_apoio = 100.#[mm] distance between the support points

D_empH = 5.#[mm] horizontal stabilizer spar external diameter
Esp_empH = 1.#[mm] horizontal stabilizer spar thickness 
d_apoioH = 5.#[mm] distance between the support points

D_empV = 5.#[mm] vertical stabilizer spar external diameter
Esp_empV = 1.#[mm] vertical stabilizer spar thickness 

safe = 1.2#[] safety coefficient 
sigmaF = 180.#[MPa] maximum yield tension
sigmaC = 80.#[MPa] maximum shear tension
pmt = 1.6e3#[kg/m3] wing spar specific mass 

if Pont:
    EE = cp/m_a
    Ponto = 15*EE+cp



#------------ Engine Data ----------------------

H = propeller(hel)
A = H.A
B = H.B
C = H.C
D = H.D
rho = H.rho

#----- Drag function input vector definition  ----------
r =[p,u,S,b,area_lg]
profile = [airfoil_asa,airfoil_estH,airfoil_estV]
Asa = [b,c_asa,S,h_asa,d_apoio]
fus = [L_fus,D_fus,A_fus]
Emp_H = [b_estH,c_estH,S_H,0,d_apoioH]
Emp_V = [b_estV,c_estV,S_V,0]
margem_CD = 0#[%] additonal drag for correction


def frange(start, stop, step):
     i = start
     while i < stop:
         yield i
         i += step

################### Thrust Performance #####################
passo = 0.1  
if Cruise:
    td=[]#[N] engine thurst
    pd=[]#[W] engine power
    vd=[]#[m/s] airfram velocity
    tr=[]#[N] requeried thrust
    pr = []#[W] requeried power
    cl=[]#[] thoerical lift coefficient
    cd=[]#[] thoerical drag coefficient
    clxcd=[]#[] lift/drag ratio
    rc = [] # [m/s] climb ratio
    const = Vi  
    l=[]
    d = []
    va = []
    c_rho = p/rho
    lmi = .9*precisao/100
    lma = 1.1*precisao/100
    while const<Vf: 
        V = const
        Td = (A*V**3+B*V**2+C*V+D)*c_rho # Engine thrust equation 
        td.append(Td)
        Pd = Td*V
        pd.append(Pd)
        vd.append(V)
        CL = (2*w)/(p*(V**2)*S)# requeried lift coefficient
        if CD0:
            CD = CDo + k*CL**2# total drag coefficent
        else:
            CD = Cd(r,profile,V,CL,Asa,Emp_H,Emp_V,area_lg,fus,margem_CD,'Voo',plotCD)
        Tr = (1/2)*p*(V**2)*S*CD[0]# requeried thrust
        tr.append(Tr)
        Pr = Tr*V
        pr.append(Pr)
        rv = round((Td/Tr),1)
        if lmi<rv<lma:
            va.append(V)
        cl.append(CL)
        cd.append(CD)
        ClxCd = CL/CD[0] 
        clxcd.append(ClxCd)  
        RC = (Td*V-Tr*V)/w 
        if RC<=0:
            rc.append(0)
        else:
            rc.append(RC)
        
        const = const+passo
################## Takeoff Distance ############################
if takeoff:
    
    Vto = Vstall*(1+seg/100)
    if CD0:
        phi = ((16*h_asa/b)**2)/(1+(16*h_asa/b)**2) # Ground Effect
        CDto = CDo+ phi*k*CLto**2
    else:
        CDtO =Cd(r,profile,Vto,CLto,Asa,Emp_H,Emp_V,area_lg,fus,margem_CD,'Tk',plotCD)
        CDto = CDtO[0]
    thRC =np.rad2deg(np.arctan(1/(CLto/CDto)))
    CLxCD = CLto/CDto
    if Cruise:
        RCmax = max(rc)
    def integrand(Vs,g,A,B,C,D,w,ua,p,CDto,CLto,S):
        aux1 = (A*Vs**3+B*Vs**2+C*Vs+D)*c_rho/w
        dyn_p = 0.5*p*(Vs**2)
        aux2 = ua-(dyn_p*(CDto-ua*CLto)/(w/S))
        return Vs/(g*(aux1+aux2))
    St=quad(integrand,0,Vto, args=(g,A,B,C,D,w,ua,p,CDto,CLto,S))
    print('\n---------------- Takeoff Distance --------------------')

    print('\nTakeoff Distance: ',round(St[0],2),'m')  

    print('\nMaximum Climb Angle: ', round(thRC,2),'°')
    
    print('\nCL/CD Ratio ',round(CLxCD,2))
    if Cruise:
        print('\nMaximum Climb Ratio: ', round(RCmax,2), 'm/s')
    
################## Diagram V-n ################################
if V_n:

    v_max = v_Max
    
    CL_alpha_inc=np.arctan(Alfa_vn/Cl_vn) 
    CL_alpha_D=(180/np.pi)*CL_alpha_inc
    v_d=c3*v_max#[m/s]diving velocity
    v_estol=Vstall 
    v_manobra=v_estol*n_max**0.5#[m/s] maneuver velocity
    n_lim_pos=n_max
    n_lim_neg=-c1*n_lim_pos
    n_ult_pos=c2*n_lim_pos
    n_ult_neg=c2*n_lim_neg
    mi=2*w/(S*p*g*c_barra*CL_alpha_inc)  
    K_g=0.88*mi/(5.3+mi)
    vv=((-2*n_lim_neg*w)/(p*S*Clmax))**0.5
    vv2=((2*n_ult_pos*w)/(p*S*Clmax))**0.5
    vv3=((-2*n_ult_neg*w)/(p*S*Clmax))**0.5 
    v_lim_pos = []
    v_lim_neg = []
    n_pos=[]
    n_neg = []
    n_neg2 = []
    v_ult_pos = []
    n_pos2 =[]
    v_ult_neg = []
    v_raj = []
    n_raj11 = []
    n_raj12 = []
    n_raj13 = []
    n_raj21 = []
    n_raj22 = []
    n_raj23 = []
    
    for vlimpos in frange(0,v_manobra+pvn,pvn):
        v_lim_pos.append(vlimpos)
        aux1 = S*Clmax/w  
        aux2 = 0.5*p*vlimpos**2
        npos = aux2*aux1
        n_pos.append(npos)
    
    for vlimneg in frange(0,vv+pvn,pvn):
        v_lim_neg.append(vlimneg)
        nneg=-0.5*p*vlimneg**2*S*Clmax/w
        n_neg.append(nneg)
    
    for vultpos in frange(v_manobra,vv2+pvn,pvn):
        v_ult_pos.append(vultpos)
        npos2=0.5*p*vultpos**2*S*Clmax/w
        n_pos2.append(npos2)
        
    for vultneg in frange(vv,vv3+pvn,pvn):
        v_ult_neg.append(vultneg)
        nneg2=-0.5*p*vultneg**2*S*Clmax/w
        n_neg2.append(nneg2)
        
    for vraj in frange (0,v_d,pvn):
        v_raj.append(vraj)
        nraj11=1+((p*vraj*CL_alpha_inc*K_g*U_g1*S)/(2*w))
        n_raj11.append(nraj11)
        
        nraj21=1-((p*vraj*CL_alpha_inc*K_g*U_g1*S)/(2*w))
        n_raj21.append(nraj21)
        
        nraj12=1+((p*vraj*CL_alpha_inc*K_g*U_g2*S)/(2*w))
        n_raj12.append(nraj12)
        
        nraj22=1-((p*vraj*CL_alpha_inc*K_g*U_g2*S)/(2*w))
        n_raj22.append(nraj22)
        
        nraj13=1+((p*vraj*CL_alpha_inc*K_g*U_g3*S)/(2*w))
        n_raj13.append(nraj13)
        
        nraj23=1-((p*vraj*CL_alpha_inc*K_g*U_g3*S)/(2*w))
        n_raj23.append(nraj23)

    
    x1=[]
    x2=[]
    x3=[]
    x4=[]
    x5=[]
    x6=[]
    x7=[]
    x8=[]
    x9=[]
    x10=[]    
    y1=[]
    y2=[]
    y3=[]
    y4=[]
    y5=[]
    y6=[]
    y7=[]
    y8=[]
    y9=[]
    y10=[]
    
    for X1 in frange(v_manobra,v_d+pvn,pvn):
        y1.append(n_lim_pos)
        x1.append(X1)        
        
    for X2 in frange(vv,v_d+pvn,pvn):
        x2.append(X2)        
        y2.append(n_lim_neg)
        
    x3=[v_d,v_d]
    y3=[n_lim_pos,n_lim_neg]
    
    for X4 in frange(vv2,v_d+pvn,pvn):
        x4.append(X4)
        y4.append(n_ult_pos)
        
    for X5 in frange(vv3,v_d+pvn,pvn):
        x5.append(X5)
        y5.append(n_ult_neg)
        
    x6=[v_d,v_d]
    y6=[n_lim_pos,n_ult_pos]
    
    x7=[v_d,v_d]
    y7=[n_lim_neg,n_ult_pos]
    
    x8=[v_estol,v_estol]
    y8=[0,1]
    
    
    x9=[v_manobra,v_manobra]
    y9=[0,n_lim_pos]
    
    
    x10 = [v_c,v_c]
    y10 = [0,n_lim_pos]

if Pont:
    print('\n----------------------- Score -----------------')    
    print('\nEmpty Mass: ',m_a,'kg')
    print('\nLoad: ',cp,'kg')
    print('\nTotal Mass ',m,'kg')
    print('\nStrucutral Efficiency: ',round(EE,2))
    print('\nScore: ', round(Ponto,2),'Points')

    

#=========== Stall Velocity Limit ===============
if Cruise:
    erro = 100-precisao
    xv = [Vstall,Vstall]
    yv=[]
    if min(td)<=min(tr):
        yv.append(min(td))
    else:
         yv.append(min(tr))  
             
    if max(td)>=max(tr):
        yv.append(max(td))
    else:
         yv.append(max(tr))   
         
if Cruise:
    erro = 100-precisao
    xp = [Vstall,Vstall]
    yp=[]
    if min(pd)<=min(pr):
        yp.append(min(pd))
    else:
         yp.append(min(pr))  
             
    if max(pd)>=max(pr):
        yp.append(max(pd))
    else:
         yp.append(max(pr))
         
#====================================================
    print('---------------- Cruise Flight --------------------\n')
    if len(va)>0:
        if va[0]<Vstall:
            print('Minimum Velocity:', round(va[0],2),'m/s error +- ',round(erro,2),'%')
        print('Maximum Velocity: ', round(va[len(va)-1],2),'m/s error +- ',round(erro,2),'%\n')
    if V_n:
        print('Diving Velocity: ',round(v_d,2),'m/s')
    plt.figure(1)
    plt.plot(vd,td)
    plt.ylabel('Engine Thrust [N]')
    plt.xlabel('Velocity [m/s]')
    plt.title('Engine Thrust X Velocity')
    plt.grid() 
    
    
    plt.figure(2)
    plt.plot(vd,tr)
    plt.ylabel('Requeried Thrust[N]')
    plt.xlabel('Velocity [m/s]')
    plt.title('Rqueried Thrust X Velocity')
    plt.grid()
   
    
    plt.figure(3)
    plt.plot(vd,td,label = "Engine Thrust")
    plt.plot(vd,tr,label = "Rqueried Thurst")
    plt.plot(xv,yv, label = "Stall Velocity")
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Thtust [N]')
    plt.title('Thrust Performance')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.grid()

    
    plt.figure(4)
    plt.plot(vd,pd,label = "Engine Power")
    plt.plot(vd,pr,label = "Requeried Power")
    plt.plot(xp,yp, label = "Stall Velocity")
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Power [W]')
    plt.title('Power Performance')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.grid()
 
    
    plt.figure(5)
    plt.plot(vd,clxcd)
    plt.xlabel('Performance [m/s]')
    plt.ylabel('CL/CD')
    plt.title('CL/CD Performance')
    plt.grid()

    
    plt.figure(6)
    plt.plot(vd,rc)
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Climb Ratio')
    plt.title('Climb Performance')
    plt.grid()


if V_n:
    plt.figure(7)
    plt.plot(v_lim_pos,n_pos)
    plt.plot(x1,y1,'black')
    plt.plot(x2,y2,'black')
    plt.plot(x3,y3,'black')
    plt.plot(x4,y4,'red')
    plt.plot(x5,y5,'red')
    plt.plot(x6,y6,'red') 
    plt.plot(x7,y7,'red') 
    plt.plot(x8,y8,'yellow', label = "Stall Velocity")
    plt.plot(x9,y9,'navy', label ="Maneuver Velocty")
    plt.plot(x10,y10,'magenta', label ="Cruise Velocity") 
    plt.plot(v_lim_neg,n_neg,'k')
    plt.plot(v_ult_pos,n_pos2,'red')
    plt.plot(v_ult_neg,n_neg2,'red')
    plt.plot(v_raj,n_raj11,'g--',label= "Gust Velocity 1")
    plt.plot(v_raj,n_raj21,'g--',label= "Gust Velocity 1")
    plt.plot(v_raj,n_raj12,'c--',label= "Gust Velocity 2")
    plt.plot(v_raj,n_raj22,'c--',label= "Gust Velocity 2")
    plt.plot(v_raj,n_raj13,'b--',label= "Gust Velocity 3")
    plt.plot(v_raj,n_raj23,'b--',label= "Gust Velocity 3")
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Load Factor')
    plt.title('Diagram V-n')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.grid()

plt.show()
