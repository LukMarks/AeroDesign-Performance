import numpy as np
import os


canard = False




def Cd(r,profile,V,Cl,Asa,Emp_H,Emp_V,la,fus,s,Voo,plotCD):
    
    def xt(airfoil):
        x = []
        y = []
        yu = []
        xu = []
        yl =[]
        xl = []
        prof_name = airfoil+'.dat'
        f=open(prof_name,'r')
        #check the order of the points
        Tks = 0
        for line in f:
            X = float(line.strip().split()[0])
            Y = float(line.strip().split()[1])
            x.append(X)
            y.append(Y)



        for i in range(1,len(x)):
            if x[i]<x[i-1]:
                cp = i+1

        for i in range(0,cp):
            xu.append(x[i])
            yu.append(y[i])

        y.reverse()
        for i in range(cp, len(x)):
            xl.append(x[i])
            yl.append(y[i])

        for i in range(0,len(yl)):
                for j in range(0,len(yu)):
                    tks = abs(yu[j]-yl[i])
                    if Tks<tks:
                        Tks=tks
                        ref1 = y[i]
                        ref2 = y[j]
                        x_t = xl[i]
                       
        '''
        for g in range (0,len(x)):
            if y[g] == ref1:
                x_t = xl[g]
            elif y[g] == ref2:
                x_t = x[g]
        '''
        

        return Tks,x_t

    a_asa = xt(profile[0])
    a_estH = xt(profile[1])
    a_estV = xt(profile[2])



        #---------- CD_wing(pag. 23 ~ 28) --------------
    def CD_wing(r,V,Cl,a,Asa,Nome,Voo):
        #Rwf = 1.08 #[] fator de indução fig 5.11
        	#Encontrado por meio do R_Nf = p*U_l*c_we/u, o número de Reynolds da asa
        p = r[0]#1.09 #[kg/m^3] air mass specific
        u = r[1]#1.78*10**(-5) # air dynamic viscosity 
        U_l = V#14. #[m/s] U_l = stall valocity
        S = Asa[2]#0.365714 #[m] wing area
        b = Asa[0]#1.6 #wing span
        c_we = Asa[1]#S/b #c_we = wet wing mean geometric chord
        V_som = 343#[m/s] sound speed
        M = V/V_som
        Sref = r[2]
        R_Nw = p*U_l*c_we/u 

        Rls = 1.07 

        if Nome == 'asa':
                Rwf = 1.08
        else:
                Rwf = 1.03
        
        
        aux3 = (np.log10(R_Nw))**2.58
        aux4 = (1+0.144*M**2)**0.58
        cfw	= 0.455/(aux3*aux4)

       # x_t #maximum thickness point
        x_t = a[1]
        if x_t>= 0.3:
        	L_l = 1.2         
        else:
        	L_l = 2.0

        tc = a[0]
        if tc>=0.05:
            K_w = 1.9767 + 0.5333*tc
        else:
            K_w = 2.0

        S_exp = S #[m²] S_exposed 
        Swet_w = K_w*S_exp #[m²] wet area
        
        aux5 = L_l*tc
        aux6 = 100*tc**4
        CD_0w = Rwf*Rls*cfw*(1+aux5+aux6)*(Swet_w/Sref)


        C_L =Cl #Lift coefficient
        C_Lw = 1.0*C_L 
        b = Asa[0]
        S = Asa[2]
        
        AR = (b**2)/S 
        e = 0.85 

        '''
        #--- Teste do e(fator de eficiência de Oswald/da envergadura) segundo valores do livro(part IV pág.27 ---#
        C_la= #procurar na Ref.:HOak, D.E., USAF Stability and Control Datcom, Flight Control Division, Air Force Flight Dynamics Laboratory
        a_semic= 0#ângulo de enflexamento no meio da corda
        Beta = (1-M**2)**(1/2)
        k = (C_la)/(2*np.pi/Beta)
        C_Law = 2*np.pi*A/((2+(A**2*Beta**2/k**2)*(1+(np.tan(v_semic))**2/Beta**2)+4)**(1/2))
        e = 1.1*(C_Law/A)/(R*(C_Law/A)+(1-R)*np.pi) #Fator de eficiência do Oswald
        R
        '''
        CD_Lw =(C_Lw**2)/(np.pi*AR*e)
        if Voo == 'Tk':
            h_asa = Asa[3]
            if Nome == 'asa':
                phi = ((16*h_asa/b)**2)/(1+(16*h_asa/b)**2)
                CD_wing = CD_0w + phi*CD_Lw
                
            else:
                CD_wing = CD_0w 
                
        else:
                
            if Nome == 'asa':
                CD_wing = CD_0w + CD_Lw
                #print('\ncomeca aqui')
                #print(CD_0w, CD_Lw,C_Lw)
                #print('asa',CD_wing)
            else:
                CD_wing = CD_0w
                #print('emp',CD_wing)
        
        return CD_wing
    nome = 'asa'
    nome2 = 'emp'
    wing = CD_wing(r,V,Cl,a_asa,Asa,nome,Voo)
    EH = CD_wing(r,V,Cl,a_estH,Emp_H,nome2,Voo)
    EV = CD_wing(r,V,Cl,a_estV,Emp_V,nome2,Voo)
    empennage = EH+EV
  
    #---------- CD_fuselage (pág. 44 ~ 47)----------

    def CD_fus(r,V,fus,Asa):

        p = r[0]#[kg/m^3] air mass specific
        u = r[1]#air dynamic viscosity 
        U_l = V#[m/s] flight velocity
        S = fus[2]#0.365714 #[m] area
        Sref = r[2]
        V_som = 343#[m/s] sound speed
        M = V/V_som

        Rwf =1.08 #declarado na função de CD_wing
        lf = fus[0] #[m] comprimento da fuselagem (e tail boom) Tab 5.1
        df = fus[1]#[m] altura da fuselagem
        R_Nf = p*U_l*lf/u
        aux3 = np.log10(R_Nf)**2.588
        aux4 = (1+0.144*M**2)**0.58
        cff	= 0.455/(aux3*aux4)
        Swet_f = S#0.05 #Área molhada da fuselagem
        aux1 = 60/((lf/df)**3) 
        aux2 = 0.0025*(lf/df)
        CD_0fus = Rwf*cff*(1+aux1+aux2)*(Swet_f/Sref)
        
        '''        
        S_bfus = 0#[m^2]Área da base da fuselagem
   
        if S_bfus == 0:
            CD_bfus = 0#somente para aviões com corte no tail boom ou tail boom de vigas, formando uma base
        else:
            CD_bfus = (0.029*(db/df)**3/(CD_0fusbase*(S/S_fus))**(1/2))*Sfus/S

        alfa = 0#Ângulo de arfagem
        Mc = M*np.sin(alfa)#sen(alfa) #Número de Mach de um escoamento cruzado
        n = 0.5 #razão do arrasto de um cilindro fiinito para o arrasto de um cilindro infinito, pegar valor no livro Fig. 5.19 por meio do valor de lf/df
        c_dc = #arrasto de escoamento cruzado experimental Livro FIg 5.20
        S = r[2]#0.365#[m]Área
        S_bfus = 0#[m^2]Área da base da fuselagem
        S_plffus = #Livro Fig 4.17
        CD_Lfus = 2*(alfa**2)*S_bfus/S+n*c_dc*abs(alfa)**3*S_plffus/S
        '''
        CD_bfus = 0.087# olhar figura 5.20
        C_fus = CD_0fus# + CD_bfus #+ CD_Lfus

        return C_fus
    fuselage = CD_fus(r,V,fus,Asa)
    
        #---------- CD_empennage (pág. 66 ~ 69)---------
    #colcor parte da asa

        #---------- CD_trim --------------
    def CD_trim():
        
        cd_trim = 0

        return cd_trim
    trim = CD_trim()

        #---------- CD_landingGear -------
    def CD_LandGear(r,La):
        a_lg = La
        Prop_area =2.884 * 0.092903#[m²] converter ft² para m²
        CD_LG = (0.016/Prop_area)*a_lg
        return CD_LG
    LandingGear = CD_LandGear(r,la)
    
    	#--------- CD_interf ------------
    def CD_interf():
        
        Cd_interf = 0
        return Cd_interf
    interf = CD_interf()



    	#--------- CD_misc --------------
    def CD_misc():
        Cd_misc = 0

        return Cd_misc
    misc = CD_misc()

    CD = (wing+fuselage+empennage+LandingGear + interf+misc+trim)*(1+s/100)
    
    if plotCD == True:    
        
        print('\nVelocity: ',round(V,2),'[m/s]')
        print('\nWing Drag = ',round(wing,3))
        print('Horizontal Estb. Drag = ',round(EH,3))
        print('Vertical Estb. Drag = ',round(EV,3))
        print('Fuselage Drag = ',round(fuselage,3))
        print('Landing Gear Drag = ',round(LandingGear,6))
        print('Total Aircraft Drag = ',round(CD,3))
        print('')
    return CD,EH,EV
#fim

'''
def CD_emp (r,V,Cl,V_stall,V,a):
        Rwf_emp = 1.0 #valor do livro p.66
        p = r[0]#1.09 #[kg/m^3] massa específica do ar
        u = r[1]#1.78*10**(-5) #u viscosidade dinâmica do ar
        U_l = V_stall#14. #[m/s] U_l = velocidade de estol
        S = r[2]#0.365714 #[m] Área da asa
        b = r[3]#1.6 #envergadura da asa

        #------ Empenagem horizontal -------#
        i_h = 2#número de superfícies de empenagens horizontais
        Rls_ht = 1.08#[] fator de correção de superficie Fig. 5.12 (Livro Fig4.2) determinar pelo ângulo de enflexamento

        x_t = a[1]
                if x_t>= 0.3:
                	L_l = 1.2         #L_l = # [] parametro de espessura local fig5.15
                else:
                	L_l = 2.0

                tc = a[0]#[] espessura  relativa
                if tc>=0.05:
                    K_w = 1.9767 + 0.5333*tc
                else:
                    K_w = 2.0

        S_exp = S #[m²] area de referencia (no caso do Aerodesign S_exposed = S)
        Swet_w = K_w*S_exp #[m²] area molhada

        Swet_ht = = K_w*S_exp #[m²]#Área molhada do profundor
        c_ht = S_ht/b_ht #c_we = corda geométrica média da parte exposta da asa
        U_l = V# 14.0#[m/s]velocidade de vôo
        M = U_l/346
        R_Nht = p*U_l*c_ht/u #valor a ser usado para fig 5.11
        Cfht = 0.455/((np.log10(R_Nht))**2.58)*((1+ 0.144*M**2)**0.58)

        CD_0ht = Rwf_emp*Rls_ht*Cfht*(1+L_l*tc+100*100/(tc**4))*(Swet_ht/S)

        #------ Empenagem Vertical -------#
        i_v = #número de superfícies de empenagens verticais
        Rls_vt =#[] fator de correção de superficie Fig. 5.12 (Livro Fig4.2) determinar pelo ângulo de enflexamento
        Swet_vt = #Área molhada do profundor
        c_vt = S_vt/b_vt #c_we = corda geométrica média da parte exposta da asa
        U_l = V #14.0#[m/s]velocidade de vôo
        M = U_l/346
        R_Nvt = p*U_l*c_vt/u #valor a ser usado para fig 5.11
        Cfvt = 0.455/((np.log10(R_Nvt))**2.58)*((1+ 0.144*M**2)**0.58)

        CD_0vt = Rwf_emp*Rls_vt*Cfvt*(1+L_l*tc+100*100/(tc**4))*(Swet_vt/S)
        if canard ==  False:
            CD_Lemp = 0

        CD_emp = CD_0ht + CD_0vt + CD_Lemp
        return CD_emp
    empennage = CD_emp ()
'''
