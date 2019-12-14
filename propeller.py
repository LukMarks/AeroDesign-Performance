
class propeller:
    def __init__(self,hel):
        self.hel = hel
        if hel == '12X4':
            self.A = -0.0014
            self.B = 0.0246
            self.C = -0.6623
            self.D = 27.292
            self.rho = 1.19
        elif hel == '12X5':
            self.A = -0.0066
            self.B = 0.1805
            self.C = -2.0333
            self.D = 28.543
            self.rho = 1.19
        elif hel == '13X4':
            self.A = -0.0004#(original) #0.002 
            self.B = 0.0235#(original)#-0.0601 
            self.C = -0.9538#(original)#0.0159 
            self.D = 32.905#(original)#28.533 
            self.rho = 1.181
        elif hel == '13X5':
            self.A = -0.0001#(original) #0.002 
            self.B = 0.0031#(original)#-0.0601 
            self.C = -0.4994#(original)#0.0159 
            self.D = 27.284#(original)#28.533    
            self.rho = 1.181
        elif hel == '14X4':
            self.A = 0.002
            self.B = -0.0677
            self.C = -0.0256
            self.D = 28.123
            self.rho = 1.19
        elif hel == '14X5':
            self.A = -0.0049
            self.B = 0.1778
            self.C = -2.6099
            self.D = 33.926
            self.rho = 1.19
        else:
            print("\nError: CAN'T FIND A PROPELLER\n")
            
        