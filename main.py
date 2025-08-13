import numpy as np
import matplotlib.pyplot as plt

class reaction_rates(): 

    @staticmethod
    def calc_decay(Y, hl):

        '''
        simple radioactive decay, A --> B --> C

        Y: [Y_A, Y_B, Y_C], abundance (number fraction). 
        hl: [tau_A, tau_B]
        '''

        decay_const  = np.log(2) / hl # in s^-1

        reac_flux1 = decay_const[0]* Y[0]   
        reac_flux2 = decay_const[1]* Y[1]

        dYA_dt = -reac_flux1
        dYB_dt = reac_flux1 - reac_flux2
        dYC_dt = reac_flux2

        return np.array([dYA_dt, dYB_dt, dYC_dt])

    def __init__(self):
        pass

    def SN(self):

        tau = {'ni56': 5.251e5,
               'co56': 6.672e6,
               'ni57': 1.28e5,
               'co57': 2.35e7,
               'ti44': 1.86e9,
               'sc44': 1.43e4
               }
        
        Y0  = [1.0, 0.0, 0.0]                    
        hl  = [tau['ni56'], tau['co56']]      

        return self.calc_decay(Y0, hl)
            
r = reaction_rates()
print(r.SN())