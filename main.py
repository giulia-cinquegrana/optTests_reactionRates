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

        tau = {'ni56': 5.251e5,
        'co56': 6.672e6,
        'ni57': 1.28e5,
        'co57': 2.35e7,
        'ti44': 1.86e9,
        'sc44': 1.43e4
        }

        self.Y = np.array([1.0, 0.0, 0.0])
        self.t = 0.0
        self.hl = np.array([tau['ni56'], tau['co56']])

        pass
    
    def step(self, dt, adaptive=True, tol_abs=1e-18, tol_rel=0.02):
        '''
        atol: pick the smallest abundance change you care about, e.g. 10^-20, 10^-18.
        '''
        # go one step
        # currently euler (p=1)
        p = 1
        s = 0.9  # safety factor
        
        while True:

            dYdt = self.calc_decay(self.Y, self.hl)

            Y_trial = self.Y + dt * dYdt
            delta = np.abs(Y_trial - self.Y)
            scale_i = tol_abs + tol_rel * np.maximum(np.abs(self.Y),np.abs(Y_trial))
            error_i = np.abs(delta) / scale_i
            error_max = float(np.max(error_i))

            if adaptive and (error_max > 1):
                fr = max(0.2, min(2.0, s * error_max**(-1/p)))
                dt *= fr
                continue

            # where error_max < 1, accept
            self.Y = np.clip(Y_trial, 0.0, None)
            self.t += dt

            if adaptive:
                # if last step easy (small error_max), grow dt.
                # if last was tricky (large error_max), shrink dt.
                # max_growth = increased by factor of 2
                # min_shrink = halved
                growth = min(2.0, max(0.5, s / max(error_max, 1e-16))) # max(error_max, 1e-16) = avoid divide by zero if error_max is perfect
                dt_next = dt * growth
            else:
                dt_next = dt

            return dt, dt_next, error_max

r = reaction_rates()
print(r.step(0.1))