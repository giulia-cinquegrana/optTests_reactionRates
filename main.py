import numpy as np
import matplotlib.pyplot as plt

class reaction_rates: 

    def __init__(self, chains, hl, Y_init):

        self.t = 0.0
        self.chains = list(chains)
        self.idx = {label:i for i, label in enumerate(self.chains)}
        self.hl = np.asarray(hl, dtype=float)
        self.Y = np.asarray(Y_init, dtype=float)

        pass

    def calc_decay(self):

        '''
        simple radioactive decay, A --> B --> C

        Y: [Y_A, Y_B, Y_C], abundance (number fraction). 
        hl: [tau_A, tau_B]
        rf: reaction flux
        '''

        decay_const  = np.log(2) / self.hl

        rf1 = decay_const[:,0]* self.Y[:,0]   
        rf2 = decay_const[:,1]* self.Y[:,1]

        dY = np.zeros_like(self.Y)

        dY[:,0] = -rf1
        dY[:,1] = rf1 - rf2
        dY[:,2] = rf2

        return dY
    
    def step(self, dt, adaptive=True, tol_abs=1e-18, tol_rel=0.02, Y_gate=1e-12):
        '''
        travel one step
        euler, p=1
        atol: pick smallest abundance change you care about, e.g. 10^-20, 10^-18.
        s = safety factor
        '''

        p = 1
        s = 0.9
        
        while True:

            dY = self.calc_decay()

            Y_trial = self.Y + dt * dY

            # fractional change per species, referenced to current amount
            frac = np.abs(dt * dY) / np.maximum(self.Y, Y_gate)
            # only species that currently exist control step
            active = self.Y > Y_gate
            error_max = float(frac[active].max()) if np.any(active) else float(frac.max())

            if adaptive and (error_max > tol_rel):
                factor = max(0.2, min(2.0, s * (tol_rel / error_max)))  # shrink and retry
                dt *= factor
                continue

            # accept when error_max < 1
            self.Y = np.clip(Y_trial, 0.0, None)
            self.t += dt

            if adaptive:
                # if last step easy (small error_max), grow dt.
                # or if large error_max, shrink dt.
                growth = min(2.0, max(0.5, s * (tol_rel / max(error_max, 1e-30))))
                dt_next = dt * growth
            else:
                dt_next = dt

            return dt, dt_next, error_max
        
    def run(self, t_max, n_max, dt_init=None, dt_out=None):

        '''
        t_max = max simulation time
        n_max = max number of steps
        dt_init = initial time step
        '''

        n = 0 # step number
        smol_n = 1e-12

        ln2 = np.log(2.0)

        if dt_init is None:
            dt_init = 0.02 * (float(np.min(self.hl[:, 0])) / ln2)  

        if dt_out is None:
            dt_out = 0.1 * float(np.min(self.hl[:, 0]))  

        #-----storage
        t_out = [self.t]
        abund_out = [self.Y.copy()]
        t_next_out = self.t + dt_out
        dt = float(dt_init)
        
        #-----main
        while self.t < t_max and n < n_max:

            dt_try = min(dt, t_max - self.t, t_next_out - self.t)
            if dt_try <= 0.0:
                # if already at/past next output boundary, bump boundary and keep going
                t_out.append(self.t)
                abund_out.append(self.Y.copy())
                t_next_out += dt_out
                continue

            dt_acc, dt_next, err = self.step(dt_try, adaptive=True)
            n += 1

            while self.t >= t_next_out - smol_n*max(1.0, abs(t_next_out)):
                t_out.append(self.t)
                abund_out.append(self.Y.copy())
                t_next_out += dt_out

        if t_out[-1] != self.t:
            t_out.append(self.t)
            abund_out.append(self.Y.copy())

        return np.asarray(t_out), np.stack(abund_out, axis=0)

if __name__ == "__main__":

    chains = ["56", "57"]
    idx = {c:i for i,c in enumerate(chains)}
    N = len(chains)

    hl = np.empty((N,2), dtype=float)
    hl[idx["56"]] = [5.251e5, 6.672e6]
    hl[idx["57"]] = [1.28e5,  2.35e7]

    Y0 = np.zeros((N,3), dtype=float)
    Y0[idx["56"],0] = 1.0
    R0 = 0.03
    Y0[idx["57"],0] = R0 * Y0[idx["56"],0]

    sim = reaction_rates(chains, hl, Y0)

    t_max = 5.0 * hl[idx["57"],1]
    T, Ysnap = sim.run(t_max, n_max=50_000)

    k56, k57 = idx["56"], idx["57"]
    t_days = T / 86400.0
    Y56 = Ysnap[:, k56, :]
    Y57 = Ysnap[:, k57, :]

    ratio = Y57[:,2] / np.maximum(Y56[:,2], 1e-30)

    fig, ax = plt.subplots(1,2, figsize=(10,4))
    ax[0].plot(t_days, Y56[:,0], label="56Ni")
    ax[0].plot(t_days, Y56[:,1], label="56Co")
    ax[0].plot(t_days, Y56[:,2], label="56Fe")
    ax[0].set_xlabel("time [days]"); ax[0].set_ylabel("Y (56-chain)"); ax[0].legend()

    ax[1].plot(t_days, ratio)
    ax[1].set_xlabel("time [days]"); ax[1].set_ylabel("57Fe / 56Fe")
    fig.tight_layout(); plt.show()