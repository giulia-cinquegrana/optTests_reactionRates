# optTests_reactionRates

ractive decay chains - batch solver. currently python. 

-independent A->B->C chains (e.g., Ni56->Co56->Fe56 and Ni57->Co57->Fe57)

    lambda_A(/B), tau_A(/B): decay constant and half life of parent(/intermediate)
    
    lambda_A(/B) = ln2 / tau_A(/B)
    dY_A = - lambda_A * Y_A
    dY_B = lambda_A * Y_A - lambda_B * Y_B
    dY_C = lambda_B * Y_B

-adaptive dt (currently euler).
    -per-step fractional-change control 
    -e.g., define the max change for a given species per timestep you want to resolve. 

