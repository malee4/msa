# from D-Wave Reservoir Management Demo

# import the packages
from dimod import BinaryQuadraticModel
from dwave.system import DWaveSampler, EmbeddingComposite

# insert the inputs here
pumps = [0, 1, 2, 3]
costs = [[36, 27],
    [56, 65],
    [48, 36],
    [52, 16]]
flow = [2, 7, 3, 8]
demand = 20
time = [0, 1]

# build variable for each pump
x = [[f'P{p}_AM', f'P{p}_PM'] for p in pumps]

# initiaize empty BQM ('BINARY' indicates 0 and 1 values)
bqm = BinaryQuadraticModel('BINARY') 

# add in objective
for p in pumps:
    for t in time:
        # add the variable if it doesn't yet exist. If it does, it adds bias to what is already there
        bqm.add_variable(x[p][t], costs[p][t])

# add in constraint 1 
# all pumps run at least once per day
for p in pumps:
    c1 = [(x[p][t], 1) for t in time] # bias = 1
    bqm.add_linear_inequality_constraint(
        c1,
        lb = 1,
        ub = len(pumps),
        lagrange_multiplier = 13,
        label = 'c1_pump_'+str(p))
# lb = lower bound
# ub = upper bound

# add in constraint 2
# at most 3 pumps can run at a time
for t in time:
    c2 = [(x[p][t], 1) for p in pumps]
    bqm.add_linear_inequality_constraint(
        c2,
        constant = -3,
        lagrange_multiplier = 1,
        label = 'c2_time_'+str(t)
    )

# add in constraint 3
# total daily flow satisfies daily demand
c3 = [(x[p][t], flow[p]) for t in time for p in pumps]
bqm.add_linear_equality_constraint(
    c3,
    constant = -demand,
    lagrange_multiplier = 28
)

# running directly on quantum computer
sampler = EmbeddingComposite(DWaveSampler())
sampleset = sampler.sample(bqm, num_reads = 1000)

sample = sampleset.first.sample
total_flow = 0
total_cost = 0

print("\n\tAM\tPM")
for p in pumps:
    printout = 'P' + str(p)
    for time in range(2):
        printout += "\t" + str(sample[x[p][time]])
        total_flow += sample[x[p][time]]*flow[p]
        total_cost += sample[x[p][time]]*costs[p][time]
    print(printout)

print("\nTotal flow:\t", total_flow)
print("Total cost:\t", total_cost, "\n")