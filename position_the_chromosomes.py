import pendolino
import pickle
import numpy
import math
import argparse


def dist(pos_1, pos_2):
    return math.sqrt(sum((pos_2 - pos_1) ** 2))


def max_d(chain, center):
    max_dist = 0
    for elem in chain:
        d = dist(elem, center)
        if d > max_dist:
            max_dist = d
    return max_dist


def loader(name):
    with open(name, "r") as f:
        return pickle.load(f)


parser = argparse.ArgumentParser()
parser.add_argument("-n", type=int, help="number of steps of decondensation", default=50000000)
args = parser.parse_args()
N = args.n

chain1 = loader("chain1")
chain2 = loader("chain2")
chain3 = loader("chain3")
chain4 = loader("chain4")
chain5 = loader("chain5")
chain6 = loader("chain6")

lens = [len(chain1), len(chain2), len(chain3), len(chain4), len(chain5), len(chain6)]
c1, c2, c3, c4, c5, c6 = pendolino.position_chains(chain1, chain2, chain3, chain4, chain5, chain6, 27)

a = []
a.extend(c1)
a.extend(c2)
a.extend(c3)
a.extend(c4)
a.extend(c5)
a.extend(c6)

# print "md", max_d(numpy.array(a), (0, 0, 0))

state = pendolino.getStateWithLamins(46, "lamin")

chain = numpy.array(a)
state = pendolino.addChainsToState(chain, state)

lens = [len(c1), len(c2), len(c3), len(c4), len(c5), len(c6)]

chain = pendolino.metropolis_out(chain, state, lens, N)
beg = 0
for i in range(1, 7):
    end = sum(lens[0:i])
    pendolino.write_pdb(chain[beg:end], "final_chain_" + str(i))
    beg = end

with open("final_chains", "w+") as f:
    pickle.dump(chain, f)
