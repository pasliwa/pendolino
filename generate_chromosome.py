import pendolino as p
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("-n", type=int, help="number of parts", default=1000)
parser.add_argument("-r", type=float, help="req_av", default=5.1)
parser.add_argument("-m", type=str, help="name", default="chain")
parser.add_argument("-d", type=int, help="delta", default=3)

args = parser.parse_args()
N = args.n
req = args.r
m = args.m
delta = args.d

chain = p.random_init(N)

center = p.find_center(chain)
av = p.av_dist(chain, center)
dista = p.d_dista(chain)
maximal = p.find_maximal_distance(chain)

p.write_pdb(chain, m + "A")

nchain, distances = p.metropolis(chain, req, delta, 1)

ncenter = p.find_center(nchain)
nav = p.av_dist(nchain, ncenter)
ndista = p.d_dista(nchain)
nmaximal = p.find_maximal_distance(nchain)

p.write_pdb(chain, m + "B")
print "Coordinates of the center Av. dist. from center\tMaximal dist. in consequent dimensions\tMaximal distance between two parts"
print center, av, dista, maximal
print ncenter, nav, ndista, nmaximal

with open(m, "w+") as f:
    pickle.dump(nchain, f)
