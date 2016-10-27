import numpy
import math
import random

MOVES = numpy.array(
    [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0], [1, 0, 1], [0, 1, 1], [1, 1, 0], [-1, 0, -1],
     [0, -1, -1], [-1, -1, 0], [1, 0, -1], [0, 1, -1], [1, -1, 0], [-1, 0, 1], [0, -1, 1], [-1, 1, 0]])

R = 25
# 2 x radius + a fringe, because lamin barrier has to be hermetic
BOUND = 2 * R + 2
DIST = 3
SQRT_TWO = math.sqrt(2)

EMPTY = 0
BINDER = 1
LAMIN = 2
BSITE_R = 3
BSITE_L = 4
REGDNA = 5


def find_components(movement):  # finds the two vector components of a diagonal movement
    first = numpy.zeros(3, dtype=numpy.int)
    second = numpy.zeros(3, dtype=numpy.int)
    if movement[0] != 0:
        if movement[1] != 0:
            return first + [movement[0], 0, 0], second + [0, movement[1], 0]
        else:
            return first + [movement[0], 0, 0], second + [0, 0, movement[2]]
    else:
        return first + [0, movement[1], 0], second + [0, 0, movement[2]]


def is_collision(chain, new_atom):
    return any(numpy.equal(chain, new_atom).all(axis=1))


def dist_from_mi(x, y, z, mi):
    return math.sqrt((x - mi) ** 2 + (y - mi) ** 2 + (z - mi) ** 2)


def dist(pos_1, pos_2):
    return math.sqrt(sum((pos_2 - pos_1) ** 2))


def find_center(chain):
    return numpy.mean(chain, axis=0)


def getStateWithLamins(bound, f):
    ext_bound = bound + 6
    state = numpy.zeros((ext_bound, ext_bound, ext_bound), dtype=numpy.int)
    MIDDLE = int(ext_bound / 2)
    r = MIDDLE - 2
    lam_name = f.split('.pdb')[0] + '_laminfin.pdb'
    save_lam = open(lam_name, "w+")
    save_lam.write("HEADER LAMINA")
    at_nr = 1

    for x in range(ext_bound):
        for y in range(ext_bound):
            for z in range(ext_bound):
                dist_sq = (x - MIDDLE) ** 2 + (y - MIDDLE) ** 2 + (z - MIDDLE) ** 2
                if ((r + 1) ** 2 <= dist_sq) and (dist_sq <= (r + 2) ** 2):
                    state[x, y, z] = LAMIN
                    # print border
                    line = "\nATOM  " + str(at_nr).rjust(5) + " " + "P".center(4) + " " + "LAM" + "  " + str(
                        at_nr).rjust(4) + "    " + str(round(x * DIST, 3)).rjust(8) + str(round(y * DIST, 3)).rjust(
                        8) + str(round(z * DIST, 3)).rjust(8) + "  0.00 00.00"
                    at_nr += 1
                    save_lam.write(line)

    save_lam.close()
    return state


def is_intersecting(chain, previous, movement):
    first, second = find_components(movement)

    # find the position in chain of previous + first; if such atom doesn't exist
    # position_of_prev_first will be an empty array
    # the same applies for position_of_prev_second
    position_of_prev_first = numpy.where((chain == previous + first).all(axis=1))[0]
    position_of_prev_second = numpy.where((chain == previous + second).all(axis=1))[0]
    # print position_of_prev_first, position_of_prev_second
    # check if both atoms are in the chain
    if position_of_prev_first.size + position_of_prev_second.size != 2:
        return False  # if they are not, then there cannot be an intersection
    else:  # both atoms are in the chain
        if ((position_of_prev_first[0] == position_of_prev_second[0] + 1) or
                (position_of_prev_first[0] == position_of_prev_second[0] - 1)):
            return True  # the positions are consecutive, hence there is an intersection
        else:
            return False  # the positions are not consecutive, there is no intersection


def intersection(new_pos, preceding, following, chain):
    if (dist(new_pos, preceding) == 1.0) or (dist(new_pos, following) == 1.0):
        return False
    else:
        if dist(new_pos, preceding) == math.sqrt(2):
            if is_intersecting(chain, preceding, new_pos - preceding):
                return True
        if dist(new_pos, following) == math.sqrt(2):
            if is_intersecting(chain, new_pos, following - new_pos):
                return True
        return False


def check_collisions(chain):
    t = False
    for rec in chain:

        if is_collision(chain, rec):
            w = numpy.where((chain == rec).all(axis=1))[0]
            if w.size > 1:
                t = True
                print "COLLISON!", w
    return t


def modify(chain, k=-1):
    i = random.randint(0, len(chain) - 1)
    if k != -1:
        i = k
    move = random.choice(MOVES)
    new = chain[i] + move

    if is_collision(chain, new):
        # print "COLLISION!!!"
        # print new, numpy.where((chain == new).all(axis=1))[0]
        return

    if (i != 0) and (i != len(chain) - 1):
        if (dist(new, chain[i - 1]) <= math.sqrt(2)) and (dist(new, chain[i + 1]) <= math.sqrt(2)):
            if intersection(new, chain[i - 1], chain[i + 1], chain):
                return
            else:
                return i, move
        else:
            return
    elif i == 0:
        if dist(new, chain[1]) <= math.sqrt(2):
            if dist(new, chain[1]) == 1.0:
                return i, move
            else:
                if is_intersecting(chain, new, chain[1] - new):
                    return
                else:
                    return i, move
        else:
            return
    else:  # i == len(chain) - 1
        if dist(new, chain[len(chain) - 2]) <= math.sqrt(2):
            if dist(new, chain[len(chain) - 2]) == 1.0:
                return i, move
            else:
                if is_intersecting(chain, chain[len(chain) - 2], new - chain[len(chain) - 2]):
                    return
                else:
                    return i, move
        else:
            return


def modify_out(chain, state, lens):
    i = random.randint(0, len(chain) - 1)
    move = random.choice(MOVES)
    new = chain[i] + move

    if state[tuple(new)] != 0:
        return

    beginnings = [0]
    current_pos = 0
    for l in lens[:-1]:
        current_pos += l
        beginnings.append(current_pos)

    endings = []
    current_pos = 0
    for l in lens:
        current_pos += l
        endings.append(current_pos - 1)
    # print beginnings
    # print endings
    if (i not in beginnings) and (i not in endings):
        if (dist(new, chain[i - 1]) <= math.sqrt(2)) and (dist(new, chain[i + 1]) <= math.sqrt(2)):
            if intersection(new, chain[i - 1], chain[i + 1], chain):
                return
            else:
                return i, move
        else:
            return
    elif i in beginnings:
        if dist(new, chain[i + 1]) <= math.sqrt(2):
            if dist(new, chain[i + 1]) == 1.0:
                return i, move
            else:
                if is_intersecting(chain, new, chain[i + 1] - new):
                    return
                else:
                    return i, move
        else:
            return
    else:  # i in endings
        if dist(new, chain[i - 1]) <= math.sqrt(2):
            if dist(new, chain[i - 1]) == 1.0:
                return i, move
            else:
                if is_intersecting(chain, chain[i - 1], new - chain[i - 1]):
                    return
                else:
                    return i, move
        else:
            return


def av_dist(chain, center):
    dist_sum = 0
    for atom in chain:
        # print atom, dist(atom, center)
        dist_sum += dist(atom, center)

    return float(dist_sum) / float(len(chain))


def metropolis_out(chain, state, lens, iters):
    print "METROPOLIS_OUT"
    count = 0

    for j in range(iters):
        modification = modify_out(chain, state, lens)
        if modification:
            # print modification
            i, move = modification
            state[tuple(chain[i])] = EMPTY
            chain[i] = chain[i] + move
            state[tuple(chain[i])] = REGDNA
        if (j % 400000) == 0:
            print j
        if (j % 10000000) == 0:
            beg = 0
            for i in range(1, 7):
                end = sum(lens[0:i])
                write_pdb(chain[beg:end], "chain_" + str(i) + "_step_" + str(j))
                beg = end
    return chain


def addChainsToState(chain, state):
    for elem in chain:
        if state[tuple(elem)] != 0:
            print "Error - place occupied by two objects - please check", elem
        else:
            state[tuple(elem)] = REGDNA
    return state


def metropolis(chain, required_av_dist, delta,
               mode):  # the larger the delta, the less likely to make distance-increasing moves
    # print "METROPOLIS"
    center = find_center(chain)
    count = 0
    c = 0
    distances = []
    print "Initial average distance of a part from the center of the molecule:", av_dist(chain, center)
    w = -1
    while True:
        if count < (len(chain) * 150) and mode == 1:
            modification = modify(chain, count % len(chain))
            w = 1
        else:
            modification = modify(chain)
            w = 2

        count += 1
        if modification:
            i, move = modification
            chi = chain[i] + move
            # print chi, chain[i]
            old_dist, new_dist = dist(chain[i], center), dist(chi, center)
            # print old_dist, new_dist
            # print old_dist - new_dist
            c += 1
            # print old_dist, new_dist

            if new_dist < old_dist:
                distances.append(old_dist - new_dist)
                chain[i] = chi
            else:  # distance increased
                if random.uniform(0.0, 1.0) < math.exp((old_dist - new_dist) * delta):
                    distances.append(old_dist - new_dist)
                    chain[i] = chi

        if count % 500000 == 0:
            newdist = av_dist(chain, center)
            print "Current average distance:", newdist  # , w, c, count
            if newdist < required_av_dist:
                print "SUCCESS - required average distance reached"
                break
                # if count % 300000 == 0:
                # print "writing", count
                # write_pdb(chain, "test" + str(count))
    return chain, distances


def position_chains(chain1, chain2, chain3, chain4, chain5, chain6, r):
    # print find_center(chain1).astype(numpy.int)
    k = 14
    chain1 = chain1 - find_center(chain1).astype(numpy.int) + [0, 0, -k]
    chain2 = chain2 - find_center(chain2).astype(numpy.int) + [k, 0, 0]
    chain3 = chain3 - find_center(chain3).astype(numpy.int) + [0, k, 0]
    chain4 = chain4 - find_center(chain4).astype(numpy.int) + [0, 0, k]
    chain5 = chain5 - find_center(chain5).astype(numpy.int) + [-k, 0, 0]
    chain6 = chain6 - find_center(chain6).astype(numpy.int) + [0, -k, 0]

    a = []
    a.extend(chain1)
    a.extend(chain2)
    a.extend(chain3)
    a.extend(chain4)
    a.extend(chain5)
    a.extend(chain6)

    a = numpy.array(a)

    # print check_collisions(a)
    chain1 += [r]
    chain2 += [r]
    chain3 += [r]
    chain4 += [r]
    chain5 += [r]
    chain6 += [r]

    write_pdb(chain1, "positioned_chain_1")
    write_pdb(chain2, "positioned_chain_2")
    write_pdb(chain3, "positioned_chain_3")
    write_pdb(chain4, "positioned_chain_4")
    write_pdb(chain5, "positioned_chain_5")
    write_pdb(chain6, "positioned_chain_6")
    return chain1, chain2, chain3, chain4, chain5, chain6


def random_init(chain_length, num_of_tries=100):
    chain = numpy.zeros((chain_length, 3), dtype=numpy.int)
    for i in range(1, chain_length):  # the first atom is in (0, 0, 0)
        tries = 0
        while True:
            tries += 1
            if tries >= num_of_tries:
                print "Try again, unable to initialize"
                return
            previous = chain[i - 1]
            movement = random.choice(MOVES)
            current = previous + movement
            if is_collision(chain, current):  # place cannot be occupied
                continue
            if current[2] < 0:  # atom has to be in the right half-space
                continue
            if sum(abs(movement)) > 1:  # when moving diagonally check for intersection
                if is_intersecting(chain, previous, movement):
                    continue

            chain[i] = current
            break

    return chain


def write_pdb(coordinates, name):
    def pdb_line(at_nr, at_name, desc, chain_n, res_nr, pos):
        return "\nATOM  " + str(at_nr).rjust(5) + " " + at_name.center(4) + " " + desc + " " + chain_n + str(
            res_nr).rjust(4) + "    " + str(round(pos[0] * DIST, 3)).rjust(8) + str(round(pos[1] * DIST, 3)).rjust(
            8) + str(round(pos[2] * DIST, 3)).rjust(8) + "  0.00 00.00"

    f = open(str(name) + '.pdb', 'w')
    f.write('HEADER 0 ' + str(len(coordinates)) + 'step 0\n')
    f.write('TITLE chromosome;bonds=5000')
    for i, atom in enumerate(coordinates):
        f.write(pdb_line(i + 1, 'C', 'UNB', 'A', i + 1, atom))
    for i in range(len(coordinates)):
        if i == len(coordinates) - 1:
            break
        f.write('\nCONECT' + str(i + 1).rjust(5) + str(i + 2).rjust(5))
    f.write('\nEND')


def find_maximal_distance(chain):
    maxim = 0
    for rec in chain:
        for re in chain:
            d = dist(re, rec)
            if d > maxim:
                maxim = d
    return maxim


def d_dista(chain):
    m = [0, 0, 0]
    for rec in chain:
        for re in chain:
            dis = abs(rec - re)
            if dis[0] > m[0]:
                m[0] = dis[0]
            if dis[1] > m[1]:
                m[1] = dis[1]
            if dis[2] > m[2]:
                m[2] = dis[2]
    return m
