import math
import numpy as np
import datetime
import sys

def is_inside_torus(x, y, z, a, c):
    if math.pow((c - math.sqrt(x * x + y * y)), 2) + z * z <= a * a and x*x + y*y >= c*c:
        return True
    else:
        return False

def PDBister(allCharArray, valueCharArray, positionRight):
    if (positionRight > len(allCharArray)):
        print ("Wrong index!")
        return ('WRONG INDEX!')
    positionRight += 1
    positionLeft = positionRight - len(valueCharArray)
    num = 0
    for i in range (positionLeft, positionRight):
        allCharArray[i] = valueCharArray[num]
        num += 1
    return allCharArray

def damemb(rDam, dmax, r1, r2, r3, delta, zCorona, name):
    pdb = open(name, 'w')
    header = "Spherical search volume created by DAMEMB   ---   "
    header += str(datetime.datetime.now())
    header += "\nDmax = " + str(dmax)
    header += "\nDAM packing radius = " + str(int(rDam))
    header += "\nTransmembrane radius = " + str(r1)
    header += "\nHydrophobic radius = " + str(r2)
    header += "\nHydrophilic radius = " + str(r3)
    header += "\nVertical position of the corona = " + str(zCorona)
    header += "\nSymmetry used: P1"
    header += "\nNumber of phases: 3"
    pdb.write(header)
    halfAxis = (int)(dmax/2.0)
    ax = np.arange(-halfAxis, halfAxis, rDam)
    ax = np.around(ax, decimals = 3)
    number = 0
    chain = 0
    for x in ax:
        chain += 1
        for y in ax:
            for z in ax:
                #uncomment it for a spherical searching volume
                #if x*x + y*y + z*z <= halfAxis*halfAxis:
                    number += 1
                    st = list('\nATOM                                                                        ')
                    st = PDBister(st, list(str(number)), 11)
                    st = PDBister(st, list('ASP'), 20)
                    st = PDBister(st, list(str(chain)), 25)
                    st = PDBister(st, list('1.00'), 60)
                    st = PDBister(st, list('20.00'), 66)
                    st = PDBister(st, list(str(x)), 38)
                    st = PDBister(st, list(str(y)), 46)
                    st = PDBister(st, list(str(z)), 54)
                    if x * x + y * y <= (r1 - delta) * (r1 - delta) and z <= (zCorona + r1) and z >= (zCorona - r1):
                        # region of protein
                        # 1 - fixed, 0 - free
                        st[68] = '1'
                        # ordinal number from iAllo
                        st[70] = '1'
                        # number of allowed phases
                        st[72] = '1'
                        # iAllo
                        st[73:77] = '1     '
                        # atom type in the damesv-mod.pdb
                        st = PDBister(st, list('CA'), 15)
                    elif z > (zCorona + r1) or z < (zCorona - r1):
                        # region of protein or solvent beyond the corona
                        # 1 - fixed, 0 - free
                        st[68] = '0'
                        # ordinal number from iAllo
                        st[70] = '2'
                        # number of allowed phases
                        st[72] = '2'
                        # iAllo
                        st[73:77] = '01    '
                        # atom type in the damesv-mod.pdb
                        st = PDBister(st, list('SA'), 15)
                    elif is_inside_torus(x, y, z - zCorona, 2 * delta, r1 - delta):
                        # region of protein or tails
                        # 1 - fixed, 0 - free
                        st[68] = '0'
                        # ordinal number from iAllo
                        st[70] = '1'
                        # number of allowed phases
                        st[72] = '2'
                        # iAllo
                        st[73:77] = '12    '
                        # atom type in the damesv-mod.pdb
                        st = PDBister(st, list('SA'), 15)
                    elif is_inside_torus(x, y, z - zCorona, (r2 - r1) / 1.0, r1 - delta):
                        # region of tails
                        # 1 - fixed, 0 - free
                        st[68] = '1'
                        # ordinal number from iAllo
                        st[70] = '1'
                        # number of allowed phases
                        st[72] = '1'
                        # iAllo
                        st[73:77] = '2     '
                        # atom type in the damesv-mod.pdb
                        st = PDBister(st, list('N'), 14)
                    elif is_inside_torus(x, y, z - zCorona, (r2 - r1) / 1.0 + 2 * delta, r1 - delta):
                        # region of heads or tails
                        # 1 - fixed, 0 - free
                        st[68] = '0'
                        # ordinal number from iAllo
                        st[70] = '1'
                        # number of allowed phases
                        st[72] = '2'
                        # iAllo
                        st[73:77] = '23    '
                        # atom type in the damesv-mod.pdb
                        st = PDBister(st, list('SA'), 15)
                    elif is_inside_torus(x, y, z - zCorona, (r3 - r1) / 1.0, r1 - delta):
                        # region of heads
                        # 1 - fixed, 0 - free
                        st[68] = '1'
                        # ordinal number from iAllo
                        st[70] = '1'
                        # number of allowed phases
                        st[72] = '1'
                        # iAllo
                        st[73:77] = '3     '
                        # atom type in the damesv-mod.pdb
                        st = PDBister(st, list('O1'), 15)
                    elif is_inside_torus(x, y, z - zCorona, (r3 - r1) / 1.0 + 2 * delta, r1 - delta):
                        # region of heads or solvent
                        # 1 - fixed, 0 - free
                        st[68] = '0'
                        # ordinal number from iAllo
                        st[70] = '2'
                        # number of allowed phases
                        st[72] = '2'
                        # iAllo
                        st[73:77] = '03    '
                        # atom type in the damesv-mod.pdb
                        st = PDBister(st, list('SA'), 15)
                    else:
                        # region of solvent
                        # 1 - fixed, 0 - free
                        st[68] = '1'
                        # ordinal number from iAllo
                        st[70] = '1'
                        # number of allowed phases
                        st[72] = '1'
                        # iAllo
                        st[73:77] = '0123  '
                        # atom type in the damesv-mod.pdb
                        st = PDBister(st, list('H'), 14)
                    #print("".join(st))
                    pdb.write("".join(st))

    pdb.close()
    print(str(name) + ": " + str(number) + " beads have been created")
if len(sys.argv) == 9:
    damemb(float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), \
    float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), str(sys.argv[8]))
elif len(sys.argv) == 2 and (sys.argv[1]) == '-h':
    print('''  damemb [<R_Dam>] [<Dmax>] [R_memb]
  [R_tails] [R_heads] [Delta] [Z_corona] [Out_name]
  The program creates starting searching volume
  for MONSA for membrane proteins with amphiphilic corona
  around their hydrophilic part.
  Parameters are: 
  [1] - Radius of a dummy atom, R_dam
  [2] - Side of a square OR diameter of a sphere, Dmax
  [3] - Distance from the origin to the end of protein phase
  [4] - Distance from the origin to the beginning of tail phase
  [5] - Distance from the origin to the beginning of head phase
  [6] - One half of a boundary region thickness between phases, Delta
  [7] - Vertical position of a corona, zCorona
  [8] - Output filename

  Example:                                       

  python damemb.py 3.7 132.4 14.3 34.3 48.2 3.0 20.1 output.pdb''')
else:
    print("Wrong number of provided parameters!")
    print("Required: 8")
    print("Provided: " + str(len(sys.argv) - 1))
