#from sklearn import datasets, svm, metrics, linear_model
from sklearn.neural_network import MLPRegressor
import random
import os, sys
import numpy as np
import NanoCore as nc
from link import link


def mycmp(a, b): return cmp(a[0], b[0])


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def predict(atom_category, input_list, broad=0.05, erange=10):

    # descriptor options
    r_cut = 10
    sigma = 0.02
    N_grid = int(r_cut * 50) + 1
    dr = r_cut / N_grid

    # Initialize Data set path
    #path = "../../Graphene_6x6_doping/" + atom_category # test
    path = "../Graphene_6x6_doping/" + atom_category    # django

    # Initialize default value
    data_set = []
    data_Ef = []
    target_list = []
    data_set_struct = []
    data_set_struct2 = []

    # reference structure for serial numbers
    ref_xyz = 'ml/ref_coord2.xyz' # django
    #ref_xyz = 'ref_coord2.xyz'    # test
    ref_at = nc.io.read_xyz(ref_xyz)
    tol = 0.01 # coordination difference criterion

    # number of dopant atoms
    temp = input_list.split(',')
    n_dop = len(temp)

    # Parse Data
    for i in os.listdir(path):
        for j in os.listdir(path + '/' + i):
            if ( j == 'DOS_%3.2f_%2.2i' % (broad, erange) ) and ( int(i.split('_')[2]) == n_dop ):
            #if j == 'DOS_%3.2f_%2.2i' % (broad, erange):
                print i

                # Parse Ef, DOS data
                f = open(path+'/'+i+'/'+j, 'r')
                f_lines = f.readlines()
                data_Ef.append( float(f_lines[9].split()[3]) )
                result = []
                for line in f_lines[13:]:
                    data = line.split('\n')[0].split(' ')
                    temp = []
                    for k in data:
                        if k != '':
                            temp.append(k)
                    result.append(np.float64(temp[1]))
                f.close()
                data_set.append(result)

                # Parse Struct data
                at = nc.siesta.read_fdf(path+'/'+i+'/STRUCT.fdf')

                ##################
                # new descriptor #
                ##################

                # at1 : original cell including C only
                # at2 : 5x5 supercell including B/N only
                at1 = at.copy()
                at2 = at * (5,5,1)
                at2.select_elements('C')
                at2.delete()
                at2.set_serials(1)
                at1.select_elements('C')
                at1.select_reverse()
                at1.delete()
                at1.set_serials(1)

                # shift the original cell to the center of (5x5) cell
                cell = at1.get_cell()
                v_tr = 2*np.array(cell[0]) + 2*np.array(cell[1])
                at1.select_all()
                at1.translate(*v_tr)

                # output array
                prdf = np.zeros(N_grid)

                # pick first atom from unitcell
                for atom1 in at1:
                    c1 = atom1.get_position(); s1 = atom1.get_symbol()
                
                    #pick second atom
                    for atom2 in at2:
                        c2 = atom2.get_position(); s2 = atom2.get_symbol()
                        distance = abs(c2-c1)
                        prdf += gaussian( np.linspace(0, r_cut, N_grid), distance, sigma ) * (1/distance**2)
                data_set_struct2.append(prdf)
                #print "prdf =", prdf

                ##################
                # new descriptor #
                ##################

                ### logical descriptor ###
                i_str = 0
                struct_result = []

                for atm in at:
                    symb = atm.get_symbol()
                    x, y, z = atm.get_position()

                    # find the real serial number
                    serial = 1
                    for atm1 in ref_at:
                        x1, y1, z1 = atm1.get_position()
                        if abs(x1-x) < tol and abs(y1-y) < tol and abs(z1-z) < tol:
                            break
                        else:
                            serial += 1

                    # C or B/N?
                    if symb == 'C':
                        struct_result.append([serial, 2])
                    else:
                        struct_result.append([serial, 1])
                    i_str += 1

                struct_result.sort(mycmp)
                struct_result = list( np.array( np.array(struct_result)[:,1] ) )
                data_set_struct.append(struct_result)
                ### logical descriptor ###

                # Add target_list
                target_list.append(i)


    # Create a classifier: a simple linear model
    #classifier = linear_model.Ridge(alpha=0.5)

    #hidden_layer_sizes = (10, 20, 40, 20, 10)
    hidden_layer_sizes = (N_grid, 2*N_grid, N_grid)
    classifier = MLPRegressor(solver='adam', alpha=1e-5,
                              hidden_layer_sizes=hidden_layer_sizes, 
                              random_state=1)

    # Input
    c_list = []
    for i in range(72):
        c_list.append('2')
    temp = input_list.split(',')
    n_dop = len(temp)
    for item in temp:
        c_list[int(item)] = '1'
    link.sort(mycmp)
    link_arr = np.array(link)[:,1]

    c_list2 = np.array(np.zeros(72), dtype=int)
    i_c = 0
    for l in link_arr:
        c_list2[l-1] = c_list[i_c]
        i_c += 1
    c_list = list( np.array(c_list2, dtype=int) )
    print "input list:", c_list

    #############
    # new Input #
    #############
    i_c = 0
    at_input = []
    for serial in c_list:
        x,y,z = ref_at[i_c].get_position()
        symb = ''
        if serial == 2:
            symb = 'C'
        elif serial == 1:
            symb = str(atom_category)
        #print str(symb)
        atom = nc.Atom(symb, [x,y,z])
        at_input.append(atom)
        i_c += 1
    cell = ref_at.get_cell()
    xx = cell[0][0]; xy = cell[0][1]; yy = cell[1][1]
    cell[0][0] = yy
    cell[1][1] = xx
    cell[0][1] = 0.
    cell[1][0] = xy
    at_input = nc.AtomsSystem(at_input, cell=cell)

    #print at_input

    at2_input = at_input * (5,5,1)
    at2_input.select_elements('C')
    at2_input.delete()
    at2_input.set_serials(1)

    at_input.select_elements('C')
    at_input.select_reverse()
    at_input.delete()
    at_input.set_serials(1)
    N_alpha = len(at)

    v_tr = 2*np.array(cell[0]) + 2*np.array(cell[1])
    at_input.select_all()
    at_input.translate(*v_tr)

    #print at_input
    #print at2_input

    # output array
    prdf_input = np.zeros(N_grid)

    # pick first atom from unitcell
    for atom1 in at_input:
        c1 = atom1.get_position(); s1 = atom1.get_symbol()
    
        #pick second atom
        for atom2 in at2_input:
            c2 = atom2.get_position(); s2 = atom2.get_symbol()
            distance = abs(c2-c1)
            prdf_input += gaussian( np.linspace(0, r_cut, N_grid), distance, sigma ) * (1/distance**2)

    #################
    # end new Input #
    #################

    # identical configuration?
    i = 0
    ex_i = -1
    for struct in data_set_struct:
        if struct == c_list:
            ex_i = i
            break
        i += 1
    print "ex_i =", ex_i

    #print "PRDF_inp =", prdf_input, np.sum(prdf_input)
    #print "PRDF_org =", data_set_struct2[ex_i], np.sum(data_set_struct2[ex_i])
    print "d(PRDF) =", np.sum(np.abs(prdf_input - data_set_struct2[ex_i]))

    expected_data = np.zeros(1001)

    if ex_i == -1:
        print "OUT: No identical configuration"
        classifier.fit( data_set_struct2, data_set)

    elif ex_i == 0:
        print "OUT: Exclude the identical configuration [0]"
        classifier.fit( data_set_struct2, data_set)
        expected_data = data_set[ex_i]

    else:
        print "OUT: Exclude the identical configuration"
        classifier.fit( list(data_set_struct2[:ex_i]) + list(data_set_struct2[ex_i+1:]),
                        list(data_set[:ex_i]) + list(data_set[ex_i+1:]) )
        expected_data = data_set[ex_i]

    #predicted = classifier.predict( np.array(c_list) )[0]
    predicted = classifier.predict( [prdf_input] )[0]
    print "predicted,", predicted

    n_factor = np.sum(expected_data)
    diff = np.sqrt( (np.array(predicted)-np.array(expected_data))**2 )
    print "Accuracy =", 100-(sum(diff)/n_factor*100)

    index = -1
    for i in range(len(data_set_struct)):
        if data_set_struct[i] == c_list:
            index = i

    # Find input list in shuffled list
    #return [predicted, expected_data, prdf_input, data_set_struct2[ex_i]]
    return [predicted, expected_data]

