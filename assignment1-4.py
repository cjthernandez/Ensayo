# Homework 1 Numerics for Bioinformatics
# Deadline: October 31st
# Homework 4

# The library numpy is imported
import numpy as np

# We create vectors with our variables X_i, constants k_i, and what the educsts
# and products of our reaction system are.
variables = np.array(['X1', 'X2', 'X3', 'X4'])
constants = np.array(['k1', 'k2', 'k3', 'k4', 'k5'])
educts = np.array(['0', 'X1', 'X1+X4', 'X2', 'X2'])
products = np.array(['X1', '0', 'X2', '0', 'X1+X3'])


def stochiometricMatrix(reacEduct, reacProduct):
    # We go over every reaction from our system
    for i in range(len(reacEduct)):
        # We create the educt and product vectors for each reaction
        ed_vec = np.array([[0], [0], [0], [0]])
        pr_vec = np.array([[0], [0], [0], [0]])
        # We get the current reaction's educts and products
        ed = reacEduct[i]
        pr = reacProduct[i]
        # We go over each variable X_i of our system
        for j in range(len(variables)):
            # We get the current variable X_i
            var = variables[j]
            # We check if the variable X_i is either in the educts or products
            # of the current reaction
            if (ed == 0):
                ed_vec[j] = 0
            elif (var in ed):
                ed_vec[j] = 1
            if (pr == 0):
                pr_vec[j] = 0
            elif (var in pr):
                pr_vec[j] = 1
            # We estimate the stochiometric vector for the current reaction
            st_vec = pr_vec - ed_vec
        # We construct the stochiometric Matrix
        if (i == 0):
            st_mat = st_vec
        else:
            st_mat = np.hstack((st_mat, st_vec))
    return st_mat


def propensityVector(stochMat, konst, vari):
    # We estimate the form of the stochiometric matrix
    (r, c) = np.shape(stochMat)
    # An empty vector is created in which the propensity functions will be
    # stored
    propen_vec = []
    for j in range(c):
        # We go through the columns of the stochiometric matrix
        cur_col = stochMat[:, j]
        # We count the amount of -1, which is equivalent to the amount of
        # species on the educts
        n1 = list(cur_col).count(-1)
        # We identify the 0th order reactions and construct the reaction rate
        if (n1 == 0):
            propen_vec.append(str(konst[j]))
        # We identify the 1st order reactions and construct the reaction rate
        elif (n1 == 1):
            w = np.where(cur_col == -1)[0]
            propen_vec.append(str(konst[j] + "*" + vari[w[0]]))
        # We identify the 2nd order reactions and construct the reaction rate
        elif (n1 == 2):
            w = np.where(cur_col == -1)[0]
            propen_vec.append(str(konst[j] + '*' + vari[w[0]] + "*"
                                  + vari[w[1]]))
    return propen_vec


def evaluateState(para, curre, stoch_mat, prop_vec):
    # A vector to store the values of the reaction rates is created
    r_vec = np.ones((len(prop_vec)))
    # We go through the dictionary of parameters
    for p, vp in para.items():
        # We go through the propensity vector in order to see in which
        # reaction rate equations the parameters are
        for g in range(len(prop_vec)):
            r = prop_vec[g]
            # We multiply the value of the parameter in the equation where it
            # appears
            if (p in r):
                r_vec[g] *= vp
    # We go through the dictionary of parameters
    for c, vc in curre.items():
        # We go through the propensity vector in order to see in which
        # reaction rate equations the values of the current state are
        for g in range(len(prop_vec)):
            r = prop_vec[g]
            # We multiply the value of the parameter in the equation where it
            # appears
            if (c in r):
                r_vec[g] *= vc
    # We perform the matrix multiplication of our stochiometric matrix and the
    # propensity vector with the equations evaluated in the current state with
    # the parameter values given
    return np.matmul(stoch_mat, r_vec)


# Part a)
st_matrix = stochiometricMatrix(educts, products)
prop_vector = propensityVector(st_matrix, constants, variables)

np.savetxt('SMatrix.txt', fmt='%.0f', X=st_matrix, delimiter=',')
np.savetxt('RVector.txt', fmt='%s', X=prop_vector, delimiter=',')

# Part b)
# We create dictionaries for the values of parameters and the current state
paramet = {"k1": 5, "k2": 3, "k3": 12, "k4": 7, "k5": 3}
curr_state = {"X1": 5, "X2": 25, "X3": 15, "X4": 5}

ODEValues = evaluateState(paramet, curr_state, st_matrix, prop_vector)

np.savetxt('ODEValues.txt', fmt='%.0f', X=ODEValues, delimiter=',')
