import numpy as np
import shape_bib as sovi
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

memory_length = 4

def f(x): return (x[0]**2)*(x[1]**2)

#negative ableitung
def f_div(x):
    return np.array([2.*x[0]*(x[1]**2), 2.*(x[0]**2)*(x[1])])

x_start = np.array([20., 20.])

bfgs_memory = sovi.bfgs_memory(np.zeros([memory_length, 2 ]),
                                 np.zeros([memory_length, 2 ]), memory_length, 0)

current_val = x_start
alpha = np.zeros(bfgs_memory.length)
diff_grad = np.zeros(2)
first_diff_grad = np.zeros(2)

bfgs_memory.update_grad(f_div(x_start))


while(np.linalg.norm(bfgs_memory.gradient[0]) > 0.0001):
    #SCHRITT CALC
    q = bfgs_memory.gradient[0]
    if (bfgs_memory.step_nr + 1 >= bfgs_memory.length):
        # bei voll besetzter memory werden alle Eintraege verwendet

        for i in range(bfgs_memory.length - 1):
            # Vorwaertsschleife
            i = i + 1
            diff_grad = (bfgs_memory.gradient[i - 1] - bfgs_memory.gradient[i])
            alpha[i] = np.dot(bfgs_memory.deformation[i - 1], q) / np.dot(diff_grad, bfgs_memory.deformation[i - 1])
            q = q - float(alpha[i]) * diff_grad

            # Reskalierung von q
        first_diff_grad = (bfgs_memory.gradient[0] - bfgs_memory.gradient[1])
        gamma = np.dot(first_diff_grad, bfgs_memory.deformation[0]) / np.dot(first_diff_grad, first_diff_grad)
        q = gamma * q

        for i in range(bfgs_memory.length - 1):
            # Rueckwaertsschleife
            i = i + 1
            diff_grad = (bfgs_memory.gradient[-(i + 1)] - bfgs_memory.gradient[-i])
            beta = np.dot(diff_grad, q) / np.dot(diff_grad, bfgs_memory.deformation[-(i + 1)])
            q = q + (float(alpha[-i]) - beta) * bfgs_memory.deformation[-(i + 1)]

    elif (bfgs_memory.step_nr == 0):
        # der erste BFGS-Schritt ist ein Gradientenschritt, U Gradient ist in negativer Richtung
        q = q

    else:
        # bei nicht voll besetzter memory werden lediglich die besetzten Eintraege verwendet
        for i in range(bfgs_memory.step_nr):
            # Vorwaertsschleife
            i = i + 1
            diff_grad = (bfgs_memory.gradient[i - 1] - bfgs_memory.gradient[i])
            alpha[i] = np.dot(bfgs_memory.deformation[i - 1], q) / np.dot(diff_grad, bfgs_memory.deformation[i - 1])
            q = q - float(alpha[i]) * diff_grad

        # Reskalierung von q
        first_diff_grad = (bfgs_memory.gradient[0] - bfgs_memory.gradient[1])
        gamma = np.dot(first_diff_grad, bfgs_memory.deformation[0]) / np.dot(first_diff_grad, first_diff_grad)
        q = gamma * q

        for i in range(bfgs_memory.step_nr):
            # Rueckwaertsschleife
            shift = (bfgs_memory.length - 1) - bfgs_memory.step_nr
            i = i + 1
            diff_grad = (bfgs_memory.gradient[-(i + 1) - shift] - bfgs_memory.gradient[-i - shift])
            beta = np.dot(diff_grad, q) / np.dot(diff_grad, bfgs_memory.deformation[-(i + 1) - shift])
            q = q + (float(alpha[-i - shift]) - beta) * bfgs_memory.deformation[-(i + 1) - shift]

    q = -1.*q


    current_val = current_val + q


    bfgs_memory.step_nr = bfgs_memory.step_nr + 1
    bfgs_memory.update_grad(f_div(current_val))
    bfgs_memory.update_defo(q)
    print("curv cond: {0:3e}".format(np.dot(bfgs_memory.gradient[0]-bfgs_memory.gradient[1], bfgs_memory.deformation[0])))
    print(current_val)


