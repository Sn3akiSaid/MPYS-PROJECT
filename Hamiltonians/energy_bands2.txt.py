# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 23:43:32 2024

@author: charl
"""
import numpy as np
import matplotlib.pyplot as plt


BANDS = 18
HOPPING_SITES = 1155
L = [0.5, 0, 0.5]
A = [0, 0, 0.5]
H = [1/3, 1/3, 0.5]


def apply_weights(data, weights):
    weights = np.ndarray.flatten(weights)
    counter = 0
    for i in range(len(data[:, 0])):

        if i == 0:
            a = 1
        elif data[i, 2] != data[i-1, 2]:
            counter += 1
        else:
            b = 1
        data[i, 5] = data[i, 5]/weights[counter]
        data[i, 6] = data[i, 6]/weights[counter]

    return data


def matrix_gen(data):
    matrix = np.empty((BANDS, BANDS))
    counter = 0

    for i in range(len(data[:, 0])):
        counter += 1
        if data[i, 3] == 1:
            counter_i = 0
        else:
            counter_i += 1
        if data[i, 4] == 1:
            counter_j = 0
        if data[i, 3] == 1 and data[i, 4] == 1:
            temp = np.empty((BANDS, BANDS), dtype=np.ndarray)

        temp[counter_i, counter_j] = np.array((data[i, 5], data[i, 6]))

        if counter_i == BANDS-1:
            counter_j += 1
        if data[i, 3] == 18 and data[i, 4] == 18:
            matrix = np.vstack((matrix, temp))

    matrix = np.delete(matrix, range(BANDS), axis=0)
    return matrix


def k_gen():
    samples = 10
    l_a = np.zeros(samples, dtype=np.ndarray)
    templa = np.array((np.linspace(L[0], A[0], samples), np.linspace(L[1], A[1], samples),
                       np.linspace(L[2], A[2], samples)))

    a_h = np.zeros(samples, dtype=np.ndarray)
    tempah = np.array((np.linspace(A[0], H[0], samples), np.linspace(A[1], H[1], samples),
                       np.linspace(A[2], H[2], samples)))

    h_l = np.zeros(samples, dtype=np.ndarray)
    temphl = np.array((np.linspace(H[0], L[0], samples), np.linspace(H[1], L[1], samples),
                       np.linspace(H[2], L[2], samples)))

    for i in range(samples):
        l_a[i] = np.array((templa[0, i], templa[1, i], templa[2, i]))
        a_h[i] = np.array((tempah[0, i], tempah[1, i], tempah[2, i]))
        h_l[i] = np.array((temphl[0, i], temphl[1, i], temphl[2, i]))

    k_values = np.hstack((l_a, a_h))
    k_values = np.hstack((k_values, h_l))
    return (k_values)


def r_val_calc(data):
    r_vals = np.zeros(HOPPING_SITES, dtype=np.ndarray)
    counter = 0
    print(data[5, 2])
    print(data[4, 2])
    for i in range(len(data[:, 0])):
        if i == 0:
            a = 1
        elif data[i, 2] != data[i-1, 2]:
            r_vals[counter] = np.array(
                (data[i-1, 0], data[i-1, 1], data[i-1, 2]))
            counter += 1
        else:
            b = 1


def main():
    weights = np.genfromtxt('trivial_weightings.txt')
    top_data = np.genfromtxt('topological_bulk.txt')
    triv_data = np.genfromtxt('trivial_bulk.txt')

    top_weighted = apply_weights(top_data, weights)
    triv_weighted = apply_weights(triv_data, weights)

    top_matrix = matrix_gen(top_weighted)
    triv_matrix = matrix_gen(triv_weighted)
    k_values = k_gen()
    r_values = r_val_calc(top_data)

    # hk_values = np.empty(BANDS)
    # counter = 0
    # for i in k_values:
    #     temp =
    return (r_values)


a = main()
b = k_gen()
