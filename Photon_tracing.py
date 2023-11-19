# Definindo objetos e interseções com retas
import numpy as np
from Interseções import *
from numpy.linalg import norm
from math import sqrt, inf, sin, cos, acos
from random import uniform


def russian_roulete(ponto, direcao, min, max, lista_photons):
    inter, p, c, rr = intersecao(ponto, direcao, min, max)
    if inter:
        dice = uniform(0, 1)
        if dice < rr[0]:
            a = 2
            # difusão

        elif dice < rr[0] + rr[1]:
                a = 2
                # relexão

        elif dice < rr[0] + rr[1] + rr[2]:
            a = 2
            #refrata

        else:
            lista_photons.append()


def emitir_photons(ponto, num_photons, lista_photons, lista_objetos):
    for i in range(num_photons):
        while True:
            x = uniform(-1, 1)
            y = uniform(-1, 1)
            z = uniform(-1, 1)
            if not (x**2 + y**2 + z**2) > 1:
                break
        direcao = np.array([x, y, z])
        p, c = intersecao(np.array(ponto), direcao, 0.001, inf, lista_objetos)
        if any(p):
            lista_photons.append(p)
