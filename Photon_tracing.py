from numpy import cross
from Interseções import *
from numpy.linalg import norm
from math import sqrt, inf, sin, cos, acos, pi
from random import uniform


def russian_roulette(ponto, direcao, up,  min, max, lista_objetos, lista_photons):

    inter, p, n, rr = intersecao(ponto, direcao, min, max, lista_objetos)
    n = np.array(n)
    if inter:
        dice = uniform(0, 1)
        if dice < rr[0]:

            # gerando valores aleatorios entre 0 e 1
            sigma1 = uniform(0, 1)
            sigma2 = uniform(0, 1)
            a = cross(n, n)
            # Gerar base ortonormal no local do ponto



            teta = 1/(cos(sqrt(sigma1)))  # angulo com a normal
            phi = 2 * pi * sigma2  # rotação ao redor da normal

            # difusão

        elif dice < rr[0] + rr[1]:
                a = 2
                # reflexão

        elif dice < rr[0] + rr[1] + rr[2]:
            a = 2
            #refrata

        else:
            lista_photons.append(p)


def emitir_photons(ponto, num_photons, up, lista_photons, lista_objetos):

    for i in range(num_photons):

        while True:
            x = uniform(-1, 1)
            y = uniform(-1, 1)
            z = uniform(-1, 1)

            if not (x**2 + y**2 + z**2) > 1:
                break

        direcao = np.array([x, y, z])
        russian_roulette(np.array(ponto), direcao, np.array(up), 0.001, inf, lista_objetos, lista_photons)  #intersecao(np.array(ponto), direcao, 0.001, inf, lista_objetos)

        # if any(p):
        #     lista_photons.append(p)
