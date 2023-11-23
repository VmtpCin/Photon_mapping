import numpy as np
from Interseções import *
from math import sqrt, inf, sin, cos, acos, pi
from random import uniform

def normalize(v):
    norma = norm(v)
    if norma == 0:
        return v
    return v / norma


def russian_roulette(ponto, direcao, up, cor,  ir1, min, max, lista_objetos, lista_photons):

    inter, p, n, rr, ir2, tipo = intersecao(ponto, direcao, min, max, lista_objetos)
    n = np.array(n)
    if tipo == "esf" and cor == [255, 255, 255]:
        cor = [255, 0, 0]
    elif tipo == "plano" and cor == [255, 255, 255]:
        cor = [0, 255, 0]
    elif cor == [255, 255, 255]:
        cor = [0, 0, 255]
    if inter:
        dice = uniform(0, 1)
        if dice < rr[0]:  # difusão

            # gerando valores aleatorios entre 0 e 1
            sigma1 = uniform(0, 1)
            sigma2 = uniform(0, 1)

            # Gerar base ortonormal no local do ponto
            w = np.cross(n, up)
            v = np.cross(w, up)
            v = norm(v)
            u = np.cross(v, w)

            # Gerar o sample point
            teta = 1 / (cos(sqrt(sigma1)))  # angulo com a normal
            phi = 2 * pi * sigma2  # rotação ao redor da normal
            sp = (sin(teta) * cos(phi) * u) + sin(teta) * sin(phi) * v + cos(teta) * w
            novo_diretor = (sp[0] * u) + (sp[1] * v) + (sp[2] * w)

            # Testa photon refletido
            russian_roulette(p, novo_diretor, up, cor, ir1, 0.001, inf, lista_objetos, lista_photons)

        elif dice < rr[0] + rr[1]:  # reflexão
            novo_diretor = 2 * n * np.dot(n, -1 * direcao) - (-1 * direcao)
            russian_roulette(p, novo_diretor, up, cor, ir1, 0.001, inf, lista_objetos, lista_photons)

        elif dice < rr[0] + rr[1] + rr[2]:  # refração
            wo = normalize(-direcao)
            cos_teta_i = np.dot(wo, n)
            eta = ir2
            if cos_teta_i < 0:
                cos_teta_i = -cos_teta_i
                n = -n
                eta = 1/ir2
            temp = 1 - (1 - cos_teta_i * cos_teta_i) / (eta * eta)
            if temp > 0:
                cos_teta2 = sqrt(temp)
                novo_diretor = direcao / eta - (cos_teta2 - cos_teta_i / eta) * n
                russian_roulette(p, novo_diretor, up, cor, ir2, 0.001, inf, lista_objetos, lista_photons)

            # if ir1 == ir2:
            #     ir2 = 1
            #     n = -n
            # snell = ir2/ir1
            # cos_externo = np.dot(n, direcao)
            # x = ((1/snell**2) * (1 - cos_externo**2))
            # cos_interno = sqrt(1 - ((1/snell**2) * (1 - cos_externo**2)))
            # novo_diretor = ((1/snell) * (direcao)) - (cos_interno - (((1/snell) * cos_externo) * n))
            # teste = np.dot(n, novo_diretor)


        else:  # absorção
            lista_photons.append([p, cor])


def emitir_photons(ponto, num_photons, up, lista_photons, lista_objetos):

    for i in range(num_photons):

        while True:
            x = uniform(-1, 1)
            y = uniform(-1, 1)
            z = uniform(-1, 1)

            if not (x**2 + y**2 + z**2) > 1:
                break

        direcao = np.array([x, y, z])
        russian_roulette(np.array(ponto), direcao, np.array(up), [255, 255, 255], 1, 0.001, inf, lista_objetos, lista_photons)  #intersecao(np.array(ponto), direcao, 0.001, inf, lista_objetos)

        # if any(p):
        #     lista_photons.append(p)
