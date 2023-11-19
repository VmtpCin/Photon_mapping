# Definindo objetos e interseções com retas
import numpy as np
from Interseções import *
from numpy.linalg import norm
import cv2 as cv
from math import sqrt, inf, sin, cos, acos
from random import uniform


def normalize(v):
    norma = norm(v)
    if norma == 0:
        return v
    return v / norma


def vetor_deslocamento(x, y, k, m, b, u):
    dl = ((2 * x) / (k - 1)) * b
    dv = ((2 * y) / (m - 1)) * u
    return dl, dv


def raycast(lista_objetos):
    camera = np.array([3, 0, 0])
    alvo = np.array([1, 0, 0])
    vetor_up = np.array([0, 1, 0])
    distancia = 1
    hres = 500
    vres = 500
    tamx = 0.5
    tamy = 0.5
    grid = np.zeros((vres, hres, 3), dtype=np.uint8)
    oa = normalize(alvo - camera)
    b = normalize(np.cross(oa, vetor_up))
    up = normalize(-1 * np.cross(b, oa))

    desl_l, desl_v = vetor_deslocamento(tamx, tamy, hres, vres, b, up)
    aj = oa * distancia
    ajj = tamx * b
    ajk = tamy * up
    vet_incial = oa * distancia - tamx * b - tamy * up

    vet_at = vet_incial
    for i in range(vres):
        if i > 0:
            vet_at[2] = vet_incial[2]
            vet_at = vet_at + desl_v
        for k in range(hres):
            if k > 0:
                vet_at = vet_at + desl_l
            p, grid[i, k] = intersecao(camera, vet_at, distancia, inf, lista_objetos)

    cv.imshow('i', grid)
    cv.waitKey(0)
    cv.destroyWindow('i')



def starfield_projection(photons):
    camera = np.array([3, 0, 0])
    alvo = np.array([1, 0, 0])
    vetor_up = np.array([0, 1, 0])
    distancia = 1
    hres = 500
    vres = 500
    tamx = 0.5
    tamy = 0.5
    grid = np.zeros((vres, hres, 3), dtype=np.uint8)
    oa = normalize(alvo - camera)
    b = normalize(np.cross(oa, vetor_up))
    up = normalize(-1 * np.cross(b, oa))

    desl_l, desl_v = vetor_deslocamento(tamx, tamy, hres, vres, b, up)
    vet_incial = oa * distancia - tamx * b - tamy * up
    intersecta_i, t_i, p_i = intersecao_pl(["pl", np.array(alvo), np.array(camera - alvo), [0, 0, 0]], camera, vet_incial)

    for photon in photons:
        vetor_proj = camera - photon
        intersecta, t, p = intersecao_pl(["pl", np.array(alvo), np.array(camera - alvo), [0, 0, 0]], photon, vetor_proj)
        if intersecta is True:
            vetor = p - p_i
            v_desl_l = np.divide(vetor, desl_l, out=np.zeros_like(vetor), where=desl_l!=0)
            if min(v_desl_l) == 0:
                num_des_l = max(v_desl_l)
            else:
                num_des_l = min(v_desl_l)
            v_desl_v = np.divide(vetor, desl_v, out=np.zeros_like(vetor), where=desl_v!=0)
            if min(v_desl_v) == 0:
                num_des_v = max(v_desl_v)
            else:
                num_des_v = min(v_desl_v)

            posicao_grid_y = round(num_des_l)
            posicao_grid_x = round(num_des_v)
            if posicao_grid_x <= hres and posicao_grid_y <= vres:
                grid[posicao_grid_x, posicao_grid_y] = [255, 255, 255]

    cv.imshow('i', grid)
    cv.waitKey(0)
    cv.destroyWindow('i')