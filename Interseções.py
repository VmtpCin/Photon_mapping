# Definindo objetos e interseções com retas
import numpy as np
from numpy.linalg import norm
from math import sqrt, inf, sin, cos, acos
from random import uniform

def criar_esfera(centro, raio, lista_objetos):
    lista_objetos.append(["esf", np.array(centro), raio])


def criar_plano(ponto, normal, lista_objetos):
    lista_objetos.append(["pl", np.array(ponto), np.array(normal)])


def criar_triangulo(v1, v2, v3, lista_objetos):
    lista_objetos.append(["tri", np.array(v1), np.array(v2), np.array(v3)])


def criar_malha():
    num_vertices = int(input())
    num_faces = int(input())
    lista_vertices = []
    for i in range(num_vertices):
        vertice = np.array(input().split(), dtype='int')
        lista_vertices.append(vertice)
    for i in range(num_faces):
        v1 = lista_vertices[int(input())]
        v2 = lista_vertices[int(input())]
        v3 = lista_vertices[int(input())]
        rr = np.array(input().split(), dtype='int')
        criar_triangulo(v1, v2, v3)


def intersecao_esf(esfera, og, vetor_diretor):
    oc = og - esfera[1]
    a = np.dot(vetor_diretor, vetor_diretor)
    b = 2 * np.dot(oc, vetor_diretor)
    c = np.dot(oc, oc) - (esfera[2] ** 2)
    delta = b ** 2 - (4 * a * c)

    if delta < 0:
        return False, -1, [0, 0, 0]

    it1 = (-b + sqrt(delta)) / (2 * a)
    it2 = (-b - sqrt(delta)) / (2 * a)
    if it1 < it2:
        menor_t = it1
    else:
        menor_t = it2

    x = og[0] + menor_t * vetor_diretor[0]
    y = og[1] + menor_t * vetor_diretor[1]
    z = og[2] + menor_t * vetor_diretor[2]
    ponto = np.array([x, y, z])

    if it1 < 0 and it2 < 0:
        return False, -1, [0, 0, 0]
    else:
        return True, menor_t, ponto


def intersecao_pl(plano, og, vetor_diretor):
    holder = np.dot(plano[1], vetor_diretor)
    if holder == 0:
        return False, -1, [0, 0, 0]
    t = (np.dot(plano[1], plano[2]) - np.dot(plano[1], og)) / np.dot(plano[1], vetor_diretor)
    x = og[0] + t * vetor_diretor[0]
    y = og[1] + t * vetor_diretor[1]
    z = og[2] + t * vetor_diretor[2]
    ponto = np.array([x, y, z])
    return True, t, ponto


def intersecao_tri(tri, og, vetor_diretor):
    v0 = tri[1]
    v1 = tri[2]
    v2 = tri[3]

    a = v0[0] - v1[0]
    b = v0[0] - v2[0]
    c = vetor_diretor[0]
    d = v0[0] - og[0]
    e = v0[1] - v1[1]
    f = v0[1] - v2[1]
    g = vetor_diretor[1]
    h = v0[1] - og[1]
    i = v0[2] - v1[2]
    j = v0[2] - v2[2]
    k = vetor_diretor[2]
    l = v0[2] - og[2]

    m = (f * k) - (g * j)
    n = (h * k) - (g * l)
    p = (f * l) - (h * j)
    q = (g * i) - (e * k)
    s = (e * j) - (f * i)

    holder = ((a * m) + (b * q) + (c * s))
    if holder == 0:
        inv_denom = 0
    else:
        inv_denom = 1/((a * m) + (b * q) + (c * s))

    e1 = (d * m) - (b * n) - (c * p)
    beta = e1 * inv_denom
    if beta < 0:
        return False, -1, [0, 0, 0]

    r = (e * l) - (h * i)
    e2 = (a * n) + (d * q) + (c * r)
    gamma = e2 * inv_denom
    if gamma < 0:
        return False, -1, [0, 0, 0]

    if beta + gamma > 1:
        return False, -1, [0, 0, 0]

    e3 = (a * p) - (b * r) + (d * s)
    t = e3 * inv_denom

    if t < 1:
        return False, -1, [0, 0, 0]

    x = og[0] + t * vetor_diretor[0]
    y = og[1] + t * vetor_diretor[1]
    z = og[2] + t * vetor_diretor[2]
    ponto = np.array([x, y, z])
    return True, t, ponto


def intersecao(ponto, vetor_diretor, min, max, lista_objetos):
    intersecta = False
    menor_t = 10000
    p = [0, 0, 0]
    cor = [0, 0, 0]
    for obj in lista_objetos:
        tipo_obj = obj[0]
        if tipo_obj == "esf":
            intersecta, t, ponto_temp = intersecao_esf(obj, ponto, vetor_diretor)
            if intersecta is True and t < menor_t and min < t < max:
                menor_t = t
                p = ponto_temp
                cor = [255, 0, 0]
        elif tipo_obj == "pl":
            intersecta, t,  ponto_temp = intersecao_pl(obj, ponto, vetor_diretor)
            if intersecta is True and t < menor_t and min < t < max:
                menor_t = t
                p = ponto_temp
                cor = [0, 255, 0]
        else:
            intersecta, t, ponto_temp = intersecao_tri(obj, ponto, vetor_diretor)
            if intersecta is True and t < menor_t and min < t < max:
                menor_t = t
                p = ponto_temp
                cor = [0, 0, 255]
    return p, cor

