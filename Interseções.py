# Definindo objetos e interseções com retas
import numpy as np
from numpy.linalg import norm
from math import sqrt, inf, sin, cos, acos


def criar_esfera(centro, raio, rr, ir, lista_objetos):
    lista_objetos.append(["esf", np.array(centro), raio, np.array(rr), ir])


def criar_plano(ponto, normal, rr, ir, lista_objetos):
    lista_objetos.append(["pl", np.array(ponto), np.array(normal), np.array(rr), ir])


def criar_triangulo(v1, v2, v3, rr, ir, lista_objetos):
    lista_objetos.append(["tri", np.array(v1), np.array(v2), np.array(v3), np.array(rr), ir])


def criar_malha(lista_objetos):
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
        ir = float(input())
        criar_triangulo(v1, v2, v3, rr, ir, lista_objetos)


def intersecao_esf(esfera, og, vetor_diretor):
    oc = og - esfera[1]
    a = np.dot(vetor_diretor, vetor_diretor)
    b = 2 * np.dot(oc, vetor_diretor)
    c = np.dot(oc, oc) - (esfera[2] ** 2)
    delta = (b ** 2) - (4 * a * c)

    if delta < 0:
        return False, -1, [0, 0, 0], [0, 0, 0], [0, 0, 0], 0

    it1 = (-b + sqrt(delta)) / (2 * a)
    it2 = (-b - sqrt(delta)) / (2 * a)


    if it1 < 0.01:
        it1 = inf
    if it2 < 0.01:
        it2 = inf

    if it1 < 0 and it2 < 0:
        return False, -1, [0, 0, 0], [0, 0, 0], [0, 0, 0], 0

    if it1 == inf and it2 == inf:
        return False, -1, [0, 0, 0], [0, 0, 0], [0, 0, 0], 0

    if it1 < it2:
        menor_t = it1
    else:
        menor_t = it2

    x = og[0] + menor_t * vetor_diretor[0]
    y = og[1] + menor_t * vetor_diretor[1]
    z = og[2] + menor_t * vetor_diretor[2]

    ponto = np.array([x, y, z])
    normal = np.array(ponto - esfera[1])
    rr = np.array(esfera[3])
    ir = esfera[4]


    return True, menor_t, ponto, normal, rr, ir


def intersecao_pl(plano, og, vetor_diretor):
    holder = np.dot(plano[2], vetor_diretor)

    if holder == 0:
        return False, -1, [0, 0, 0], [0, 0, 0], [0, 0, 0], 0

    t = (np.dot(plano[2], plano[1]) - np.dot(plano[2], og)) / np.dot(plano[2], vetor_diretor)

    x = og[0] + t * vetor_diretor[0]
    y = og[1] + t * vetor_diretor[1]
    z = og[2] + t * vetor_diretor[2]

    ponto = np.array([x, y, z])
    normal = np.array(plano[2])
    rr = np.array(plano[3])
    ir = plano[4]

    return True, t, ponto, normal, rr, ir


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
        return False, -1, [0, 0, 0], [0, 0, 0], [0, 0, 0], 0

    r = (e * l) - (h * i)
    e2 = (a * n) + (d * q) + (c * r)
    gamma = e2 * inv_denom
    if gamma < 0:
        return False, -1, [0, 0, 0], [0, 0, 0], [0, 0, 0], 0

    if beta + gamma > 1:
        return False, -1, [0, 0, 0], [0, 0, 0], [0, 0, 0], 0

    e3 = (a * p) - (b * r) + (d * s)
    t = e3 * inv_denom

    if t < 1:
        return False, -1, [0, 0, 0], [0, 0, 0], [0, 0, 0], 0

    x = og[0] + t * vetor_diretor[0]
    y = og[1] + t * vetor_diretor[1]
    z = og[2] + t * vetor_diretor[2]

    ponto = np.array([x, y, z])
    normal = np.array(np.cross((v1 - v0), (v2 - v0)))
    normal = norm(normal)
    rr = np.array(tri[4])
    ir = tri[5]

    return True, t, ponto, normal, rr, ir


def intersecao_rt(ponto, vetor_diretor, min, max, lista_objetos):
    intersecta = False
    menor_t = 10000
    p = [0, 0, 0]
    cor = [0, 0, 0]

    for obj in lista_objetos:

        tipo_obj = obj[0]

        if tipo_obj == "esf":
            intersecta, t, ponto_temp, n, rr, ir = intersecao_esf(obj, ponto, vetor_diretor)
            if intersecta is True and t < menor_t and min < t < max:
                menor_t = t
                cor = [255, 0, 0]

        elif tipo_obj == "pl":
            intersecta, t,  ponto_temp, n, rr, ir = intersecao_pl(obj, ponto, vetor_diretor)
            if intersecta is True and t < menor_t and min < t < max:
                menor_t = t
                cor = [0, 255, 0]

        else:
            intersecta, t, ponto_temp, n, rr, ir = intersecao_tri(obj, ponto, vetor_diretor)
            if intersecta is True and t < menor_t and min < t < max:
                menor_t = t
                cor = [0, 0, 255]
    return cor


def intersecao(ponto, vetor_diretor, min, max, lista_objetos):
    inters = False
    menor_t = 10000
    ir2 = 0
    p = [0, 0, 0]
    normal = [0, 0, 0]
    rr = [0, 0, 0]
    tipo = ""

    for obj in lista_objetos:

        tipo_obj = obj[0]

        if tipo_obj == "esf":
            inter, t, ponto_temp, normal_temp, rr_temp, ir_temp = intersecao_esf(obj, ponto, vetor_diretor)
            if inter is True and t < menor_t and min < t < max:
                menor_t = t
                p = ponto_temp
                normal = normal_temp
                ir2 = ir_temp
                rr = rr_temp
                tipo = "esf"
                inters = True

        elif tipo_obj == "pl":
            inter, t,  ponto_temp, normal_temp, rr_temp, ir_temp = intersecao_pl(obj, ponto, vetor_diretor)
            if inter is True and t < menor_t and min < t < max:
                menor_t = t
                p = ponto_temp
                ir2 = ir_temp
                normal = normal_temp
                rr = rr_temp
                tipo = "plano"
                inters = True

        else:
            inter, t, ponto_temp, normal_temp, rr_temp, ir_temp = intersecao_tri(obj, ponto, vetor_diretor)
            if inter is True and t < menor_t and min < t < max:
                menor_t = t
                p = ponto_temp
                ir2 = ir_temp
                normal = normal_temp
                tipo = "tri"
                rr = rr_temp
                inters = True

    return inters, p, normal, rr, ir2, tipo

