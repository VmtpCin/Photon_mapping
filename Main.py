from Rendering import *
from Photon_tracing import *

lista_objetos = []
lista_photons = []

criar_esfera([4, 0, 0], 0.25, [0, 0, 0], 1, lista_objetos)
criar_plano([0, -1, 0], [0, 1, 0], [0, 0, 0], 0, lista_objetos)

# criar_triangulo([0, 1, 1], [0, 1, -1], [0, 0, 0], [0, 0, 0], 0, lista_objetos)

raycast(lista_objetos)
emitir_photons([-1, 1, 0], 1000000, [0, 1, 0], lista_photons, lista_objetos)
starfield_projection(lista_photons)
