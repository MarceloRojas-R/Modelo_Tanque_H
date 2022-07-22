import numpy as np
import matplotlib.pyplot as plt


def Tanque(t, X):
    """
    :param t:
    :param X:
    :return: T = [ Tg  Pg  w  Perfil_de_temperaturas_pared]
    Modelo de tanque de alta presion con entrada de 480 bar
    """

    T = np.ones(64)

    pf = 480e5           # [Pa] presion final
    po = 6e5             # [Pa] presion inicial
    t_fill = 5 * 60      # [s]
    w = 0.3 / 60         # [kg/s]
    t_inf = 273.15 + 25  # [K]

    # propiedades fisicas del H2

    pm_g = 2.016               # [kg / kmol]
    mass_rho_g = 0.08375       # [kg / m3]
    k_g = 0.1825               # [W / (m * K)]
    rho_g = mass_rho_g / pm_g  # [kmol / m3]
    cp_g = 14.29e3             # [J/ (kg * K)]
    beta_g = 1/330             # [1 / K]

    # coeficientes de tranferencia de calor

    hi = 250  # [W / (m2 * K)]
    ho = 30   # [W / (m2 * K)]

    # dimensiones tanque

    v = 300e-3     # [m3]
    D = 0.28       # [m]
    si = 2 * 0.66  # [m2]

    # LINER propiedades termofisicas

    t1 = 5e-3    # [m] espesor
    k1 = 1.17    # [W / (m * K)]
    c1 = 1578    # [J / (kg * K)]
    Rho1 = 1286  # [kg / m3]

    # discretizacion EDP
    N = 20                       # numero de elementos
    delta_X1 = t1 / N
    alpha_1 = k1 / (Rho1 * c1)
    beta_1 = hi / (Rho1 * c1)

    # CRFP
    # Geometria y Propiedades Termofísicas
    t2 = 30e-3 # [m] espesor
    k2 = 0.66 # [W / (m * K)]
    c2 = 1075 # [J / (kg * K)]
    Rho2 = 1375 # [kg / m3]

    #discretacion EDP
    M = 40
    delta_X2 = t2 / M
    alpha_2 = k2 / (Rho2 * c2)
    beta_2 = ho / (Rho2 * c2)

    t_ref = 273.15              # [K]
    t_e = 10 + 273.15           # [K]
    he = cp_g * (t_e - t_ref)
    h = cp_g * (X[1] - t_ref)

    if t <= t_fill:
        T[1] = (pf-po)/t_fill
        T[2] = w
    else:
        T[1] = 0
        T[2] = 0
    T[0] = (1 / (X[2] * cp_g)) * (v * beta_g * X[0] * T[1] + hi * si * (X[3]-X[0])) + w * (he-h)

    T[3] = (2 * alpha_1 / delta_X1 ** 2) * (X[4]-X[3]) + (2 * beta_1 / delta_X1) * (X[0]-X[3])

    for i in range(5,N):
        T[i] = (2 * alpha_1 / delta_X1 ** 2)
