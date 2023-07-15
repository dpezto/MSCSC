"""
-------------------------------------------------------------------------------------
   Autor: Dai López Jacinto
-------------------------------------------------------------------------------------
"""

import os
from numpy import savetxt, column_stack

# -----------------------------------------------------------------------------------
# DEFINICIONES DE FUNCIONES
# -----------------------------------------------------------------------------------


def get_k(basis_i):
    """
    Toma elemento [i] de la base y regresa el índice k
    para generar el siguiente elemento
    """
    from numpy import zeros
    M = basis_i.size
    for i in range(M - 1):
        if basis_i[i] != 0:
            if all(basis_i[i + 1:M - 1] == zeros(M - 1 - i - 1, int)):
                k = i + 1
    return k


def next_element(basis_i, N, M):
    """
    Dado el elemento [i] de la base de Fock genera el elemento [i+1]
    Parámetros:
        basis_i: elemento [i] de la base
        N: partículas
        M: sitios
    Regresa:
        basis_next: siguente elemento de la base
    """
    k = get_k(basis_i)
    # ni_ = ni para 1 <= i <= k-1
    basis_next = basis_i.copy()
    # nk_ = nk - 1
    basis_next[k - 1] = basis_next[k - 1] - 1
    # nk+1_ = N - sum_i=1^k(ni_)
    basis_next[k] = N - sum(basis_next[i] for i in range(k))
    # ni_ = 0 para i >= k+2
    for i in range(k + 1, M):
        basis_next[i] = 0
    return basis_next


def fock_basis(N, M, D):
    """
    Genera la base de fock
    Parámetros:
        N: partículas
        M: sitios
        D: dimensión
    Regresa:
        basis: base completa en orden lexicográfico
    """
    from numpy import zeros
    basis = zeros((D, M), int)
    basis[0, 0] = N  # primera entrada de la base
    # Construyo el resto de la base
    for i in range(D - 1):
        basis[i + 1] = next_element(basis[i], N, M)
    return basis


def tag(A_v, M):
    """
    Función que genera el tag del vector base A_v
    Parámetros: A_v (vector base), M (sitios)
    """
    from numpy import sqrt
    return sum(sqrt(100 * i + 3) * A_v[i] for i in range(M))


def partition(arr, lo, hi):
    """
    Partición para ordenar usando quickSort
    """
    i = (lo - 1)  # índice del elemento más pequeño
    pivote = arr[hi]
    for j in range(lo, hi):
        if arr[j] <= pivote:  # Si elemento es <= que el pivote
            i = i + 1  # => incrementa el índice del menor elemento
            arr[i], arr[j] = arr[j], arr[i]
    arr[i + 1], arr[hi] = arr[hi], arr[i + 1]
    return (i + 1)


def quickSort(arr, lo, hi):
    """
    Algoritmo quickSort
    """
    if arr.size == 1:
        return arr
    if lo < hi:
        # índice de partición
        pi = partition(arr, lo, hi)
        # Ordena por separado los elementos antes y después de la partición
        quickSort(arr, lo, pi - 1)
        quickSort(arr, pi + 1, hi)


def sortit(arr):
    """
    Ordena array
    """
    quickSort(arr, 0, arr.size - 1)
    # Regresa la lista completa ordenada
    return arr


def binary_search(arr, x):
    """
    Algoritmo de bisección
    Regresa el índice 'i' de un array 'arr' ordenado donde el valor 'x'
    debe ser insertado para mantener el orden.
    """
    lo, hi = 0, arr.size - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if arr[mid] < x:  # Si x es mayor, ignora la mitad inferior
            lo = mid + 1
        elif arr[mid] > x:  # Si x es menor, ignora la mitad superior
            hi = mid - 1
        else:  # x está en medio
            return mid
    # El elemento no está presente
    return -1


def expval(Op, ψ):
    """
    Calcula el valor esperado de un operador
    Parámetros:
        Op: Operador
        ψ: en este caso el estado base
    Regresa:
        expval: el valor esperado
    """
    # como kets vec col => ψ.T
    return ψ.conjugate() @ Op @ ψ.T


def normalize(v):
    """
    Normaliza vector
    """
    from numpy import sqrt
    return v / sqrt((v**2).sum())


def b_cd_p(vec, M, N):
    """
    periodicas
    """
    from numpy import sqrt, empty
    b = empty(M, int)  # el ket de descenso
    b[M - 1] = vec[M - 1] + 1  # aplica operador creación
    b[0] = vec[0] - 1  # aplica op aniquilación
    # cte del ket de ascenso y descenso
    cte_c = sqrt(vec[M - 1] + 1)
    cte_d = sqrt(vec[0])
    cte_tot_p = cte_c * cte_d

    for i in range(M):
        if i != 0 and i != M - 1:
            b[i] = vec[i]

    if b[0] < 0:  # Si < 0 lo descarta
        cte_d = 0
    elif b[M - 1] > N:  # Si > N lo descarta
        cte_c = 0
    if cte_c == 0 or cte_d == 0:
        b = 0

    return vec, b, cte_tot_p


def b_cd(vec, indice, M, N):
    from numpy import sqrt, empty
    b = empty(M, int)  # el ket de descenso
    b[indice] = vec[indice] + 1  # aplica operador creacion
    b[indice + 1] = vec[indice + 1] - 1  # aplica op aniquilación
    # cte del ket de ascenso y descenso
    cte_c = sqrt(vec[indice] + 1)
    cte_d = sqrt(vec[indice + 1])
    cte_tot = cte_c * cte_d

    for i in range(0, M):
        if i != indice and i != indice + 1:
            b[i] = vec[i]
    if b[indice + 1] < 0:  # Si < 0 lo descarta
        cte_d = 0
    elif b[indice] > N:  # Si > N lo descartas
        cte_c = 0
    if cte_c == 0 or cte_d == 0:
        b = 0

    return vec, b, cte_tot


def hamiltonian(N, M):
    """
    Construye el hamiltoniano (H)
    Parámetros: N, M
    Regresa: H_kin, H_int, A (base de Fock), Ts (tags ordenados) e
    idx (índices que ordenan los tags)
    """
    from numpy import argsort, empty
    from math import factorial
    from scipy.sparse import coo_matrix, diags

    # coo_matrix faster for constructing large matrices
    # lil_matrix faster for changes to sparsity structure
    # csr_matrix faster for arithmetics operations

    D = factorial(N + M - 1) // (factorial(N) * factorial(M - 1))  # Dimensión
    A = fock_basis(N, M, D)  # Genero base de Fock
    T = empty(D)

    for v in range(D):  # Genero Tags
        T[v] = tag(A[v], M)
    idx = argsort(T)  # índices que ordenan T
    Ts = sortit(T)  # Ordeno los Tags

    # Diagonal/potencial -----------------------------------------------------
    inter, interaccion = empty([M, D]), empty(D)
    for i in range(M):
        for j in range(D):
            inter[i, j] = A[idx[j], i]**2 - A[idx[j], i]
    for i in range(D):
        interaccion[i] = inter[:, i].sum() / 2
    H_int = diags(interaccion).tocsr()

    # Off-diagonal/cinética --------------------------------------------------
    # TODA ESTA PARTE PUEDE MEJORAR MUCHÍSIMO
    H_kin = coo_matrix((D, D)).tolil()
    bra, ket = empty([D, M, M]), empty([D, M, M])
    ctes = empty([D, M])
    for i in range(D):
        for j in range(M-1):
            bra[i, j, :], ket[i, j, :], ctes[i, j] = b_cd(A[idx[i]], j, M, N)
            bra[i, M-1, :], ket[i, M-1, :], ctes[i, M-1] = b_cd_p(A[idx[i]], M, N)
    # Tags de los bras y kets
    Tbra, Tket = empty([D, M]), empty([D, M])
    for i in range(D):
        for j in range(M):
            Tbra[i, j] = tag(bra[i, j, :], M)
            Tket[i, j] = tag(ket[i, j, :], M)
    bra_i, ket_i = coo_matrix((D, M), int).tolil(), coo_matrix((D, M), int).tolil()
    for i in range(D):
        for j in range(M):
            # Matriz de bras y kets para H_kin
            if Tbra[i, j] != 0:
                bra_i[i, j] = binary_search(Ts, Tbra[i, j])
            else:
                bra_i[i, j] = -1
            if Tket[i, j] != 0:
                ket_i[i, j] = binary_search(Ts, Tket[i, j])
            else:
                ket_i[i, j] = -1
            # llena H_kin
            if ket_i[i, j] >= 0:  # Descarta tags = -1
                H_kin[bra_i[i, j], ket_i[i, j]] = ctes[i, j]
                # Llena complementarios
                H_kin[ket_i[i, j], bra_i[i, j]] = H_kin[bra_i[i, j], ket_i[i, j]]
    H_kin = H_kin.tocsr()

    return H_kin, H_int, A, Ts, idx


def hamiltonian_tU(H_kin, H_int, U=1, t=1):
    """
    Regresa el hamiltoniano H completo para los valores de U y t de entrada
    eq 1 Raventós
    """
    return -t * H_kin + U * H_int


def b_c(vec, indice, M, N):
    """
    b_create
    """
    from numpy import sqrt, empty
    b = empty(M, int)  # el ket de descenso
    b[indice] = vec[indice] + 1  # aplico operador creacion
    cte_c = sqrt(vec[indice] + 1)  # cte creación

    for i in range(0, M):
        if i != indice:
            b[i] = vec[i]

    if b[indice] > N:  # Si > N lo descarta
        cte_c = 0
    if cte_c == 0:
        b = 0
    if vec.sum() == 0:
        b, cte_c = 0, 0

    return b, cte_c


def b_d(vec, indice, M, N):
    """
    b_destroy
    """
    from numpy import sqrt, empty
    b = empty(M, int)  # el ket de descenso
    b[indice] = vec[indice] - 1  # aplico operador de aniquilacion
    cte_d = sqrt(vec[indice])  # cte aniquilación

    for i in range(0, M):
        if i != indice:
            b[i] = vec[i]

    if b[indice] < 0:  # Si < 0 lo descarta
        cte_d = 0
    if cte_d == 0:
        b = 0

    return vec, b, cte_d


def reduced_matrix(M, ψg, A, Ts, index):
    """
    Para la fracción de condensado

    ESTA BASURA PUEDE MEJORAR MUCHÍSIMO, IGUAL QUE LOS OPERADORES C/D
    """
    # Como es muy pequeña no vale la pena construir como sparse
    from numpy import zeros

    D = A.shape[0]
    bra, ket = zeros([D, M, M], int), zeros([D, M, M], int)
    ctes_d = zeros([D, M])  # para el de descenso
    ket_2 = zeros([D, M, M, M], int)
    # Ascenso (M casos)
    ctes_c, ctes_mr = zeros([D, M, M]), zeros([D, M, M])

    for i in range(D):
        for j in range(M):
            # primero se aplica el de aniquilacion
            bra[i, j, :], ket[i, j, :], ctes_d[i, j] = b_d(
                A[index[i]], j, M, M)
            for k in range(M):
                ket_2[i, j, k, :], ctes_c[i, j, k] = b_c(ket[i, j, :], k, M, M)
                ctes_mr[i, j, k] = ctes_c[i, j, k] * ctes_d[i, j]

    # Tags de los bra
    Tbra = zeros([D, M])
    for i in range(D):
        for j in range(M):
            Tbra[i, j] = tag(bra[i, j, :], M)
    # Posición de los bras
    bra_i = zeros([D, M], int)
    for i in range(D):
        for j in range(M):
            if Tbra[i, j] != 0:
                bra_i[i, j] = binary_search(Ts, Tbra[i, j])
            else:
                bra_i[i, j] = -1
    # Tags de los ket
    Tket = zeros([M, M, D])
    for i in range(M):
        for j in range(M):
            for k in range(D):
                Tket[i, j, k] = tag(ket_2[k, j, i, :], M)
    # Posición de los kets
    ket_i = zeros([M, M, D], int)
    for i in range(M):
        for j in range(M):
            for k in range(D):
                if Tket[i, j, k] != 0:
                    ket_i[i, j, k] = binary_search(Ts, Tket[i, j, k])
                else:
                    ket_i[i, j, k] = -1
    # LLena matriz reducida
    Mr = zeros([M, M])
    for i in range(M):
        for j in range(M):
            for k in range(D):
                if ket_i[i, j, k] >= 0:  # Descarto los -1
                    Mr[j, i] += ψg[bra_i[k, j]] * \
                        ctes_mr[k, j, i] * ψg[ket_i[i, j, k]]
    return Mr


def quantities(M, U, t, Ham, E0p, E0m):
    """
    Regresa fraciones de condensado (fc), fluctuaciones (Δn2),
    energía cinética (kin) y potencial químico (μp y μm)
    """
    from numpy import diagflat, empty
    from scipy.sparse.linalg import eigsh

    A = Ham[2]  # Base de Fock
    Ts = Ham[3]  # Tags ordenados
    idx = Ham[4]  # Índices que odenan los tags
    # Para potencial químico
    H = hamiltonian_tU(Ham[0], Ham[1], U, t)  # Para fc, Δn2 y kin
    Hp = hamiltonian_tU(E0p[0], E0p[1], U, t)
    Hm = hamiltonian_tU(E0m[0], E0m[1], U, t)

    # Eigenval y eigenvec
    λ, ψ = eigsh(H, 1, which='SA', tol=1e-5, mode='buckling')  # eigenval menor
    λp, ψp = eigsh(Hp, 1, which='SA', tol=1e-5, mode='buckling')
    λm, ψm = eigsh(Hm, 1, which='SA', tol=1e-5, mode='buckling')

    λg = λ.min()  # min energía
    ψg = normalize(ψ[:, binary_search(λ, λg)])  # edo. base normalizado
    λgp = λp.min()
    λgm = λm.min()

    # Potencial químico -------------------------------------------------------
    μp = λgp - λg  # eq 24 Raventós
    μm = λg - λgm  # eq 25 Raventós
    # Fracción de condensado --------------------------------------------------
    Mr = reduced_matrix(M, ψg, A, Ts, idx)
    λ2, ψ2 = eigsh(Mr, 1, which='LA', tol=1e-5)  # eigenval mayor
    λmax = λ2.max()  # max eigenval
    fc = λmax/M  # ya que N = M
    # Energía cinética --------------------------------------------------------
    np, ns = empty([M, A.shape[0]]), empty([M, A.shape[0]])  # [M,D]
    for i in range(M):
        for j in range(A.shape[0]):
            np[i, j] = A[idx[j], i] * ψg[j]**2
            ns[i, j] = A[idx[j], i]**2 * ψg[j]**2
    kin = λg - U/2 * (ns.sum() - np.sum())
    # Fluctuaciones -----------------------------------------------------------
    n = [diagflat(A[:, i]) for i in range(M)]  # Operador de número por sitio
    Δn2 = empty(M)
    for i in range(M):
        Δn2[i] = expval(n[i]**2, ψg) - expval(n[i], ψg)**2

    return fc, Δn2[0], kin, μp, μm


def graphs(M, s, H, E0p, E0m, mode):
    """
    Datos para las gráficas
    Parámetros: M (partículas/sitios), s (samples), H/E0p/E0m (Hamiltonianos)
    mode (ttilde ó utilde)
    Regresa: column_stack con los datos de las gráficas
    """
    from numpy import linspace, column_stack
    from numpy import empty as em

    FC, Δn, K, Mp, Mm = em(s), em(s), em(s), em(s), em(s)
    if mode == "ttilde":
        U = 1.  # H_int, Interacción, E. interacción
        t = linspace(0, 0.5, s)
        for i in range(len(t)):
            FC[i], Δn[i], K[i], Mp[i], Mm[i] = quantities(M, U, t[i], H, E0p, E0m)
        dat_ttilde = column_stack((
            t/U,  # ttilde
            FC, Δn, K,  # frac condensado, fluctuaciones, kinetic
            Mp, Mm  # chem pot + y -
        ))
        return dat_ttilde
    elif mode == "utilde":
        t = 1.  # H_kin, Tight binding, E. cinética
        U = linspace(0, 20, s)
        for i in range(len(U)):
            FC[i], Δn[i], K[i], Mp[i], Mm[i] = quantities(M, U[i], t, H, E0p, E0m)
        dat_utilde = column_stack((
            U/t,  # Utilde
            FC, Δn, K,  # frac condensado, fluctuaciones, kinetic
            Mp, Mm  # chem pot + y -
        ))
        return dat_utilde


# -----------------------------------------------------------------------------------
# PROGRAMA PRINCIPAL
# -----------------------------------------------------------------------------------
# guarda en archivos *.csv los datos de todas las gráficas
# corriendo hasta N=5 ~11s, N=6 ~70s, N=7 ~320s, N=8 ~1850s
path = ['Data', 'Plots']  # Folders de datos y plots
for i in range(len(path)):
    os.makedirs(path[i]) if not os.path.exists(path[i]) else None

samples = 100  # No. de puntos para las gráficas
ttil = open('Data/ttilde.csv', 'w+')
util = open('Data/utilde.csv', 'w+')
for N in [3, 4, 5, 6, 7, 8]:
    M = N
    # Calculo los H una sóla vez por N,M
    H = hamiltonian(M, M)
    E0p = hamiltonian(M + 1, M)  # Para μ+
    E0m = hamiltonian(M - 1, M)  # Para μ-
    # Guardar datos en *.csv
    dat_ttilde = graphs(M, samples, H, E0p, E0m, "ttilde")
    dat_utilde = graphs(M, samples, H, E0p, E0m, "utilde")
    with open('Data/ttilde.csv', 'a') as ttil:
        savetxt(ttil, dat_ttilde, fmt='% 0.5f', delimiter=',',
                header='M={} t/U, fc, Δn2, kin, μ+, μ-'.format(M))
        ttil.write('\n\n')  # índices de gnuplot
    with open('Data/utilde.csv', 'a') as util:
        savetxt(util, dat_utilde, fmt='% 0.5f', delimiter=',',
                header='M={} U/t, fc, Δn2, kin, μ+, μ-'.format(M))
        util.write('\n\n')  # índices de gnuplot
ttil.close()
util.close()
# Fig 1 Zhang, + ~45s en guardar ambos archivos
for M in [6, 10]:
    Ham = open('Data/H_M{}.csv'.format(M), 'w')
    H_ = hamiltonian(M, M)
    H = hamiltonian_tU(H_[0], H_[1])
    savetxt(Ham, column_stack((
                  # índices i y j del hamiltoniano
                  H.nonzero()[0], H.nonzero()[1],
                  # valor de la entrada [i,j]
                  H.toarray()[H.nonzero()]
                  )),
            fmt='%0.5g', delimiter=',')
    Ham.close()
