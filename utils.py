# -*- coding: utf-8 -*-

"""
This module contains functions to create and manipulate quantum states, including creation and annihilation operators, time evolution, and transition probabilities. It is designed for use in quantum mechanics simulations, particularly in the context of spin systems.
"""

# Importing necessary libraries

import numpy as np
import itertools
import scipy.linalg
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import scipy.constants as sc
import ctypes
import os
# Constants

# load utils_backup c function
script_dir = os.path.dirname(os.path.abspath(__file__))
dll_path = os.path.join(script_dir, 'utils_backup.dll')
lib_utils_c = ctypes.CDLL(dll_path)
#lib_utils_c = ctypes.CDLL('./utils.dll')

# Définir les types de structures C (basés sur votre code)
class State(ctypes.Structure):
    _fields_ = [
        ("size", ctypes.c_int),
        ("occupancy", ctypes.POINTER(ctypes.c_int))
    ]

class StateList(ctypes.Structure):
    _fields_ = [
        ("count", ctypes.c_int),
        ("states", ctypes.POINTER(State))
    ]

class Matrix(ctypes.Structure):
    _fields_ = [
        ("dim", ctypes.c_longlong),
        ("matrix", ctypes.POINTER(ctypes.POINTER(ctypes.c_double)))
    ]


spins = ['u', 'd']  # spins up and down
eV = 1.602176634e-19  # eV to Joules conversion factor

# déclaration utile pour l'interface python - c :

lib_utils_c.get_hubbard_states.argtypes = [
    ctypes.c_int,               # N
    ctypes.c_int                # dim
    ]
lib_utils_c.get_hubbard_states.restype = ctypes.POINTER(StateList)

# Configuration des types pour hubbard_hamiltonian_matrix
lib_utils_c.hubbard_hamiltonian_matrix.argtypes = [
    ctypes.c_int,                    # N
    ctypes.POINTER(Matrix),         # t_matrix
    ctypes.c_double,                 # U
    ctypes.c_int                    # dim
]
lib_utils_c.hubbard_hamiltonian_matrix.restype = ctypes.POINTER(Matrix)

# Fonctions utiles pour la conversion :
# Fonctions utilitaires pour la conversion
def numpy_to_c_matrix(np_matrix):
    """Convertit une matrice NumPy en structure C"""
    rows, cols = np_matrix.shape
    
    # Créer la structure TMatrix
    t_matrix = Matrix()
    t_matrix.size = rows
    
    # Allouer et remplir la matrice
    t_matrix.t_matrix = (ctypes.POINTER(ctypes.c_double) * rows)()
    for i in range(rows):
        t_matrix.t_matrix[i] = (ctypes.c_double * cols)()
        for j in range(cols):
            t_matrix.t_matrix[i][j] = np_matrix[i, j]
    
    return t_matrix

def c_matrix_to_numpy(hamiltonian_c):
    """Convertit une HamiltonianMatrix C en matrice NumPy"""
    dim = hamiltonian_c.contents.dim
    np_matrix = np.zeros((dim, dim))
    
    for i in range(dim):
        for j in range(dim):
            np_matrix[i, j] = hamiltonian_c.contents.matrix[i][j]
    
    return np_matrix

def c_states_to_numpy(statelist_c):
    """Convertit une StateList C en array NumPy"""
    count = statelist_c.contents.count
    states = []
    
    for i in range(count):
        state = statelist_c.contents.states[i]
        size = state.size
        occupancy = [state.occupancy[j] for j in range(size)]
        states.append(occupancy)
    
    return np.array(states)

def make_matrix(matrix_py):
    dim = len(matrix_py)
    RowArrayType = ctypes.c_double * dim  # type pour une ligne

    # Créer un tableau de lignes
    rows = []
    for row in matrix_py:
        row_array = RowArrayType(*row)
        rows.append(row_array)

    # Créer le tableau de pointeurs vers les lignes
    MatrixType = ctypes.POINTER(ctypes.c_double) * dim
    matrix_ptrs = MatrixType(*(ctypes.cast(row, ctypes.POINTER(ctypes.c_double)) for row in rows))

    # Créer l'objet Matrix final
    matrix_c = Matrix()
    matrix_c.matrix = matrix_ptrs
    matrix_c.dim = dim

    return matrix_c, rows  # Important : garder `rows` en vie


# Functions
def python_list_to_c_state(python_list) :
    """
    Convertit une liste Python en structure C State
    
    Args:
        python_list: Liste d'entiers (généralement 0 et 1) représentant l'état
        
    Returns:
        Pointeur vers la structure C State
        
    Example:
        >>> init_state = [0, 1, 1, 0, 1, 0, 1, 0]
        >>> c_state = python_list_to_c_state(init_state)
    """
    
    # Définition de la structure State en Python/ctypes
    class State(ctypes.Structure):
        _fields_ = [
            ("size", ctypes.c_int),
            ("occupancy", ctypes.POINTER(ctypes.c_longlong))
        ]
    
    # Vérifications d'entrée
    if not isinstance(python_list, (list, tuple, np.ndarray)):
        raise TypeError("L'entrée doit être une liste, tuple ou array numpy")
    
    if len(python_list) == 0:
        raise ValueError("La liste ne peut pas être vide")
    
    # Conversion en liste Python si c'est un array numpy
    if isinstance(python_list, np.ndarray):
        python_list = python_list.tolist()
    
    # Vérification que les valeurs sont des entiers
    if not all(isinstance(x, (int, np.integer)) for x in python_list):
        raise ValueError("Tous les éléments doivent être des entiers")
    
    # Création de la structure State
    state = State()
    state.size = len(python_list)
    
    # Création du tableau C pour occupancy
    occupancy_array = (ctypes.c_longlong * len(python_list))()
    for i, value in enumerate(python_list):
        occupancy_array[i] = int(value)
    
    # Attribution du pointeur
    state.occupancy = ctypes.cast(occupancy_array, ctypes.POINTER(ctypes.c_longlong))
    




def get_U(a: float, er: float = 11.7):

    """
    Returns the value of the coulomb interaction U in meV for a given lenght constant a.

    Parameters
    a : float, length
    er : float, relative permittivity (default is 11.7 for GaAs)

    Returns:
    U : float, coulomb interaction in meV
    """

    return sc.e/(2*np.pi*sc.epsilon_0*er*np.sqrt(2*np.pi)*a)/1e-3


def get_a(U: float, er: float = 11.7):

    """
    Returns the value of the length constant a in m for a given coulomb interaction U in eV.

    Parameters:
    U : float, coulomb interaction
    er : float, relative permittivity (default is 11.7 for GaAs)

    Returns:
    a : float, length constant in m
    """

    return sc.e/(2*np.pi*sc.epsilon_0*er*np.sqrt(2*np.pi)*U)


def create_bit_strings(N):

    """
    Create all N-bit binary strings

    Parameters:
    N: int, number of site

    Returns:
    ret_array: array of int, list of all arrangements 
    """

    ret_arr = list(itertools.product([0,1], repeat=N))
    ret_arr = np.array([list(arr) for arr in ret_arr])

    return ret_arr
    

def time_evol_state(H, T, u, hbar=1):

    """
    Returns an array of statevectors corresponding to the time evolution of u under H, according to |u(t)> = exp(iHt/h)|u(0)>.

    Parameters:
    H: array, hamiltonien of the system
    T: array, times of the system, often a np.linspace
    u: array, state considered
    hbar: float, reduced Planck constant (default = 1)

    Returns:
    array, array of statevectors over time
    """

    return np.array([(time_evol_operator(H, t, hbar) @ u) for t in tqdm(T, leave=False)])


def transition_probability_over_time(left_state, right_states):

    """
    This function returns an array corresponding to |<left_state|right_state>|^2 over time. It assumes right_states is an array of the T statevectors over time

    Parameters:
    left_state: array, reference state of the system
    right_states: array, array of statevectors over time

    Returns:
    array, array of probabilities of the left_state over time
    """
    
    ret_component = np.array([(np.vdot(left_state, right_state)) for right_state in right_states])
    ret_component = np.square(np.absolute(ret_component))

    return ret_component


def generate_base_states(N):

    """
    Generates the simple basis states of length N (e.g. [1,0], [0,1]) for a system of N sites

    Parameters:
    N: int, number of sites

    Returns:
    array, array of the simple basis states
    """

    return np.eye(N)


def time_evol_operator(H, t, hbar=1):

    """
    This function returns the time evolution operator U for a given Hamiltonian H and time t. U = exp(-iHt/hbar)

    Parameters:
    H: array, hamiltonian of the system
    t: float, time of the system
    hbar: float, reduced Planck constant (default = 1)

    Returns:
    array, time evolution operator
    """

    return scipy.linalg.expm(-1j * H * t / hbar)


def prob_over_time(H,T,u,v,transpose = True):

    """
    Returns the transition probability over time T for a given Hamiltonian H, initial state u, and observed state v.
    The function computes the time evolution of the state u under the Hamiltonian H and then calculates the transition probability to the state v.

    Parameters:
    H: array, Hamiltonian of the system
    T: array, time points at which to evaluate the transition probability
    u: array, initial state of the system
    v: array, observed state of the system
    transpose: bool, whether to transpose u (default is True)
    
    Returns:
    array, transition probability over time
    """

    if transpose: 
        u = np.atleast_2d(u).T
    # Ensure u is a column vector

    U = time_evol_state(H,T,u)
    ret_component = transition_probability_over_time(v,U)
    return ret_component


def get_label(u):

    """
    Returns the label of the state u in the form of a string with ↑ and ↓ symbols.
    The function assumes that u is a binary array where 1 represents an occupied state (↑ or ↓) and 0 represents an unoccupied state.

    Parameters:
    u: array, binary array representing the state

    Returns:
    str, label of the state in the form of a string with ↑ and ↓ symbols
    """

    # [:-1] to remove the last '|' character

    return ''.join(['0|' if not u[2*i] and not u[2*i+1] else 
                    '↑↓|' if u[2*i] and u[2*i+1]  else 
                    '↑ |' if u[2*i] and not u[2*i+1] else 
                    ' ↓|' for i in range(len(u) // 2)])[:-1]

def get_sampling_timestep(H_py):

    """
    Returns the sampling time steps for the Hamiltonian H.
    The function computes the eigenvalues of the Hamiltonian and returns the time steps based on the maximum eigenvalue.

    Parameters:
    H: array, Hamiltonian of the system

    Returns:
    float, sampling time step based on the energy difference of the eigenvalues
    """

    E = np.linalg.eigvals(H_py)
    Delta_E = np.abs(E.max()-E.min())

    return (np.pi*sc.hbar / Delta_E)


def get_hopping_simple_matrix(N, t):

    """
    Returns a simple hopping matrix for a 1D lattice with N sites.
    The matrix is tridiagonal with t on the off-diagonal elements.

    Parameters:
    N: int, number of sites
    t: float, hopping integral

    Returns:
    array, hopping matrix
    """

    t_matrix_py = np.zeros((N, N))
    for i in range(N-1):
        t_matrix_py[i, i+1] = t
        t_matrix_py[i+1, i] = t
    return t_matrix_py

def top_hubbard_states(T, U, t_matrix_py, init_binary_state=[0,1,1,0,1,0,1,0], top_n=4, figsize=(12,6), nbr_pts=1000, N = 2):

    """
    Plot the top_n Hubbard states with the highest transition probabilities over time.

    Parameters:
    T : float
    U : float
        On-site interaction strength.
    t_matrix_py : array

        Total time for the simulation.
        Time points at which to evaluate the transition probabilities.
        Hopping integral matrix (symmetric NxN matrix).
    init_binary_state : array
        Initial binary state of the system.
    top_n : int
        Number of top states to plot.
    figsize : tuple
        Size of the figure for plotting.
    nbr_pts : int
        Number of points in the time array.

    Returns:
    None
        Displays the plot of the top_n Hubbard states with the highest transition probabilities over time.
    """

    # U = U * eV  # Convert U from eV to Joules
    # t_matrix_py = t_matrix_py * eV  # Convert t from eV to Joules

    # # Number of sites
    # N = len(init_binary_state) // 2
    # dim = len(init_binary_state)

    #display = True

    # Hamiltonian and states
    #H = hubbard_hamiltonian_matrix(N, t_matrix, U)
    #H = lib_utils_c.hubbard_hamiltonian_matrix(N, t_matrix, U, statelist)

    # Avant l'appel à get_hubbard_states :
    #display = True
    #c_states = lib_utils_c.get_hubbard_states(N, 2*N) # it work
    #print("c_state = ", c_states)
    #print("test__2")
    #dim = c_states.contents.count
    #print("test__1")
    # #statelist = lib_utils_c.get_hubbard_states(N, dim) #
    # statelist_py = c_states_to_numpy(c_states)
    # #print(statelist_py)
    # print("test_2")
    # #Vérifier si c_states est valide
    # if not c_states:
    #     print("ERROR: get_hubbard_states returned NULL")
    #     return None, None, None

    # print(f"DEBUG: c_states pointer: {c_states}")
    
    # # Essayer de lire le contenu
    # try:
    #     dim = c_states.contents.count
    #     print(f"DEBUG: Retrieved dim={dim} from C function")
    # except Exception as e:
    #     print(f"ERROR: Cannot read c_states.contents.count: {e}")
    #     return None, None, None

    # # Continuer avec le reste du code seulement si tout va bien jusqu'ici
    # try:
    #     statelist_py = c_states_to_numpy(c_states)
    #     print(f"DEBUG: Successfully converted to numpy, shape: {statelist_py.shape}")
    # except Exception as e:
    #     print(f"ERROR: Cannot convert c_states to numpy: {e}")
    #     return None, None, None

   

    # var = ctypes.pointer(t_matrix_c)
    # print("Hello world !")

    # H_c = lib_utils_c.hubbard_hamiltonian_matrix( 
    #     ctypes.c_int(N), 
    #     var, # probleme ?
    #     ctypes.c_double(U), 
    #     ctypes.c_int(dim), 
    #     ctypes.c_int(V)
    #     )
    
    # print("Bientot")

    
    # Ta liste Python
    init_binary_state = [0, 1, 1, 0, 1, 0, 1, 0]
    size = len(init_binary_state)

    # Convertir en tableau C de int
    OccupancyArrayType = ctypes.c_int * size
    occupancy_array = OccupancyArrayType(*init_binary_state)


    #MatrixArrayType = ctypes.c_double * dim
    matrix_py = 0.0394e-3*get_hopping_simple_matrix(4, 1)
    scaled_matrix_py = [[0.0394e-3 * val for val in row] for row in matrix_py]
    t_matrix_c, _row_refs = make_matrix(scaled_matrix_py)

    N = 4
    dim = 2 * N
    t_matrix_c # probleme ?
    init_binary_state_c = State(
        size=size,
        occupancy=occupancy_array
    )
    top_n=4
    nbr_pts =1000

    lib_utils_c.top_hubbard_states_calculation(
        ctypes.c_double(T), # Any
        ctypes.c_double(U), # Any
        ctypes.pointer(t_matrix_c), # Matrix*
        ctypes.pointer(init_binary_state_c), # State*
        ctypes.c_int(top_n), # Literal[4] ?
        ctypes.c_int(nbr_pts) # Literal[1000] ?
    )
    
    H_c = lib_utils_c.hubbard_hamiltonian_matrix( 
        ctypes.c_int(N),  
        ctypes.pointer(t_matrix_c),  
        ctypes.c_double(U),  
        ctypes.c_int(dim)
    )
    H_py = c_matrix_to_numpy(H_c)
    # print(H_py)


    dt = get_sampling_timestep(H_py)
    print(dt)
    nbr_pts   = int(T/dt) # nombre théorique de points

    nbr_pts   = min(1.2*nbr_pts, 500000)         
    if nbr_pts == 500000:       
        #print(f"Trop long")
        print("Trop long")


    T = np.linspace(0, T, int(nbr_pts))
    T = T/sc.hbar

    # Initial state
    idx0 = np.where(np.all(c_states == init_binary_state, axis=1))[0]
    if idx0.size == 0:
        raise ValueError("Inital state not valid")
    psi0 = np.zeros(len(c_states), dtype=complex)
    psi0[idx0[0]] = 1.0

    # Temporal evolution
    psi_t = time_evol_state(H_py, T, psi0)
    print(psi_t)

    dim = len(statelist_py)
    probs = np.zeros((dim, len(T)), dtype=float)
    for i in tqdm(range(dim), desc="États", leave=False):
        v = np.zeros(dim, dtype=complex)
        v[i] = 1.0
        probs[i] = transition_probability_over_time(v, psi_t)

    # Max probabilities and top indices
    max_probs = probs.max(axis=1)
    top_idxs  = np.argsort(max_probs)[::-1][:top_n]

    display = True
    # Plotting only the top states 
    if display:
        plt.figure(figsize=figsize)
        for i in tqdm(top_idxs, desc="Tracé des top états"):
            label = f"|{get_label(c_states[i])}>"
            plt.plot(sc.hbar * T, probs[i], label=label)

        plt.legend(loc='best', title=f'Top {top_n} états')
        plt.xlabel('Temps')
        plt.ylim(0.3, 1)
        plt.ylabel('Probabilité')
        plt.title(f"Top {top_n} probabilités d'occupation (N={N}, U={U})")
        plt.tight_layout()
        plt.show()

    return sc.hbar*T, probs[top_idxs], c_states[top_idxs]

# Fonction de test pour get_hubbard_states
def test_get_hubbard_states(N, dim):
    """Test de la fonction get_hubbard_states"""
    try:
        print(f"Test get_hubbard_states avec N={N}")
        
        # Appeler la fonction C
        result = lib_utils_c.get_hubbard_states(N, dim)
        
        if not result:
            print("Erreur: get_hubbard_states a retourné NULL")
            return None
        
        # Convertir en NumPy
        states = c_states_to_numpy(result)
        print(f"Nombre d'états générés: {len(states)}")
        print(f"Premiers états:")
        for i, state in enumerate(states[:5]):  # Afficher les 5 premiers
            print(f"  État {i}: {state}")
        
        return states
        
    except Exception as e:
        print(f"Erreur lors du test: {e}")
        return None

# Configuration des types pour get_hubbard_states

# N = 2
# dim = 2*N
#test_get_hubbard_states(N, dim)
T = 4
U = 4
t_matrix_py = get_hopping_simple_matrix(N=3, t=0.42)
init_binary_state=[0,1,1,0,1,0,1,0]
top_n=4
figsize=(12,6)
nbr_pts=1000
N = 2

T, Tmp, _ = top_hubbard_states(T, U, t_matrix_py, init_binary_state, top_n, figsize, nbr_pts, N)

# lib_utils_c.simple.argtypes = []
# lib_utils_c.simple.restype = ctypes.c_int

# resultat = lib_utils_c.simple()
# print(resultat)

