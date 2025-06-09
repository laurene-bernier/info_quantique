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

class HamiltonianMatrix(ctypes.Structure):
    _fields_ = [
        ("dim", ctypes.c_int),
        ("matrix", ctypes.POINTER(ctypes.POINTER(ctypes.c_double)))
    ]

class TMatrix(ctypes.Structure):
    _fields_ = [
        ("size", ctypes.c_int),
        ("t_matrix", ctypes.POINTER(ctypes.POINTER(ctypes.c_double)))
    ]

spins = ['u', 'd']  # spins up and down
eV = 1.602176634e-19  # eV to Joules conversion factor

# déclaration utile pour l'interface python - c :

# Configuration des types pour get_hubbard_states
lib_utils_c.get_hubbard_states.argtypes = [
    ctypes.c_int,               # N
    ctypes.c_int                # dim
    ]
lib_utils_c.get_hubbard_states.restype = ctypes.POINTER(StateList)

# Configuration des types pour hubbard_hamiltonian_matrix
lib_utils_c.hubbard_hamiltonian_matrix.argtypes = [
    ctypes.c_int,                    # N
    ctypes.POINTER(TMatrix),         # t_matrix
    ctypes.c_double,                 # U
    ctypes.c_int,                    # dim
    ctypes.c_int                     # V
    #ctypes.POINTER(StateList)       # statelist
]
lib_utils_c.hubbard_hamiltonian_matrix.restype = ctypes.POINTER(HamiltonianMatrix)

lib_utils_c.top_hubbard_states_calculation.argtypes = [
    ctypes.c_double,                 # T
    ctypes.c_double,                 # U  
    ctypes.POINTER(TMatrix),         # t_matrix_py (needs to be converted)
    ctypes.POINTER(State),           # init_binary_state (converted State)
    ctypes.c_int,                    # top_n
    ctypes.c_double,                 # figsize[0] - you might need to adjust this
    ctypes.c_double,                 # figsize[1] - you might need to adjust this
    ctypes.c_int                     # nbr_pts
]
lib_utils_c.top_hubbard_states_calculation.restype = ctypes.c_void_p  # o

# Fonctions utiles pour la conversion :
# Fonctions utilitaires pour la conversion
def python_list_to_c_state(python_list):
    """
    Converts a Python list to a C State structure
    
    Parameters:
    python_list: list, Python list representing the state (e.g., [0,1,1,0,1,0,1,0])
    
    Returns:
    State: C State structure
    """
    # Create the State structure
    state = State()
    state.size = len(python_list)
    
    # Allocate memory for the occupancy array
    state.occupancy = (ctypes.c_int * len(python_list))()
    
    # Fill the occupancy array with values from the Python list
    for i, value in enumerate(python_list):
        state.occupancy[i] = int(value)
    
    return state

def numpy_array_to_c_state(numpy_array):
    """
    Converts a NumPy array to a C State structure
    
    Parameters:
    numpy_array: np.array, NumPy array representing the state
    
    Returns:
    State: C State structure
    """
    return python_list_to_c_state(numpy_array.tolist())

def numpy_to_c_matrix(np_matrix):
    """Convertit une matrice NumPy en structure C"""
    rows, cols = np_matrix.shape
    
    # Créer la structure TMatrix
    t_matrix = TMatrix()
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


# Functions




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

# --> en c
# def creation(state, i, spin):

#     """
#     Creation operator acting on spin of site i (0-indexed) of state

#     Parameters:
#     state: array, state of the system
#     i: int, site wanted
#     spin: char, "u" or "d"

#     Returns:
#     ret_state: array, new state after creation operator is applied
#     """

#     # Calculate the index of the spin in the state array
#     idx = 2 * i if spin == 'u' else 2*i + 1

#     # Calculate the Fermi sign factor
#     S_i = np.abs(np.sum(state[0:idx]))
#     sign_factor = (-1) ** S_i

#     if not state.any():
#         return np.zeros(len(state))
#     elif np.abs(state[idx]) == 1:
#         return np.zeros(len(state))
#     else:
#         ret_state = state.copy()
#         ret_state[idx] = 1
#         return sign_factor*ret_state
    
# --> en c
# def annihilation(state, i, spin):

#     """
#     Annihilation operator acting on spin of site i (0-indexed) of state

#     Parameters:
#     state: array, state of the system
#     i: int, site wanted
#     spin: char, "u" or "d"

#     Returns:
#     ret_state: array, new state after annihilation operator is applied
#     """

#     # Calculate the index of the spin in the state array
#     idx = 2 * i if spin == 'u' else 2*i + 1

#     # Calculate the Fermi sign factor
#     S_i = np.abs(np.sum(state[0:idx]))
#     sign_factor = (-1) ** S_i
#     if not state.any():
#         return np.zeros(len(state))
#     elif state[idx] == 0:
#         return np.zeros(len(state))
#     else:
#         ret_state = state.copy()
#         ret_state[idx] = 0
#         return sign_factor*ret_state
    
# --> en c
# def number_operator(state, i, spin):

#     """
#     Number operator acting on spin of site i (0-indexed) of state. Returns 0 or 1 whether the site i have the spin wanted or not

#     Parameters:
#     state: array, state of the system
#     i: int, site wanted
#     spin: char, "u" or "d"

#     Returns:
#     int, 0 or 1 depending on whether the site i has the spin wanted
#     """

#     idx = 2 * i if spin == 'u' else 2 * i + 1

#     return int(state[idx])
    

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


# def hopping_term_sign_factor(state, i, k, spin):

#     """
#     This function returns the sign factor for the hopping term in the Hamiltonian.

#     Parameters:
#     state: array, state of the system
#     i: int, site index of the initial state
#     k: int, site index of the final state
#     spin: char, "u" or "d"

#     Returns:
#     int, sign factor for the hopping term
#     """

#     # Hopping is equivalent to annihilation at i and creation at k
#     # The sign factor is (-1)^(S_i + S_k), where S_i and S_k are the number of spins at sites i and k respectively

#     idx_i = 2 * i if spin == 'u' else 2*i + 1
#     idx_k = 2 * k if spin == 'u' else 2*k + 1

#     min_idx = min(idx_i, idx_k)
#     max_idx = max(idx_i, idx_k)

#     S = np.sum(state[min_idx+1:max_idx])


#     return (-1) ** (S)

# --> en c
# def get_hubbard_states(N):

#     """
#     Generates all possible Hubbard states for a system of N sites.
#     Each state is represented as a binary array of length 2N, where the first N bits represent spin-up electrons and the last N bits represent spin-down electrons.

#     Parameters:
#     N: int, number of sites

#     Returns:
#     array, array of all possible Hubbard states
#     """

#     dim = 2 * N  # 2 états (↑ et ↓) par site
#     all_states = []

#     # On choisit N positions parmi 2N pour y mettre les électrons (1s)
#     for occ_indices in itertools.combinations(range(dim), N):
#         state = np.zeros(dim, dtype=int)
#         state[list(occ_indices)] = 1
#         all_states.append(state)

#     return np.array(all_states)

# --> en c
# def hubbard_hamiltonian_matrix(N, t, U, V = 0,states = None):

#     """
#     Returns the Hubbard Hamiltonian matrix for a system of N sites.
    
#     Parameters:
#     N: int, number of sites
#     t: array, hopping integral matrix (symmetric NxN matrix), t[i][j] represents the hopping amplitude between sites i and j
#     U: float, on-site interaction strength
#     V: float, nearest-neighbor interaction strength (default is 0)
#     states: array, optional, list of states to consider (if None, all possible Hubbard states are generated)

#     Returns:
#     array, Hubbard Hamiltonian matrix in the basis of all possible states
#     """
    
#     if states is None:
#         states = get_hubbard_states(N)  # Get all possible Hubbard states
#         dim = len(states)  # Dimension of the Hilbert space
#     else:
#         dim = len(states)
#     H = np.zeros((dim, dim))
    
#     # Loop over all states (rows)
#     for i in range(dim):
#         state_i = states[i]
        
#         # Loop over all states (columns)
#         for j in range(dim):
#             state_j = states[j]
            
#             # Diagonal elements: Coulomb interaction term 
#             if i == j:
#                 for site in range(N):
#                     # Check if both up and down spins are present at the site
#                     n_up = number_operator(state_i, site, 'u')
#                     n_down = number_operator(state_i, site, 'd')
#                     H[i, j] += U * n_up * n_down
                
#                 if V != 0:
#                     for site1 in range(N-1):
#                         site2 = site1 + 1
#                         n1 = number_operator(state_i, site1, 'u') + number_operator(state_i, site1, 'd')
#                         n2 = number_operator(state_i, site2, 'u') + number_operator(state_i, site2, 'd')
#                         H[i,i] += V * n1 * n2
                
#             # Off-diagonal: Hopping terms
#             else:
#                 # Determine if states i and j differ by a single hopping event
#                 for site1 in range(N):
#                     # Hubbard nearest-neighbor hopping
#                     for site2 in (site1-1, site1+1):
#                         if 0 <= site2 < N:
#                             for spin in ['u','d']:
#                                 temp = annihilation(state_i, site1, spin)
#                                 # Check if there is a spin to move at site1 with spin
#                                 if np.any(temp):
#                                     final = creation(temp, site2, spin) # 0 if already occupied
                                   
#                                     if np.array_equal(np.abs(final), state_j):
#                                         sign = hopping_term_sign_factor(state_i, site1, site2, spin)
#                                         H[i, j] -= t[site1,site2] * sign
    
#     return H


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


def top_hubbard_states(T, U, t_matrix_py, init_binary_state=[0,1,1,0,1,0,1,0], top_n=4, figsize=(12,6), nbr_pts=1000):

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
    # print(f"DEBUG: About to call get_hubbard_states with N={N}, dim={dim}")
    
    # print(f"DEBUG: Types: N={type(N)}, dim={type(dim)}")
    # #display = True
    # c_states = lib_utils_c.get_hubbard_states(N, 2*N)
    # #print("c_state = ", c_states)

    # dim = c_states.contents.count

    # statelist = lib_utils_c.get_hubbard_states(N, dim)
    
    #print(statelist_py)

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

    V = 0
    N = 4
    t_matrix_c = numpy_to_c_matrix(t_matrix_py)

    # H_c = lib_utils_c.hubbard_hamiltonian_matrix(
    #     ctypes.c_int(N), 
    #     ctypes.byref(t_matrix_c), 
    #     ctypes.c_double(U), 
    #     ctypes.c_int(dim), 
    #     ctypes.c_int(V)
    #     )
    
    init_binary_state_c = python_list_to_c_state(init_binary_state)

    lib_utils_c.top_hubbard_states_calculation(
        ctypes.c_double(T), 
        ctypes.c_double(U), 
        ctypes.byref(t_matrix_c), 
        ctypes.byref(init_binary_state_c), 
        ctypes.c_int(top_n=4), 
        ctypes.c_int(nbr_pts=1000)
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

    statelist_py = c_states_to_numpy(statelist)

    # Initial state
    idx0 = np.where(np.all(statelist == init_binary_state, axis=1))[0]
    if idx0.size == 0:
        raise ValueError("Inital state not valid")
    psi0 = np.zeros(len(statelist), dtype=complex)
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
            label = f"|{get_label(statelist[i])}>"
            plt.plot(sc.hbar * T, probs[i], label=label)

        plt.legend(loc='best', title=f'Top {top_n} états')
        plt.xlabel('Temps')
        plt.ylim(0.3, 1)
        plt.ylabel('Probabilité')
        plt.title(f"Top {top_n} probabilités d'occupation (N={N}, U={U})")
        plt.tight_layout()
        plt.show()

    return sc.hbar*T, probs[top_idxs], statelist[top_idxs]

# Fonction de test pour get_hubbard_states
def test_get_hubbard_states(N):
    """Test de la fonction get_hubbard_states"""
    try:
        print(f"Test get_hubbard_states avec N={N}")
        
        # Appeler la fonction C
        result = lib_utils_c.get_hubbard_states(
            ctypes.c_int(N)
            )
        
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