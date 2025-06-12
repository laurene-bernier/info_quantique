# -- coding: utf-8 --
"""
This module provides functions to analyze and visualize the different aspects of a system described by the Hubbard model.
It includes functions to convert CSV output from C to NPZ, plot the top Hubbard states, and visualize results.
"""

import os
import time as time_module
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.constants as sc
from ctypes import c_int, c_double, c_char_p, c_longlong, POINTER, Structure, CDLL, byref
import ctypes


# Constantes globales
eV = 1.602176634e-19
N = 2
U = 2.45e-3
T_final = 1e-10
init_binary_state = [0,1,1,0,1,0,1,0]
top_n = 4
figsize = (12,6)
nbr_pts = 10000
csv_file = "top_hubbard_states.csv"
npz_file = "top_hubbard_states.npz"

# Définition de la structure State en Python (miroir de la structure C)
class State(Structure):
    _fields_ = [
        ("size", c_longlong),
        ("occupancy", POINTER(c_longlong))
    ]

def load_c_library():
    """Charge la bibliothèque C avec gestion d'erreurs robuste."""
    library_paths = ['./libhubbard.dll', './libhubbard.so', 'libhubbard.dll', 'libhubbard.so']
    
    for lib_path in library_paths:
        try:
            if os.path.exists(lib_path):
                lib = CDLL(lib_path)
                # Configuration des types de la fonction C
                if hasattr(lib, 'top_hubbard_states_interface'):
                    # CORRECTION: Ajout de la structure State dans les arguments
                    lib.top_hubbard_states_interface.argtypes = [
                        c_int,          # N
                        c_double,       # U
                        c_double,       # T_final
                        c_int,          # nbr_pts
                        c_int,          # top_n
                        POINTER(State), # init_state - AJOUTÉ
                        c_char_p        # filename
                    ]
                    lib.top_hubbard_states_interface.restype = None
                    print(f"✅ Bibliothèque C chargée : {lib_path}")
                    return lib
                else:
                    print(f"⚠ Fonction 'top_hubbard_states_interface' non trouvée dans {lib_path}")
        except Exception as e:
            print(f"❌ Erreur de chargement {lib_path}: {e}")
    
    print("❌ Aucune bibliothèque C trouvée. Mode simulation activé.")
    return None

def python_list_to_state(python_list):
    """
    Convertit une liste Python en structure State* compatible avec C
    
    Parameters:
    python_list : list
        Liste Python d'entiers (0 ou 1) représentant l'état d'occupation
    
    Returns:
    tuple : (State, c_array) Structure State allouée et remplie avec référence au tableau C
    """
    # Créer une nouvelle instance de State
    state = State()
    
    # Définir la taille
    state.size = len(python_list)
    
    # Convertir la liste Python en tableau C
    ArrayType = c_longlong * len(python_list)
    c_array = ArrayType(*python_list)
    
    # Assigner le pointeur vers le tableau
    state.occupancy = c_array
    
    return state, c_array  # Retourner aussi c_array pour éviter la garbage collection

def convert_csv_to_npz(csv_file: str, npz_file: str):
    df = pd.read_csv(csv_file, converters={'occupation': lambda x: x})
    print(f"📄 Colonnes dans le fichier CSV : {df.columns.tolist()}")


    # Ignore columns that cannot be converted
    df = df.drop(columns=['occupation'], errors='ignore')

    times = sorted(df['t_idx'].unique())
    data = {}

    for idx in df['idx'].unique():
        df_idx = df[df['idx'] == idx]
        if len(df_idx) == len(times):
            data[f'state_{int(idx)}'] = df_idx.sort_values('t_idx')['proba'].values
        else:
            print(f"⚠ Incomplete data for state {idx}, skipping...")

    np.savez(npz_file, time=np.array(times), **data)
    return True

def load_data_from_npz(npz_file: str):
    """Charge les données NPZ en toute sécurité."""
    try:
        if not os.path.exists(npz_file):
            print(f"❌ Fichier NPZ non trouvé : {npz_file}")
            return None, None
        
        # CORRECTION: Essayer d'abord sans allow_pickle, puis avec si nécessaire
        try:
            npz_data = np.load(npz_file, allow_pickle=False)
        except ValueError as ve:
            if "pickled" in str(ve).lower():
                print("⚠ Fichier contient des données pickled, chargement avec allow_pickle=True")
                npz_data = np.load(npz_file, allow_pickle=True)
            else:
                raise ve
            
        time = npz_data['time']
        data = {k: npz_data[k] for k in npz_data.files if k != 'time'}
        print(f"✅ Données NPZ chargées : {len(data)} états")
        return time, data
    except Exception as e:
        print(f"❌ Erreur de chargement de {npz_file}: {e}")
        return None, None

def clean_old_files():
    """Nettoie les anciens fichiers pour éviter les conflits."""
    files_to_clean = [csv_file, npz_file]
    for file_path in files_to_clean:
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                print(f"🗑 Ancien fichier nettoyé: {file_path}")
            except Exception as e:
                print(f"⚠ Impossible de nettoyer {file_path}: {e}")

def generate_mock_data(nbr_pts, top_n):
    """Génère des données simulées si rien n'est disponible."""
    print("🔄 Génération de données simulées...")
    time = np.linspace(0, T_final, nbr_pts)
    data = {}
    for i in range(top_n):
        freq = 0.5 + i * 0.3
        phase = i * np.pi / 4
        amplitude = 0.8 - i * 0.1
        offset = 0.4 + i * 0.05
        prob = amplitude * np.exp(-time/20) * np.cos(freq * time + phase)**2 + offset
        prob = np.clip(prob, 0, 1)
        data[f'state_{i}'] = prob
    return time, data

def plot_top_hubbard_states(time, data, top_n=5, figsize=(12,6)):
    """Trace les top_n états avec les plus fortes probabilités."""
    if not data:
        print("❌ Aucune donnée à tracer")
        return None
        
    fig, ax = plt.subplots(figsize=figsize)
    colors = plt.cm.tab10(np.linspace(0, 1, min(len(data), top_n)))

    plotted_count = 0
    for i, (state_name, probabilities) in enumerate(data.items()):
        if plotted_count >= top_n:
            break
        
        time_plot = time * sc.hbar
        label = f"État {plotted_count+1}" if state_name.startswith('state_') else state_name
        ax.plot(time_plot, probabilities, label=label, color=colors[plotted_count], linewidth=2)
        plotted_count += 1

    ax.set_xlabel("Temps (J·s)")
    ax.set_ylabel("Probabilité d'occupation")
    ax.set_title(f"Top {plotted_count} probabilités d'occupation (N={N}, U={U} eV)")
    ax.set_ylim(0, 1)
    ax.legend(loc='best', title=f'Top {plotted_count} états')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    return fig

def call_c_library(lib, init_state):
    """Appelle la bibliothèque C de manière sécurisée."""
    if lib is None:
        print("❌ Bibliothèque C non disponible")
        return False
        
    try:
        print("📞 Appel de la bibliothèque C...")
        # Conversion de l'état initial
        state_struct, c_array = python_list_to_state(init_state)
        
        # CORRECTION: Passer la structure par référence
        lib.top_hubbard_states_interface(
            N, 
            c_double(U), 
            c_double(T_final), 
            nbr_pts, 
            top_n, 
            byref(state_struct),  # CORRECTION: Utiliser byref pour passer par référence
            csv_file.encode('utf-8')
        )
        
        # Attendre que le fichier soit écrit
        print("⏳ Attente de l'écriture du fichier...")
        max_wait = 10  # Attendre maximum 10 secondes
        wait_time = 0
        while not os.path.exists(csv_file) and wait_time < max_wait:
            time_module.sleep(0.5)
            wait_time += 0.5
        
        if os.path.exists(csv_file):
            print("✅ Appel C terminé, fichier CSV généré")
            return True
        else:
            print("❌ Appel C terminé mais aucun fichier CSV généré")
            return False
        
    except Exception as e:
        print(f"❌ Erreur appel C: {e}")
        import traceback
        traceback.print_exc()
        return False

def top_hubbard_states_graph(init_state):
    """Routine principale d'analyse et de visualisation."""
    print("🧠 Début de l'analyse des états de Hubbard...")
    print(f"   Paramètres: N={N}, U={U} eV, T_final={T_final}, nbr_pts={nbr_pts}, top_n={top_n}")
    
    # CORRECTION: Nettoyer les anciens fichiers
    clean_old_files()
    
    # Chargement de la bibliothèque C
    lib = load_c_library()
    
    # Tentative d'appel à la bibliothèque C
    c_success = call_c_library(lib, init_state)
    
    # Tentative de conversion CSV → NPZ
    if c_success and os.path.exists(csv_file):
        if convert_csv_to_npz(csv_file, npz_file):
            # Nettoyage optionnel du CSV
            try:
                os.remove(csv_file)
                print("🗑 Fichier CSV nettoyé")
            except:
                pass
    
    # Chargement des données
    time, data = load_data_from_npz(npz_file)
    
    # Si pas de données, utiliser des données simulées
    if time is None or data is None or len(data) == 0:
        print("⚠ Utilisation de données simulées.")
        time, data = generate_mock_data(nbr_pts, top_n)
    
    # Tracé
    print("📊 Génération du graphique...")
    fig = plot_top_hubbard_states(time, data, top_n, figsize)
    
    if fig is not None:
        print("✅ Analyse terminée avec succès.")
        print(f"   Données tracées pour {len(data)} états")
    else:
        print("❌ Échec de la génération du graphique")
    
    return time, data, fig

# Point d'entrée principal
def main():
    """Fonction principale avec gestion d'erreurs complète."""
    try:
        print("=" * 60)
        print("🚀 ANALYSE DES ÉTATS DE HUBBARD")
        print("=" * 60)
        
        # Conversion de l'état initial
        init_state_struct, c_array = python_list_to_state(init_binary_state)
        print(f"État initial: {init_binary_state}")
        
        # Lancement de l'analyse
        time, data, fig = top_hubbard_states_graph(init_binary_state)
        
        if time is not None and data:
            print(f"✔ Analyse terminée. Données pour {len(data)} états.")
            print(f"  Durée simulée: {T_final} unités de temps")
            print(f"  Points de données: {len(time)}")
        else:
            print("❌ Échec de l'analyse.")
            
    except KeyboardInterrupt:
        print("\n🛑 Analyse interrompue par l'utilisateur")
    except Exception as e:
        print(f"❌ Erreur critique: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()