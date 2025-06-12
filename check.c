// Inclure le header avec les macros de vérification
#include "malloc_check.h"  // Le fichier précédent
#include "utils_backup.h"

// =============================================================================
// VOS FONCTIONS ADAPTÉES AVEC VÉRIFICATIONS
// =============================================================================

Matrix* create_tridiagonal_matrix_safe(int N) {
    Matrix* matrix;
    
    // Allocation sécurisée de la structure principale
    SAFE_MALLOC(matrix, sizeof(Matrix), Matrix);
    
    matrix->dim = N;
    
    // Allocation sécurisée du tableau de pointeurs
    SAFE_MALLOC(matrix->matrix, N * sizeof(double*), double*);
    
    // Allocation sécurisée de chaque ligne
    for (int i = 0; i < N; i++) {
        SAFE_MALLOC(matrix->matrix[i], N * sizeof(double), double);
        // Initialiser à zéro
        memset(matrix->matrix[i], 0, N * sizeof(double));
    }
    
    // Remplissage de la matrice tridiagonale
    for (int i = 0; i < N; i++) {
        if (i > 0)
            matrix->matrix[i][i - 1] = 1.0;
        if (i < N - 1)
            matrix->matrix[i][i + 1] = 1.0;
    }
    
    return matrix;
}

Matrix* allocate_memory_matrix_safe(int dim) {
    Matrix *H;
    
    SAFE_MALLOC(H, sizeof(Matrix), Matrix);
    H->dim = dim;
    
    SAFE_MALLOC(H->matrix, dim * sizeof(double*), double*);
    
    for(int i = 0; i < dim; i++){
        SAFE_CALLOC(H->matrix[i], dim, sizeof(double), double);
    }
    
    return H;
}

State* annihilation_safe(State* state, int i, char spin) {
    if (!state || !state->occupancy) {
        fprintf(stderr, "ERREUR: État d'entrée invalide dans annihilation\n");
        return NULL;
    }
    
    int idx = (spin == 'u') ? 2 * i : 2 * i + 1;
    
    // Vérifier les limites
    if (idx >= state->size) {
        fprintf(stderr, "ERREUR: Index %d hors limites (taille: %lld)\n", idx, state->size);
        return NULL;
    }
    
    // Calcul du facteur de signe de Fermi
    long long summed_state = 0;
    for(long long j = 0; j < idx; j++) {
        summed_state += state->occupancy[j];
    }
    int sign_factor = (summed_state % 2 == 0) ? 1 : -1;
    
    State *new_state;
    SAFE_MALLOC(new_state, sizeof(State), State);
    
    new_state->size = state->size;
    SAFE_MALLOC(new_state->occupancy, new_state->size * sizeof(long long), long long);
    
    // Vérifier si l'annihilation est possible
    if (state->occupancy[idx] == 0) {
        // Retourner un état nul
        memset(new_state->occupancy, 0, new_state->size * sizeof(long long));
        return new_state;
    }
    
    // Copier l'état et appliquer l'annihilation
    memcpy(new_state->occupancy, state->occupancy, state->size * sizeof(long long));
    new_state->occupancy[idx] = 0;
    
    // Appliquer le signe de Fermi
    for (int j = 0; j < state->size; j++) {
        new_state->occupancy[j] *= sign_factor;
    }
    
    return new_state;
}

State* creation_safe(State* state, int i, char spin) {
    if (!state || !state->occupancy) {
        fprintf(stderr, "ERREUR: État d'entrée invalide dans creation\n");
        return NULL;
    }
    
    int idx = (spin == 'u') ? 2 * i : 2 * i + 1;
    
    // Vérifier les limites
    if (idx >= state->size) {
        fprintf(stderr, "ERREUR: Index %d hors limites (taille: %lld)\n", idx, state->size);
        return NULL;
    }
    
    // Calcul du facteur de signe
    long long summed_state = 0;
    for(long long j = 0; j < idx; j++) {
        summed_state += state->occupancy[j];
    }
    int sign_factor = (summed_state % 2 == 0) ? 1 : -1;
    
    State *new_state;
    SAFE_MALLOC(new_state, sizeof(State), State);
    
    new_state->size = state->size;
    SAFE_MALLOC(new_state->occupancy, new_state->size * sizeof(long long), long long);
    
    // Vérifier le principe d'exclusion de Pauli
    if (state->occupancy[idx] == 1) {
        // Retourner un état nul
        memset(new_state->occupancy, 0, new_state->size * sizeof(long long));
        return new_state;
    }
    
    // Copier l'état et appliquer la création
    memcpy(new_state->occupancy, state->occupancy, state->size * sizeof(long long));
    new_state->occupancy[idx] = 1;
    
    // Appliquer le signe de Fermi
    for (int j = 0; j < state->size; j++) {
        new_state->occupancy[j] *= sign_factor;
    }
    
    return new_state;
}

StateList* get_hubbard_states_safe(int N) {
    int dim = 2 * N;
    CombinationList *combs = combinations_iterative(N, dim);
    if (!combs) {
        fprintf(stderr, "ERREUR: Impossible de générer les combinaisons\n");
        return NULL;
    }
    
    StateList *state_list;
    SAFE_MALLOC(state_list, sizeof(StateList), StateList);
    
    state_list->count = combs->count;
    SAFE_MALLOC(state_list->states, combs->count * sizeof(State), State);
    
    for (long long i = 0; i < combs->count; i++) {
        state_list->states[i].size = dim;
        SAFE_CALLOC(state_list->states[i].occupancy, dim, sizeof(long long), long long);
        
        // Mettre 1 aux positions occupées
        for (int j = 0; j < combs->combinations[i].size; j++) {
            long long index = combs->combinations[i].indices[j];
            if (index < dim) {  // Vérification de sécurité
                state_list->states[i].occupancy[index] = 1;
            }
        }
    }
    
    free_combination_list(combs);
    return state_list;
}

// =============================================================================
// VERSION AVEC NETTOYAGE AUTOMATIQUE
// =============================================================================

Matrix* hubbard_hamiltonian_matrix_safe(int N, Matrix* t_matrix, double U, int dim) {
    // Initialiser le système de nettoyage
    CleanupList* cleanup = init_cleanup_list();
    if (!cleanup) {
        fprintf(stderr, "ERREUR: Impossible d'initialiser le nettoyage automatique\n");
        return NULL;
    }
    
    StateList* statelist = get_hubbard_states_safe(N);
    if (!statelist) {
        printf("Erreur: impossible de générer les états\n");
        cleanup_all(cleanup);
        return NULL;
    }
    
    dim = statelist->count;
    int hilbert_dim = (int)(statelist->count);
    
    Matrix* H = allocate_memory_matrix_safe(hilbert_dim);
    if (!H) {
        printf("Erreur: impossible d'allouer la matrice Hamiltonienne\n");
        free_state_list(statelist);
        cleanup_all(cleanup);
        return NULL;
    }
    
    // Le reste du calcul...
    for(int i = 0; i < hilbert_dim; i++){
        State* state_i = &statelist->states[i];
        
        for(int j = 0; j < hilbert_dim; j++){
            State* state_j = &statelist->states[j];
            
            if(i == j){
                // Terme diagonal: interaction Coulombienne
                for (int site = 0; site < N; site++){
                    int n_up = number_operator(state_i, site, 'u');
                    int n_down = number_operator(state_i, site, 'd');
                    H->matrix[i][j] += U * n_up * n_down;
                }
            } else {
                // Termes de saut
                for(int site1 = 0; site1 < N; site1++){
                    int site2_list[2] = {site1-1, site1+1};
                    
                    for(int s = 0; s < 2; s++){
                        if(site2_list[s] >= 0 && site2_list[s] < N){
                            char spins[2] = {'u', 'd'};
                            
                            for (int l = 0; l < 2; l++) {
                                char spin = spins[l];
                                
                                State* temp = annihilation_safe(state_i, site1, spin);
                                if(temp && any(temp)){
                                    int site2 = site2_list[s];
                                    State* final = creation_safe(temp, site2, spin);
                                    
                                    if(final && state_equal(abs_state(final), state_j)){
                                        int sign = hopping_term_sign_factor(state_i, site1, site2, spin);
                                        H->matrix[i][j] -= t_matrix->matrix[site1][site2] * sign;
                                    }
                                    
                                    // Libérer les états temporaires
                                    if(final) {
                                        free(final->occupancy);
                                        free(final);
                                    }
                                }
                                
                                if(temp) {
                                    free(temp->occupancy);
                                    free(temp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    free_state_list(statelist);
    cleanup_all(cleanup);
    return H;
}

// =============================================================================
// MAIN AVEC GESTION D'ERREURS COMPLÈTE
// =============================================================================

int main_safe() {
    int N = 3;
    int temps = 4;
    double U = 2;
    int dim = 2 * N;
    
    printf("Début du programme avec N=%d\n", N);
    
    // Générer les états de Hubbard
    StateList* statelist = get_hubbard_states_safe(N);
    if (!statelist) {
        fprintf(stderr, "ERREUR: Impossible de générer les états de Hubbard\n");
        return EXIT_FAILURE;
    }
    printf("États générés avec succès: %lld états\n", statelist->count);
    
    // Créer la matrice de saut
    Matrix* t_matrix = create_tridiagonal_matrix_safe(N);
    if (!t_matrix) {
        fprintf(stderr, "ERREUR: Impossible de créer la matrice de saut\n");
        free_state_list(statelist);
        return EXIT_FAILURE;
    }
    printf("Matrice de saut créée avec succès\n");
    
    // Créer l'état initial
    State* init_binary_state;
    ALLOC_STATE(init_binary_state, dim);
    
    // Initialiser avec un motif par défaut (alternance spin up/down)
    for (int i = 0; i < dim; i += 2) {
        init_binary_state->occupancy[i] = 1; // spin up
    }
    printf("État initial créé\n");
    
    // Calculer l'Hamiltonien
    Matrix* H = hubbard_hamiltonian_matrix_safe(N, t_matrix, U, dim);
    if (!H) {
        fprintf(stderr, "ERREUR: Impossible de calculer l'Hamiltonien\n");
        free(init_binary_state->occupancy);
        free(init_binary_state);
        free_memory_matrix(t_matrix, N);
        free_state_list(statelist);
        return EXIT_FAILURE;
    }
    
    printf("Hamiltonien calculé avec succès\n");
    print_matrix(H);
    
    // Nettoyage final
    free(init_binary_state->occupancy);
    free(init_binary_state);
    free_memory_matrix(t_matrix, N);
    free_memory_matrix(H, H->dim);
    free_state_list(statelist);
    
    printf("Programme terminé avec succès\n");
    return EXIT_SUCCESS;
}