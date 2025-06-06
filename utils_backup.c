#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils_backup.h"

/**
    * @brief Returns the Hubbard Hamiltonian matrix for a system of N sites.
    
    Parameters:
    * @param N : int, number of sites
    * @param t : array, hopping integral matrix (symmetric NxN matrix), t[i][j] represents the hopping amplitude between sites i and j
    * @param U : float, on-site interaction strength
    * @param V : float, nearest-neighbor interaction strength (default is 0)
    * @param states : array, optional, list of states to consider (if None, all possible Hubbard states are generated)

    * @return array, Hubbard Hamiltonian matrix in the basis of all possible states
*/

// temporaire :
//int N; 
//int dim;
//double U;
//int V = 0;
double eV = 1.602176634e-19;
//char spin;

// Note à moi-mème : qq recherches sur des manières d'optimiser son code en c/conseils trouvés sur internet 
// -
// -


// Definition of basic useful function :
// 1) for memory gestion :
HamiltonianMatrix* allocate_memory_hamiltonian(int dim) { //function to allocate memory for a dynamic array of the Hamiltonian matrix :
    HamiltonianMatrix *H = malloc(sizeof(HamiltonianMatrix)); // = new HamiltonianMatrix;
    H->dim = dim;
    H->matrix = malloc(dim * sizeof(double*));
    for(int i = 0; i <dim; i++){
        H->matrix[i] = calloc(dim, sizeof(double));
    }
    return H;
}

void free_memory_hamiltonian(HamiltonianMatrix* H, int dim) { //function to free memory of the Hamiltonian matrix :
    if(H){ // test if H is a NULL pointer --> false, true if not
        for(int i= 0; i < dim; i++){
            free(H->matrix[i]);
        }
        free(H->matrix);
        free(H);
    } 
}

void free_combination_list(CombinationList *list) { // to free combination list
    if (list) {
        for (int i = 0; i < list->max_count; i++) {
            free(list->combinations[i].indices);
        }
        free(list->combinations);
        free(list);
    }
}

void free_state_list(StateList *list) { // To free state list
    if (list) {
        for (int i = 0; i < list->count; i++) {
            free(list->states[i].occupancy);
        }
        free(list->states);
        free(list);
    }
}


// 2) basic :
// definition of basic function
HamiltonianMatrix* initialize_matrix_with_zeros(HamiltonianMatrix* H, int dim){ // does what it is said in the title
    for(int i = 0; i<dim; i++){
        for(int j = 0; j<dim; j++){
            H->matrix[i][j] = 0;
        }
    }
    return H;
}

bool state_equal(State* state1, State* state2){ // test if state are equal in size and in content
    if(!state1 || !state2) return false;
    if(state1->size == state2->size) {
        for(int i = 0; i < state1->size; i++) {
            if (state1->occupancy[i] != state2->occupancy[i]) {
                return false;
            }
        }
        return true;
    }
}


State* abs_state(State* state_not_in_absolute_value){ // transform each value of state in an absolute value
    State* state_in_absolute_value = malloc(sizeof(State));
    state_in_absolute_value->size = state_not_in_absolute_value->size;
    state_in_absolute_value->occupancy = malloc(state_not_in_absolute_value->size * sizeof(int)); // need to allocate memory to occupancy
    for(int i = 0; i < state_not_in_absolute_value->size; i++){
        state_in_absolute_value->occupancy[i] = abs(state_not_in_absolute_value->occupancy[i]);
    }
    return state_in_absolute_value;
}

bool any(State* state){
    if (!state || !state->occupancy) return false;
    for (int i = 0; i < state->size; i++) {
        if (state->occupancy[i] != 0) return true;
    }
    return false;
}

State* make_a_vector_of_zero_state_lengthed(int dim, State* state){
    for(int i = 0; i<dim; i++){
            state->occupancy[i] = 0;
        }
    return state;
}

int min(int a, int b) { return (a < b) ? a : b; }
int max(int a, int b) { return (a > b) ? a : b; }

// 3) for the denombrement in get_hubard_state :
int binomial_coefficient(int n, int k) { // Calculation of the binomial coefficient C(n,k)
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    
    // Utiliser la propriété C(n,k) = C(n,n-k) pour optimiser
    if (k > n - k) k = n - k;
    
    long long result = 1;
    for (int i = 0; i < k; i++) {
        result = result * (n - i) / (i + 1);
    }
    return (int)result;
}

CombinationList* init_combination_list(int n, int k) { // Initialization of the combinations list
    CombinationList *list = malloc(sizeof(CombinationList));
    list->count = 0;
    list->max_count = binomial_coefficient(n, k);
    list->combinations = malloc(list->max_count * sizeof(Combination));
    
    // Initialiser chaque combinaison
    for (int i = 0; i < list->max_count; i++) {
        list->combinations[i].indices = malloc(k * sizeof(int));
        list->combinations[i].size = k;
    }
    
    return list;
}

CombinationList* combinations_iterative(int n, int k) { // Itérative version  (plus efficace pour certains cas)
    if (k > n || k < 0) return NULL;
    //printf("k & n : %d &  %d", k, n);
    CombinationList *list = init_combination_list(n, k);
    int *combination = malloc(k * sizeof(int));
    //printf(list);
    // Initialiser la première combinaison [0, 1, 2, ..., k-1]
    for (int i = 0; i < k; i++) {
        combination[i] = i;
        //printf(combination[i]);
    }
    
    do {
        // Copier la combinaison actuelle dans la liste
        for (int i = 0; i < k; i++) {
            list->combinations[list->count].indices[i] = combination[i];
        }
        list->count++;
        
        // Générer la combinaison suivante
        int i = k - 1;
        while (i >= 0 && combination[i] == n - k + i) {
            i--;
        }
        
        if (i < 0) break; // Plus de combinaisons
        
        combination[i]++;
        for (int j = i + 1; j < k; j++) {
            combination[j] = combination[j-1] + 1;
        }
        
    } while (true);
    
    free(combination);
    return list;
}

void print_combinations(CombinationList *list) { // To print combinations
    printf("Nombre de combinaisons: %d\n", list->count);
    for (int i = 0; i < list->count; i++) {
        printf("Combinaison %d: [", i);
        for (int j = 0; j < list->combinations[i].size; j++) {
            printf("%d", list->combinations[i].indices[j]);
            if (j < list->combinations[i].size - 1) printf(", ");
        }
        printf("]\n");
    }
}


// Other specific function :
int number_operator(State* state, int i, char spin) { // number_operator
    return (spin == 'u') ? (int)(state->occupancy[2 * i]) : (int)(state->occupancy[2 * i + 1]); //(int)
}

State* annihilation(State* state, int i, char spin, int dim){

    /**
    * @Annihilation operator acting on spin of site i (0-indexed) of state

    Parameters:
    state: array, state of the system
    i: int, site wanted
    spin: char, "u" or "d"

    Returns:
    * @return ret_state: array, new state after annihilation operator is applied
    */
    //Calculate the index of the spin in the state array
    int idx;
    if(spin == 'u') idx = 2 * i;
    else idx = 2*i + 1;

    //Calculate the Fermi sign factor
    int summed_state = 0;
    for(int i = 0; i < idx; i ++) summed_state += state->occupancy[i];
    int S_i = abs(summed_state);
    int sign_factor = pow((-1), S_i);

    if(any(state) == false){
        return make_a_vector_of_zero_state_lengthed(dim, state);
    }else if(abs(state->occupancy[idx]) == 0){
        return make_a_vector_of_zero_state_lengthed(dim, state);
    }else{
        State* ret_state = malloc(sizeof(state));
        ret_state->occupancy[idx]= 0;
        State* signed_ret_state = malloc(sizeof(State));
        signed_ret_state->size = dim;
        signed_ret_state->occupancy = malloc(dim * sizeof(int));
        for(int i = 0; i<dim; i++){ //dim = sizeof(state) normalement
            signed_ret_state->occupancy[i] = sign_factor * ret_state->occupancy[i];
        }
        return signed_ret_state;
    }
}
// attention --> copy() ?
State* creation(State* state, int i, char spin, int dim){
    
    /** @brief Creation operator acting on spin of site i (0-indexed) of state

    Parameters:
    * @param state: array, state of the system
    * @param i: int, site wanted
    * @param spin: char, "u" or "d"

    Returns:
    * @return ret_state: array, new state after creation operator is applied
    */

    // Calculate the index of the spin in the state array
    int idx;
    if(spin == 'u') idx = 2 * i;
    else idx = 2*i + 1;

    // Calculate the Fermi sign factor
    int summed_state = 0;
    for(int i = 0; i < idx; i ++) summed_state += state->occupancy[i];
    int S_i = abs(summed_state);
    int sign_factor = pow((-1), S_i);

    if(any(state) == false){
        return make_a_vector_of_zero_state_lengthed(dim, state);
    }else if(abs(state->occupancy[idx]) == 1){
        return make_a_vector_of_zero_state_lengthed(dim, state);
    }else{
        State* ret_state = state;
        ret_state->occupancy[idx] = 1;
        State* signed_ret_state;
        for(int i = 0; i<dim; i++){ //dim = sizeof(state) normalement
            signed_ret_state->occupancy[i] = sign_factor * ret_state->occupancy[i];
        }
        return signed_ret_state;
    }
}

int hopping_term_sign_factor(State* state_i, int i, int k, char spin){
    int idx_i = 0;
    int idx_k = 0;
    if(spin == 'u') idx_i = (2 * i);
    else idx_i = 2*i + 1;
    if(spin == 'u') idx_k = (2 * k); 
    else idx_k = 2*k + 1;

    int min_idx = min(idx_i, idx_k);
    int max_idx = max(idx_i, idx_k);

    int S = 0;
    for(int i = min_idx+1; i < max_idx; i++) S += state_i->occupancy[i];

    return (-1)^(S);
}

// Main function :
StateList* get_hubbard_states(int N, int dim) { // get_hubbard_states(N, dim)
    CombinationList *combs = combinations_iterative(dim, N); // dim = nombre de particule & N = nombre de site
    //if (!combs) return NULL;
    
    StateList *state_list = malloc(sizeof(StateList));
    state_list->count = combs->count;
    state_list->states = malloc(combs->count * sizeof(State));
    
    for (int i = 0; i < combs->count; i++) {
        // Initialiser l'état avec des zéros
        state_list->states[i].size = dim;
        state_list->states[i].occupancy = calloc(dim, sizeof(int));
        
        // Mettre 1 aux positions occupées
        for (int j = 0; j < combs->combinations[i].size; j++) {
            int index = combs->combinations[i].indices[j];
            state_list->states[i].occupancy[index] = 1;
        }
    }

    if (!state_list || !state_list->states) {
        free(state_list); // au cas où state_list est non NULL
        free_combination_list(combs);
        //return NULL;
    }
    
    free_combination_list(combs);
    return state_list;
}

HamiltonianMatrix* hubbard_hamiltonian_matrix(int N, t_matrix* t_matrix, double U, int dim, int V){
    StateList* statelist = get_hubbard_states(N, dim); // Get all possible Hubbard states
    dim = statelist->count;  // Dimension of the Hilbert space
    

    HamiltonianMatrix* H = allocate_memory_hamiltonian(dim);
    H = initialize_matrix_with_zeros(H, dim);
    
    
    // Loop over all states (rows)
    for(int i = 0; i < dim; i++){
        State* state_i = &statelist->states[i];
        
        // Loop over all states (columns)
        for(int j = 0; j < dim; j++){
            State* state_j = &statelist->states[j];
            
            // Diagonal elements: Coulomb interaction term 
            if(i == j){
                for (int site =0; site < N; site ++){
                    // Check if both up and down spins are present at the site
                    int n_up = number_operator(state_i, site, 'u');
                    int n_down = number_operator(state_i, site, 'd');
                    H->matrix[i][j] += U * n_up * n_down;
                }
                
                if(V != 0){
                    for (int site1 = 0; site1 < N; site1++){
                        int site2 = site1 + 1;
                        int n1 = number_operator(state_i, site1, 'u') + number_operator(state_i, site1, 'd');
                        int n2 = number_operator(state_i, site2, 'u') + number_operator(state_i, site2, 'd');
                        H->matrix[i][i] += V * n1 * n2;
                    }
                }
            }
                
            // Off-diagonal: Hopping terms
            else{
                // Determine if statelist i and j differ by a single hopping event
                for(int site1 = 0; site1 < N; site1++){
                    // Hubbard nearest-neighbor hopping
                    int site2_list[2] = {site1-1, site1+1};
                    for(int s = 0; s < 2; s++){
                        if(0 <= site2_list[s] && site2_list[s] < N){
                            char spins[2] = {'u', 'd'};
                            for (int l = 0; l < 2; l++) {
                                char spin = spins[l];
                                State* temp = annihilation(state_i, site1, spin, dim);
                                // Check if there is a spin to move at site1 with spin
                                if(temp != NULL){
                                    int site2 = site2_list[s];
                                    State* final = creation(temp, site2, spin, dim); // 0 if already occupied

                                    if(state_equal(abs_state(final), state_j)){ // problem
                                        int sign = hopping_term_sign_factor(state_i, site1, site2, spin); // --> a voir antisymétries des fermions
                                        H->matrix[i][j] -= t_matrix->t_matrix[site1][site2] * sign;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return H;
}

// --> problème de mémoire : crash 3 