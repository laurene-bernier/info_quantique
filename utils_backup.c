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

bool state_equal(State* state1, State* state2) {
    if (!state1 || !state2 || !state1->occupancy || !state2->occupancy) 
        return false;
    
    if (state1->size != state2->size) 
        return false;
    
    for(int i = 0; i < state1->size; i++) {
        if (state1->occupancy[i] != state2->occupancy[i]) {
            return false;
        }
    }
    return true;
}


State* abs_state(State* state_not_in_absolute_value){ // transform each value of state in an absolute value
    State* state_in_absolute_value = malloc(sizeof(State));
    state_in_absolute_value->size = state_not_in_absolute_value->size;
    state_in_absolute_value->occupancy = malloc(state_in_absolute_value->size * sizeof(int)); // need to allocate memory to occupancy // state_not_in_absolute_value->size * ?
    for(int i = 0; i < state_not_in_absolute_value->size; i++){
        state_in_absolute_value->occupancy[i] = abs(state_not_in_absolute_value->occupancy[i]);
    }
    return state_in_absolute_value;
}

bool any(State* state){
    //state_in_absolute_value->occupancy = malloc(state_not_in_absolute_value->size * sizeof(int));
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
long long binomial_coefficient(int n, int k) { // Calculation of the binomial coefficient C(n,k)
    // Version corrigée avec gestion d'erreurs
    // Vérifications d'entrée
    if (n < 0 || k < 0) return -1;  // Erreur : valeur négative
    if (k > n) return 0;            // Cas valide : C(n,k) = 0 si k > n
    if (k == 0 || k == n) return 1; // Cas de base
    
    // Optimisation : utiliser la propriété C(n,k) = C(n,n-k)
    if (k > n - k) k = n - k;
    
    long long result = 1;
    
    // Calcul sécurisé pour éviter les débordements
    for (int i = 0; i < k; i++) {
        // Vérifier le débordement potentiel
        if (result > LLONG_MAX / (n - i)) {
            return -1; // Erreur : débordement détecté
        }
        
        result = result * (n - i) / (i + 1);
    }
    
    return result;
}

CombinationList* init_combination_list(int k, int n) { // Initialization of the combinations list
    CombinationList *list = malloc(sizeof(CombinationList));
    list->count = 0;
    list->max_count = binomial_coefficient(n, k);
    long long abs_count = abs(list->max_count * sizeof(Combination));
    list->combinations = malloc(abs_count);
    
    // Initialiser chaque combinaison
    for (int i = 0; i < list->max_count; i++) {
        list->combinations[i].indices = malloc(list->max_count * sizeof(int));
        list->combinations[i].size = k;
    }
    
    return list;
}

CombinationList* combinations_iterative(int k, int n) { // Itérative version  (plus efficace pour certains cas)
    printf("Aie");
    if (k > n || k < 0) return NULL;
    //printf("k & n : %d &  %d", k, n);
    printf("salut,");
    CombinationList *list = init_combination_list(k, n);
    printf(" je ");
    int *combination = malloc(k * sizeof(int));
    printf("m'");
    //printf(list);
    // Initialiser la première combinaison [0, 1, 2, ..., k-1]
    for (int i = 0; i < k; i++) {
        combination[i] = i;
        //printf(combination[i]);
    }
    printf("appelle ");
    
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
    printf("Paul");
    free(combination);
    printf("Jambon");
    return list;
}

void print_combinations(CombinationList *list) { // To print combinations
    printf("Nombre de combinaisons: %lld\n", list->count);
    for (int i = 0; i < list->count; i++) {
        printf("Combinaison %d: [", i);
        for (int j = 0; j < list->combinations[i].size; j++) {
            printf("%lld", list->combinations[i].indices[j]);
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
    long long summed_state = 0;
    for(long long i = 0; i < idx; i ++) summed_state += state->occupancy[i];
    long long S_i = abs(summed_state);
    int sign_factor = pow((-1), S_i);

    if(any(state) == false){
        return make_a_vector_of_zero_state_lengthed(dim, state);
    }else if(abs(state->occupancy[idx]) == 0){
        return make_a_vector_of_zero_state_lengthed(dim, state);
    }else{
        State* ret_state = malloc(sizeof(State));
        ret_state->size = dim;
        ret_state->occupancy = malloc(ret_state->size * sizeof(int));
        ret_state->occupancy[idx]= 0;

        State* signed_ret_state = malloc(sizeof(State));
        signed_ret_state->size = dim;
        //signed_ret_state->occupancy = malloc(sizeof(int));
        signed_ret_state->occupancy = malloc(signed_ret_state->size * sizeof(int));
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
        State* ret_state = malloc(sizeof(State));
        ret_state->size = dim;
        ret_state->occupancy = malloc(ret_state->size * sizeof(int)); //ou alors juste ret_state->size * sizeof(in) ?
        for(long long i = 0; i<state->size; i++) ret_state->occupancy[i] = state->occupancy[i]; // instead of state.copy() ?
        ret_state->occupancy[idx] = 1;

        State* signed_ret_state = malloc(sizeof(State));
        signed_ret_state->size = dim;
        signed_ret_state->occupancy = malloc(signed_ret_state->size * sizeof(long long));
        for(long long i = 0; i<dim; i++) signed_ret_state->occupancy[i] = sign_factor * ret_state->occupancy[i];//dim = sizeof(state) normalement
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

// void print_pointeur(CombinationList* listcombi){
//     for(int i = 0; i<listcombi->max_count; i++){
//         printf(listcombi->combinations);
//     }
// }

// Main function :
StateList* get_hubbard_states(int N, int dim) { // get_hubbard_states(N, dim)
    //printf("test_0_ghs");
    if(dim>N) printf("jaime");
    if(dim<N) printf("nope");
    CombinationList *combs = combinations_iterative(N, dim); // dim = nombre de particule & N = nombre de site
    //printf(combs->combinations);
    //printf("test_1_ghs");
    if (!combs) return NULL;
    //printf("test_2_ghs");
    StateList *state_list = malloc(sizeof(StateList));
    //printf("test_3_ghs");
    state_list->count = combs->count;
    //printf("test_4_ghs");
    state_list->states = malloc(combs->count * sizeof(State));
    //printf(" test_5_ghs");
    
    for (long long i = 0; i < combs->count; i++) {
        // Initialiser l'état avec des zéros
        state_list->states[i].size = dim;
        //printf(" test_6_ghs");
        state_list->states[i].occupancy = calloc(dim, sizeof(int));
        //printf(" test_7_ghs");
        // Mettre 1 aux positions occupées
        for (int j = 0; j < combs->combinations[i].size; j++) {
            long long index = combs->combinations[i].indices[j];
            //printf(" test_8_ghs");
            state_list->states[i].occupancy[index] = 1;
            //printf(" test_9_ghs");
        }
    }

    if (!state_list || !state_list->states) {

        free(state_list); // au cas où state_list est non NULL

        free_combination_list(combs);
        return NULL;
    }
    //printf("laliloulalal");
    free_combination_list(combs);
    return state_list;
}


HamiltonianMatrix* hubbard_hamiltonian_matrix(int N, t_matrix* t_matrix, double U, int dim, int V){
    printf("2_Hello world !");
    StateList* statelist = get_hubbard_states(N, dim); // Get all possible Hubbard states
    printf("3_Hello world !");
    //dim = statelist->count;  // Dimension of the Hilbert space
    if (!statelist) {
        printf("Erreur: impossible de générer les états\n");
        return NULL;
    }
    int hilbert_dim = statelist->count;  // Utilisez une nouvelle variable
    
    HamiltonianMatrix* H = allocate_memory_hamiltonian(hilbert_dim);
    H = initialize_matrix_with_zeros(H, hilbert_dim);
    printf("4_Hello world !");
    
    // Loop over all states (rows)
    for(int i = 0; i < hilbert_dim; i++){
        State* state_i = &statelist->states[i];
        printf("6_Hello world !");
        // Loop over all states (columns)
        for(int j = 0; j < hilbert_dim; j++){
            State* state_j = &statelist->states[j];
            printf("7_Hello world !");
            // Diagonal elements: Coulomb interaction term 
            if(i == j){
                for (int site =0; site < N; site ++){
                    // Check if both up and down spins are present at the site
                    int n_up = number_operator(state_i, site, 'u');
                    int n_down = number_operator(state_i, site, 'd');
                    H->matrix[i][j] += U * n_up * n_down;
                    printf("poppopoopoopokok");
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
                                printf("2525252525!");
                                // Check if there is a spin to move at site1 with spin
                                if(temp != NULL){
                                    int site2 = site2_list[s];
                                    State* final = creation(temp, site2, spin, dim); // 0 if already occupied
                                    printf("creationzeeze");
                                    if(state_equal(abs_state(final), state_j)){ // problem
                                        printf("sdfg !");
                                        int sign = hopping_term_sign_factor(state_i, site1, site2, spin); // --> a voir antisymétries des fermions
                                        printf("hopping !");
                                        H->matrix[i][j] -= t_matrix->t_matrix[site1][site2] * sign; //problème ???
                                        printf("Hello_____world !");
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

//int init_binary_state = {0,1,1,0,1,0,1,0};
int top_n = 1;
int nbr_pts = 10;
int fig_width = 8;
int fig_heigth = 4;

// --> problème de mémoire : crash 3 
void top_hubbard_states_calculation(int temps, int U, t_matrix* t_matrix, State* init_binary_state, int top_n, int nbr_pts){
    U = U * eV;  // Convert U from eV to Joules
     // Convert t from eV to Joules
    printf("avantoutchose");
    int N = (init_binary_state->size) / 2;
    printf("foie");
    t_matrix->t_dim = N;
    
    t_matrix->t_matrix = malloc(t_matrix->t_dim * sizeof(t_matrix));
    for(int i = 0; i < t_matrix->t_dim; i++) {
        t_matrix->t_matrix[i] = malloc(t_matrix->t_dim * sizeof(double));
        if (t_matrix->t_matrix[i] == NULL) {
            printf("ERREUR: Allocation échouée pour ligne %d\n", i);
            // Libérer ce qui a été alloué
            for(int k = 0; k < i; k++) {
                free(t_matrix->t_matrix[k]);
            }
            free(t_matrix->t_matrix);
            return;
        }
    }
    printf("parcequ'en fait");
    for(int i = 0; i < t_matrix->t_dim; i++) {
        for(int j = 0; j < t_matrix->t_dim; j++) {
            t_matrix->t_matrix[i][j] *= eV;
        }
    }
    
    // // Number of sites
    N = (init_binary_state->size)/2; 
    int dim = init_binary_state->size;
    StateList* states = malloc(sizeof(StateList));
    states->count = dim;
    states->states = malloc(states->count * sizeof(State));
    states = get_hubbard_states(N, dim);
    //states->count = 
    dim = states->count;
    int V = 0;

    //var = ctypes.pointer(t_matrix_c)
    printf("111111111111111111!");

    HamiltonianMatrix* H = hubbard_hamiltonian_matrix(N, t_matrix, U, dim, V);
    printf("Bientot");
    // int V = 0;
    // HamiltonianMatrix* H = hubbard_hamiltonian_matrix(N, t_matrix, U, dim, V);
}

int simple(){
    return 1+1;
}

int main(){
    int V = 0;
    int N = 4;
    int dim = 2 * N;
    int temps = 4;
    double U = 2;
    // // t_matrix* t_matrix;
    // // printf("1_Hello world !");
    printf("CVBVFDFGH");
    t_matrix* t_matrix = malloc(sizeof(t_matrix));
    t_matrix->t_dim = N;
    t_matrix->t_matrix = malloc(N * sizeof(double*));
    printf("AZERETRTHJ");
    //hubbard_hamiltonian_matrix( N, t_matrix, U, dim, V);

    State* init_binary_state = malloc(sizeof(State));
    init_binary_state->size = 8;
    init_binary_state->occupancy = malloc(init_binary_state->size* sizeof(State));
    printf("japprecie ");
    int values[] = {0, 1, 1, 0, 1, 0, 1, 0};
    //printf("Taille: %d\n", init_binary_state->size);
    printf("les ");
    // Copier les valeurs
    for(int i = 0; i < init_binary_state->size; i++) {
        init_binary_state->occupancy[i] = values[i];
    }
    printf("japprecie ");
    N = (init_binary_state->size)/2; 
    dim = init_binary_state->size;
    StateList* states = malloc(sizeof(StateList));
    states = get_hubbard_states(N, dim);
    printf("le ");
    //states->count = 
    dim = states->count;
    V = 0;
    printf("OLOMOLM");
    top_hubbard_states_calculation(temps, U,  t_matrix,  init_binary_state,  top_n,  nbr_pts);
    printf("cordon bleu");
    // // HamiltonianMatrix* Hreturned = hubbard_hamiltonian_matrix(N, t_matrix, U, dim, V);
    // // printf("5_Hello world !");

    // int n = 5;
    // int k = 2;
    //combinations_iterative(k, n);
    //get_hubbard_states(N, dim);
    //printf(states->states);
    return 0;
}