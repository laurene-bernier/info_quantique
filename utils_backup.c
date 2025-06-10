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

// Fonction d'affichage d'un seul état
void print_state(State *state) { //const
    printf("État (taille = %lld) : ", state->size);
    for (long long i = 0; i < state->size; i++) {
        // Affiche l’occupation du site i (entier long long)
        printf("%lld ", state->occupancy[i]);
    }
    printf("\n");
}

// Fonction d'affichage d'une liste d'états
void print_state_list(const StateList *list) {
    printf("Liste d'états : (nombre = %lld)\n", list->count);
    for (long long i = 0; i < list->count; i++) {
        printf("État #%lld : ", i);
        print_state(&list->states[i]);
    }
}

void print_matrix(Matrix *H) {
    for (long long i = 0; i < H->dim; i++) {
        for (long long j = 0; j < H->dim; j++) {
            printf("%10.4f ", H->matrix[i][j]);
        }
        printf("\n");
    }
}

Matrix* create_tridiagonal_matrix(int N) {
    // Allocation d'un tableau de pointeurs
    Matrix* matrix = malloc(N * sizeof(Matrix));
    if (!matrix) {
        perror("malloc failed");
        exit(EXIT_FAILURE);
    }
    matrix->dim = N;
    matrix->matrix = malloc(matrix->dim * sizeof(double));
    for (int i = 0; i < N; i++) {
        matrix->matrix[i] = calloc(N, sizeof(double));  // Initialise à 0
    }
    // Remplissage de la matrice tridiagonale
    for (int i = 0; i < N; i++) {
        if (i > 0)
            matrix->matrix[i][i - 1] = 1.0; // diagonale inférieure
        if (i < N - 1)
            matrix->matrix[i][i + 1] = 1.0; // diagonale supérieure
    }
    return matrix;
}


// Definition of basic useful function :
// 1) for memory gestion :
Matrix* allocate_memory_matrix(int dim) { //function to allocate memory for a dynamic array of the Hamiltonian matrix :
    Matrix *H = malloc(sizeof(Matrix)); // = new Matrix;
    H->dim = dim;
    H->matrix = malloc(dim * sizeof(double*));
    for(int i = 0; i <dim; i++){
        H->matrix[i] = calloc(dim, sizeof(double));
    }
    return H;
}

void free_memory_matrix(Matrix* H, int dim) { //function to free memory of the Hamiltonian matrix :
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
Matrix* initialize_matrix_with_zeros(int dim){ // does what it is said in the title


    Matrix* H = allocate_memory_matrix(dim);
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

State* make_a_vector_of_zero_state_lengthed(int dim){
    State* state_nul = malloc(sizeof(State));
    state_nul->size = dim;
    state_nul->occupancy = calloc(dim, sizeof(long long));
    for(int i = 0; i<dim; i++){
        state_nul->occupancy[i] = 0;
    }
    return state_nul;
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
    if (k > n || k < 0) return NULL;
    //printf("k & n : %d &  %d", k, n);
    CombinationList *list = init_combination_list(k, n);
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

State* annihilation(State* state, int i, char spin){
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
    int idx = 0;
    if(spin == 'u') idx = 2 * i;
    else idx = 2*i + 1;


    //Calculate the Fermi sign factor
    long long summed_state = 0;

    for(long long i = 0; i < idx; i ++) summed_state += state->occupancy[i];

    long long S_i = abs(summed_state);
    int sign_factor = (S_i % 2 == 0) ? 1 : -1;

    State *new_state = malloc(sizeof(State));
    new_state->size = state->size;
    int dim = state->size;
    new_state->occupancy = malloc(state->size * sizeof(int));
    if(!any(state)){
        printf("any");
        return make_a_vector_of_zero_state_lengthed(dim);
    }else if(abs(state->occupancy[idx]) == 0){
        printf("dim : %d", dim);
        print_state(make_a_vector_of_zero_state_lengthed(dim));
        return make_a_vector_of_zero_state_lengthed(dim);
    }else{
        printf("else");
        printf("test_1");

        printf("test_2");
        for (int j = 0; j < state->size; j++) {
            new_state->occupancy[j] = state->occupancy[j];
        }
        printf("test_3");
        new_state->occupancy[idx] = 0;

        // applique le signe de Fermi
        for (int j = 0; j < state->size; j++) {
            new_state->occupancy[j] *= sign_factor;
        }

        //print_state(ret_state);
        //ret_state->occupancy[idx]= 0;

        return new_state;
    }
}

State* creation(State* state, int i, char spin){
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
    int idx = 0;
    if(spin == 'u') idx = 2 * i;
    else idx = 2*i + 1;


    //Calculate the Fermi sign factor
    long long summed_state = 0;

    for(long long i = 0; i < idx; i ++) summed_state += state->occupancy[i];

    long long S_i = abs(summed_state);
    int sign_factor = (S_i % 2 == 0) ? 1 : -1;

    State *new_state = malloc(sizeof(State));
    new_state->size = state->size;
    int dim = state->size;
    new_state->occupancy = malloc(state->size * sizeof(int));
    if(!any(state)){
        printf("any");
        return make_a_vector_of_zero_state_lengthed(dim);
    }else if(abs(state->occupancy[idx]) == 1){
        printf("dim : %d", dim);
        print_state(make_a_vector_of_zero_state_lengthed(dim));
        return make_a_vector_of_zero_state_lengthed(dim);
    }else{
        for (int j = 0; j < state->size; j++) {
            new_state->occupancy[j] = state->occupancy[j];
        }
        new_state->occupancy[idx] = 1;

        // applique le signe de Fermi
        for (int j = 0; j < state->size; j++) {
            new_state->occupancy[j] *= sign_factor;
        }
        return new_state;
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
StateList* get_hubbard_states(int N) { // get_hubbard_states(N, dim
    int dim = 2 * N;
    CombinationList *combs = combinations_iterative(N, dim); // dim = nombre de particule & N = nombre de site
    if (!combs) {
        free_combination_list(combs);
        return NULL;
    }
    StateList *state_list = malloc(sizeof(StateList));
    state_list->count = combs->count;
    state_list->states = malloc(combs->count * sizeof(State));

    for (long long i = 0; i < combs->count; i++) {
        // Initialiser l'état avec des zéros
        state_list->states[i].size = dim;
        state_list->states[i].occupancy = calloc(dim, sizeof(long long));
        // Mettre 1 aux positions occupées
        for (int j = 0; j < combs->combinations[i].size; j++) {
            long long index = combs->combinations[i].indices[j];
            state_list->states[i].occupancy[index] = 1;
        }
    }

    if (!state_list || !state_list->states) {
        free(state_list); // au cas où state_list est non NULL
        free_combination_list(combs);
        return NULL;
    }

    free_combination_list(combs);
    return state_list;
}


Matrix* hubbard_hamiltonian_matrix(int N, Matrix* t_matrix, double U, int dim, int V){
    StateList* statelist = get_hubbard_states(N); // Get all possible Hubbard states
    //dim = statelist->count;  // Dimension of the Hilbert space
    if (!statelist) {
        printf("Erreur: impossible de générer les états\n");
        free_state_list(statelist);
        return NULL;
    }
    int hilbert_dim = (int)(statelist->count);  // Utilisez une nouvelle variable

    Matrix* H = initialize_matrix_with_zeros(hilbert_dim);
    // Loop over all states (rows)
    for(int i = 0; i < hilbert_dim; i++){
        printf("i=%d ", i);
        State* state_i = &statelist->states[i];
        // Loop over all states (columns)
        for(int j = 0; j < hilbert_dim; j++){
            printf("j=%d ", j);
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
                printf("hopping");
                // Determine if statelist i and j differ by a single hopping event
                for(int site1 = 0; site1 < N; site1++){
                    // Hubbard nearest-neighbor hopping
                    int site2_list[2] = {site1-1, site1+1};
                    for(int s = 0; s < 2; s++){
                        if(0 <= site2_list[s] && site2_list[s] < N){
                            char spins[2] = {'u', 'd'};
                            for (int l = 0; l < 2; l++) {
                                char spin = spins[l];

                                // State* temp = annihilation(state_i, site1, spin);
                                // // Check if there is a spin to move at site1 with spin
                                // if(any(temp)){
                                //     int site2 = site2_list[s];

                                //     State* final = creation(temp, site2, spin, dim); // 0 if already occupie

                                //     printf("final :");
                                //     print_state(final);
                                //     if(state_equal(abs_state(final), state_j)){

                                //         printf("sdfg !");
                                //         int sign = hopping_term_sign_factor(state_i, site1, site2, spin); // --> a voir antisymétries des fermions
                                //         printf("hopping !");
                                //         H->matrix[i][j] -= t_matrix->matrix[site1][site2] * sign; //problème ???
                                //         printf("Hello___world !");
                                //     }
                                //}
                            }
                        }
                    }
                }
            }
            printf("j = %d", j);
        }
    }
    printf("test_111");
    free_state_list(statelist);
    return H;
}

//int init_binary_state = {0,1,1,0,1,0,1,0};
int top_n = 1;
int nbr_pts = 10;

void top_hubbard_states_calculation(int temps, int U, Matrix* t_matrix, State* init_binary_state, int nbr_pts){
    U = U * eV;  // Convert U from eV to Joules
    // Convert t from eV to Joules
    printf("avantoutchose");
    int N = (init_binary_state->size) / 2;
    printf("foie");
    t_matrix->dim = N;

    t_matrix->matrix = malloc(t_matrix->dim * sizeof(t_matrix));
    for(int i = 0; i < t_matrix->dim; i++) {
        t_matrix->matrix[i] = malloc(t_matrix->dim * sizeof(double));
        if (t_matrix->matrix[i] == NULL) {
            printf("ERREUR: Allocation échouée pour ligne %d\n", i);
            // Libérer ce qui a été alloué
            for(int k = 0; k < i; k++) {
                free(t_matrix->matrix[k]);
            }
            free(t_matrix->matrix);
            return;
        }
    }
    printf("parcequ'en fait");
    for(int i = 0; i < t_matrix->dim; i++) {
        for(int j = 0; j < t_matrix->dim; j++) {
            t_matrix->matrix[i][j] *= eV;
        }
    }

    // // Number of sites
    N = (init_binary_state->size)/2;
    int dim = init_binary_state->size;
    StateList* states = malloc(sizeof(StateList));
    states->count = dim;
    states->states = malloc(states->count * sizeof(State));
    states = get_hubbard_states(N);
    //states->count =
    dim = states->count;
    int V = 0;

    Matrix* H = hubbard_hamiltonian_matrix(N, t_matrix, U, dim, V);
    printf("Bientot");

    // if(H){
    //     free_memory_matrix(H, dim); // free memory
    // }
    free_state_list(states);
    // int V = 0;
    // Matrix* H = hubbard_hamiltonian_matrix(N, t_matrix, U, dim, V);
}

int main(){
    int V = 0;
    int N = 2;
    //int dim = 2 * N;
    int temps = 4;
    double U = 2;
    int dim = 2 * N;

    StateList* statelist = malloc(sizeof(StateList));
    statelist = get_hubbard_states(N); // segmentation fault
    //print_state_list(statelist);

    Matrix* t_matrix = create_tridiagonal_matrix(N);

    State state = statelist->states[0];

    //print_state(&state);

    //State* state_nul;
    //state_nul = make_a_vector_of_zero_state_lengthed(3);
    //print_state(state_nul);
    free_state_list(statelist);

    //State* state_anni;
    //state_anni = creation(&state,  0,  'u');

    //Matrix* H = hubbard_hamiltonian_matrix(N, t_matrix, U, dim, V); // in process
    
    //print_matrix(H);

    State* a;
    State* b;
    a->occupancy = malloc(sizeof(long long));
    b->occupancy = malloc(sizeof(long long));
    a->size = b->size;
    bool bolen = state_equal(a, b);
    printf("bool : %d", bolen);
    //print_state(state_anni);

    //free_state_list(statelist);
    return 0;
}
