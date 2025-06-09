#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


// Structure declaration :

// Structure for t_matrix
typedef struct {
    long long **t_matrix;
    long long t_dim;
} t_matrix;

// Structure pour la matrice Hamiltonienne
typedef struct {
    double **matrix;
    long long dim;
} HamiltonianMatrix;

// Structure pour représenter un état quantique
typedef struct {
    long long *occupancy;  // Tableau d'occupation des sites (bits pour spin up/down)
    long long size;        // Nombre de sites
} State;

// 3 structures pour générer une combinaison :
// Structure to stock a state list
typedef struct {
    State *states;
    long long count;
} StateList;

// Structure to stock one state combination
typedef struct {
    long long *indices;
    long long size;
} Combination;

// Structure to stock all states combinations
typedef struct {
    Combination *combinations;
    long long count;
    long long max_count;
} CombinationList;



// Function declaration (in order) :

// Definition of basic useful function :
// 1) for memory gestion :
HamiltonianMatrix* allocate_memory_hamiltonian(int dim); // done
void free_memory_hamiltonian(HamiltonianMatrix* H, int dim); // done
void free_combination_list(CombinationList *list);
void free_state_list(StateList *list);


// 2) basic :
HamiltonianMatrix* initialize_matrix_with_zeros(HamiltonianMatrix* H, int dim); // done
bool state_equal(State* state1, State* state2); // done
State* abs_state(State* not_abs_state); // done
//int power(int minus_one_sign, int S_i); // done
bool any(State* state); // done
State* make_a_vector_of_zero_state_lengthed(int dim, State* state); // done


// 3) for the denombrement in get_hubard_state :
int binomial_coefficient(int n, int k); // done
CombinationList* init_combination_list(int n, int k); // done
CombinationList* combinations_iterative(int n, int k); // done
//CombinationList* combinations(int n, int k);
void print_combinations(CombinationList *list); // utile ?
StateList* generate_hubbard_states(int dim, int N);


// Other specific function :
int number_operator(State* state, int i, char spin); // in process//
State* annihilation(State* state, int i, char spin, int dim); // done
State* creation(State* state, int i, char spin, int dim); // done
int hopping_term_sign_factor(State* state_i, int i, int k, char spin);

// Main function :

// to generate a hamiltonian matrix :
HamiltonianMatrix* hubbard_hamiltonian_matrix(int N, t_matrix* t_matrix, double U, int dim, int V); // in process
//HamiltonianMatrix* hubbard_hamiltonian_matrix(int N, t_matrix* t_matrix, double U, StateList* statelist);
StateList* get_hubbard_states(int N, int dim); // in process



//State* get_hubbard_states(int N, int *num_states);


//int states_equal(State state1, State state2);

//double hopping_term_sign_factor(State state, int site1, int site2, char spin); // to do


