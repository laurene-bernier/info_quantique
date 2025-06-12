#ifndef MAIN_H
#define MAIN_H

#include <complex.h>
#include <stdbool.h>

typedef double complex cplx;

/* --- États quantiques --- */
typedef struct {
    long long    size;       // nombre de cases dans occupancy
    long long   *occupancy;  // tableau d’occupation 0/1
} State;

typedef struct {
    int size;       // Taille de l'état
    cplx* vector;   // Tableau de coefficients complexes
} StateComplexe;

typedef struct {
    State      *states;
    long long   count;
} StateList;

typedef struct {
    StateComplexe* complexe_state;
    long long count;
} StateListComplexe;

typedef struct {
    int dim;
    cplx **compMatrix;
} ComplexeMatrix;


/* --- Combinaisons pour get_hubbard_states --- */
typedef struct {
    long long *indices;
    long long  size;
} Combination;

typedef struct {
    Combination *combinations;
    long long    count;
    long long    max_count;
} CombinationList;

/* --- Matrices --- */
typedef struct {
    double    **matrix;
    long long   dim;
} Matrix;

typedef double complex cplx;

/* --- Prototypes --- */
void            print_state(const State *s);
void            print_state_list(const StateList *L);
void            print_matrix(const Matrix *M);

Matrix*         allocate_memory_matrix(long long dim);
Matrix*         initialize_matrix_with_zeros(long long dim);
void            free_memory_matrix(Matrix *M);
int hopping_term_sign_factor(const State *state_i,
                             int i,
                             int k,
                             char spin);
void print_state_complexe(StateComplexe* state);

CombinationList* combinations_iterative(int k, int n);
void            free_combination_list(CombinationList *L);

StateList*      get_hubbard_states(int N);
void            free_state_list(StateList *L);
void print_state_list_complexe(StateListComplexe* list);

ComplexeMatrix* time_evol_operator(Matrix* H, double t);
StateListComplexe* time_evol_state(Matrix* H, double* T_array, int nbr_pts, StateComplexe* u);
double* transition_probability_over_time(StateComplexe* left_state, StateListComplexe * list);

int             number_operator(const State *s, int site, char spin);
State*          annihilation(const State *s, int site, char spin);
State*          creation    (const State *s, int site, char spin);

Matrix*         hubbard_hamiltonian_matrix(int N, Matrix *tmat, double U);
void top_hubbard_states_interface(int N, double U, double T_final, int nbr_pts, int top_n, State* init_state, const char* filename);
#endif // MAIN_H