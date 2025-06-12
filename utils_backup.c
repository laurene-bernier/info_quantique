#include "utils_backup.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#define MAX_TERMS 100
#define TOL 1e-12

static const double eV = 1.602176634e-19;

void print_state(const State *s) {
    printf("État (taille = %lld) : ", s->size);
    for (long long i = 0; i < s->size; i++) {
        printf("%lld ", s->occupancy[i]);
    }
    printf("\n");
}

Matrix* allocate_memory_matrix(long long dim) {
    Matrix *H = malloc(sizeof *H);
    H->dim    = dim;
    H->matrix = malloc(dim * sizeof *H->matrix);
    for (long long i = 0; i < dim; i++) {
        H->matrix[i] = calloc(dim, sizeof *H->matrix[i]);
    }
    return H;
}

void free_memory_matrix(Matrix *H) {
    if (!H) return;
    for (long long i = 0; i < H->dim; i++)
        free(H->matrix[i]);
    free(H->matrix);
    free(H);
}


StateComplexe* basis_vector(int dim, int idx) {
    StateComplexe* sc = malloc(sizeof(StateComplexe));
    sc->size = dim;
    sc->vector = calloc(dim, sizeof(cplx));
    sc->vector[idx] = 1.0 + 0.0*I;
    return sc;
}

static long long binomial_coefficient(int n, int k) {
    if (n < 0 || k < 0) return -1;
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    long long res = 1;
    for (int i = 1; i <= k; i++)
        res = res * (n - k + i) / i;
    return res;
}

CombinationList* init_combination_list(int k, int n) {
    CombinationList *L = malloc(sizeof *L);
    L->count     = 0;
    L->max_count = binomial_coefficient(n, k);
    L->combinations = malloc(L->max_count * sizeof *L->combinations);
    for (long long i = 0; i < L->max_count; i++) {
        L->combinations[i].size    = k;
        L->combinations[i].indices = malloc(k * sizeof *L->combinations[i].indices);
    }
    return L;
}

CombinationList* combinations_iterative(int k, int n) {
    if (k > n || k < 0) return NULL;
    CombinationList *L = init_combination_list(k, n);
    long long *comb = malloc(k * sizeof *comb);
    for (int i = 0; i < k; i++) comb[i] = i;
    do {
        memcpy(L->combinations[L->count].indices, comb, k * sizeof *comb);
        L->count++;
        int pos = k - 1;
        while (pos >= 0 && comb[pos] == n - k + pos) pos--;
        if (pos < 0) break;
        comb[pos]++;
        for (int j = pos + 1; j < k; j++)
            comb[j] = comb[j-1] + 1;
    } while (1);
    free(comb);
    return L;
}

void free_combination_list(CombinationList *L) {
    if (!L) return;
    for (long long i = 0; i < L->max_count; i++)
        free(L->combinations[i].indices);
    free(L->combinations);
    free(L);
}

bool any(const State *s) {
    for (long long i = 0; i < s->size; i++)
        if (s->occupancy[i] != 0) return 1;
    return 0;
}

State* make_zero_state(long long dim) {
    State *Z = malloc(sizeof *Z);
    Z->size      = dim;
    Z->occupancy = calloc(dim, sizeof *Z->occupancy);
    return Z;
}

State* abs_state(const State *src) {
    State *A = malloc(sizeof *A);
    A->size      = src->size;
    A->occupancy = malloc(src->size * sizeof *A->occupancy);
    for (long long i = 0; i < src->size; i++)
        A->occupancy[i] = llabs(src->occupancy[i]);
    return A;
}

bool state_equal(const State *a, const State *b) {
    if (a->size != b->size) return 0;
    for (long long i = 0; i < a->size; i++)
        if (a->occupancy[i] != b->occupancy[i]) return 0;
    return 1;
}

int number_operator(const State *s, int site, char spin) {
    return spin=='u'
           ?  (int)s->occupancy[2*site]
           :  (int)s->occupancy[2*site+1];
}

State* annihilation(const State *st, int site, char spin) {
    long long idx = (spin=='u' ? 2*site : 2*site+1);
    long long sum = 0;
    for (long long i = 0; i < idx; i++) sum += st->occupancy[i];
    int sign = (sum % 2 == 0) ? +1 : -1;

    if (!any(st) || st->occupancy[idx] == 0)
        return make_zero_state(st->size);

    State *out = malloc(sizeof *out);
    out->size      = st->size;
    out->occupancy = malloc(st->size * sizeof *out->occupancy);
    memcpy(out->occupancy, st->occupancy, st->size * sizeof *out->occupancy);
    out->occupancy[idx] = 0;
    for (long long i = 0; i < out->size; i++)
        out->occupancy[i] *= sign;
    return out;
}

State* creation(const State *st, int site, char spin) {
    long long idx = (spin=='u' ? 2*site : 2*site+1);
    long long sum = 0;
    for (long long i = 0; i < idx; i++) sum += st->occupancy[i];
    int sign = (sum % 2 == 0) ? +1 : -1;

    if (!any(st) || st->occupancy[idx] == 1)
        return make_zero_state(st->size);

    State *out = malloc(sizeof *out);
    out->size      = st->size;
    out->occupancy = malloc(st->size * sizeof *out->occupancy);
    memcpy(out->occupancy, st->occupancy, st->size * sizeof *out->occupancy);
    out->occupancy[idx] = 1;
    for (long long i = 0; i < out->size; i++)
        out->occupancy[i] *= sign;
    return out;
}

double complex hermitian_dot(const double complex* a, const double complex* b, int size) {
    double complex result = 0.0 + 0.0 * I;
    for (int i = 0; i < size; i++)
        result += conj(a[i]) * b[i];
    return result;
}

int hbar = 1;

ComplexeMatrix* allocate_complex_matrix(int dim) {
    ComplexeMatrix* m = malloc(sizeof(ComplexeMatrix));
    m->dim = dim;
    m->compMatrix = malloc(dim * sizeof(cplx*));
    for (int i = 0; i < dim; i++)
        m->compMatrix[i] = calloc(dim, sizeof(cplx));
    return m;
}

ComplexeMatrix* mat_identity(int dim) {
    ComplexeMatrix* out = allocate_complex_matrix(dim);
    for (int i = 0; i < out->dim; i++)
        for (int j = 0; j < out->dim; j++)
            out->compMatrix[i][j] = (i == j) ? 1.0 + 0.0*I : 0.0 + 0.0*I;
    return out;
}

ComplexeMatrix* mat_add(ComplexeMatrix* a, ComplexeMatrix* b) {
    if (!a || !b || a->dim != b->dim) return NULL;
    int dim = a->dim;
    ComplexeMatrix* out = allocate_complex_matrix(dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            out->compMatrix[i][j] = a->compMatrix[i][j] + b->compMatrix[i][j];
    return out;
}

ComplexeMatrix* mat_mul(ComplexeMatrix* a, ComplexeMatrix* b) {
    if (!a || !b || a->dim != b->dim) return NULL;
    int dim = a->dim;
    ComplexeMatrix* out = allocate_complex_matrix(dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            for (int k = 0; k < dim; k++)
                out->compMatrix[i][j] += a->compMatrix[i][k] * b->compMatrix[k][j];
    return out;
}

ComplexeMatrix* mat_scalar_div(ComplexeMatrix* m, double scalar) {
    for (int i = 0; i < m->dim; i++)
        for (int j = 0; j < m->dim; j++)
            m->compMatrix[i][j] /= scalar;
    return m;
}

void free_complex_matrix(ComplexeMatrix* matrix) {
    if (!matrix) return;
    for (int i = 0; i < matrix->dim; i++)
        free(matrix->compMatrix[i]);
    free(matrix->compMatrix);
    free(matrix);
}

void generate_time_array(double T_final, int nbr_pts, double hbar, double* T_array) {
    double dt = T_final / (nbr_pts - 1);
    for (int i = 0; i < nbr_pts; i++)
        T_array[i] = (i * dt) / hbar;
}

ComplexeMatrix* matrix_exponential(ComplexeMatrix* A) {
    if (!A) return NULL;
    int dim = A->dim;
    ComplexeMatrix* result = mat_identity(dim);
    ComplexeMatrix* term = mat_identity(dim);
    ComplexeMatrix* temp = NULL;

    for (int k = 1; k < MAX_TERMS; k++) {
        temp = mat_mul(term, A);
        free_complex_matrix(term);
        term = temp;
        mat_scalar_div(term, k);

        ComplexeMatrix* next_result = mat_add(result, term);
        free_complex_matrix(result);
        result = next_result;

        double norm = 0.0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                norm += cabs(term->compMatrix[i][j]);
        if (norm < TOL) break;
    }
    free_complex_matrix(term);
    return result;
}

ComplexeMatrix* time_evol_operator(Matrix* H, double t) {
    ComplexeMatrix* A = malloc(sizeof(ComplexeMatrix));
    A->dim = H->dim;
    A->compMatrix = malloc(A->dim * sizeof(cplx*));
    for (int i = 0; i < A->dim; i++)
        A->compMatrix[i] = malloc(A->dim * sizeof(cplx));
    for (int i = 0; i < A->dim; i++)
        for (int j = 0; j < A->dim; j++)
            A->compMatrix[i][j] = -I * H->matrix[i][j] * t / hbar;
    ComplexeMatrix* U = matrix_exponential(A);
    free_complex_matrix(A);
    return U;
}

StateListComplexe* time_evol_state(Matrix* H, double* T_array, int nbr_pts, StateComplexe* u) {
    StateListComplexe* result_array = malloc(sizeof(StateListComplexe));
    result_array->count = nbr_pts;
    result_array->complexe_state = malloc(nbr_pts * sizeof(StateComplexe));
    for (int i = 0; i < nbr_pts; i++) {
        result_array->complexe_state[i].size = u->size;
        result_array->complexe_state[i].vector = calloc(u->size, sizeof(cplx));
    }
    ComplexeMatrix* U;
    for (int t_idx = 0; t_idx < nbr_pts; t_idx++) {
        double t = T_array[t_idx];
        U = time_evol_operator(H, t);
        for (int i = 0; i < U->dim; i++)
            for (int j = 0; j < U->dim; j++)
                result_array->complexe_state[t_idx].vector[i] += U->compMatrix[i][j] * u->vector[j];
        free_complex_matrix(U);
    }
    return result_array;
}

StateList* get_hubbard_states(int N) {
    long long dim = 2LL * N;
    CombinationList *C = combinations_iterative(N, dim);
    if (!C) return NULL;
    StateList *L = malloc(sizeof *L);
    L->count  = C->count;
    L->states = malloc(L->count * sizeof *L->states);
    for (long long i = 0; i < C->count; i++) {
        L->states[i].size      = dim;
        L->states[i].occupancy = calloc(dim, sizeof *L->states[i].occupancy);
        for (long long k = 0; k < C->combinations[i].size; k++) {
            long long idx = C->combinations[i].indices[k];
            L->states[i].occupancy[idx] = 1;
        }
    }
    free_combination_list(C);
    return L;
}

void free_state_list(StateList *L) {
    if (!L) return;
    for (long long i = 0; i < L->count; i++)
        free(L->states[i].occupancy);
    free(L->states);
    free(L);
}

Matrix* initialize_matrix_with_zeros(long long dim) {
    return allocate_memory_matrix(dim);
}

void print_matrix(const Matrix *M) {
    for (long long i = 0; i < M->dim; i++) {
        for (long long j = 0; j < M->dim; j++) {
            printf("%10.4f ", M->matrix[i][j]);
        }
        printf("\n");
    }
}

int hopping_term_sign_factor(const State *state_i, int i, int k, char spin) {
    int idx_i = (spin=='u' ? 2*i   : 2*i+1);
    int idx_k = (spin=='u' ? 2*k   : 2*k+1);
    int min_idx = idx_i < idx_k ? idx_i : idx_k;
    int max_idx = idx_i < idx_k ? idx_k : idx_i;
    long long S = 0;
    for (int j = min_idx+1; j < max_idx; j++)
        S += state_i->occupancy[j];
    return (S % 2) ? -1 : +1;
}

Matrix* create_tridiagonal_matrix(long long N) {
    Matrix *M = allocate_memory_matrix(N);
    if (!M) {
        perror("Erreur d’allocation");
        exit(EXIT_FAILURE);
    }
    for (long long i = 0; i < N; i++) {
        if (i > 0)     M->matrix[i][i-1] = 1.0;
        if (i < N-1)   M->matrix[i][i+1] = 1.0;
    }
    return M;
}

Matrix* hubbard_hamiltonian_matrix(int N, Matrix* t_matrix, double U){
    StateList* statelist = get_hubbard_states(N);
    if (!statelist) {
        free_state_list(statelist);
        return NULL;
    }
    int hilbert_dim = (int)(statelist->count);
    Matrix* H = initialize_matrix_with_zeros(hilbert_dim);
    for(int i = 0; i < hilbert_dim; i++){
        State* state_i = &statelist->states[i];
        for(int j = 0; j < hilbert_dim; j++){
            State* state_j = &statelist->states[j];
            if(i == j){
                for (int site =0; site < N; site ++){
                    int n_up = number_operator(state_i, site, 'u');
                    int n_down = number_operator(state_i, site, 'd');
                    H->matrix[i][j] += U * n_up * n_down;
                }
            } else{
                for(int site1 = 0; site1 < N; site1++){
                    char spins[2] = {'u', 'd'};
                    for (int l = 0; l < 2; l++) {
                        char spin = spins[l];
                        State* temp = annihilation(state_i, site1, spin);
                        if(any(temp)){
                            int site2_list[2] = {site1-1, site1+1};
                            for(int s = 0; s < 2; s++){
                                if(0 <= site2_list[s] && site2_list[s] < N){
                                    int site2 = site2_list[s];
                                    State* final = creation(temp, site2, spin);
                                    if(state_equal(abs_state(final), state_j)){
                                        int sign = hopping_term_sign_factor(state_i, site1, site2, spin);
                                        H->matrix[i][j] -= t_matrix->matrix[site1][site2] * sign;
                                    }
                                    free(final->occupancy);
                                    free(final);
                                }
                            }
                        }
                        free(temp->occupancy);
                        free(temp);
                    }
                }
            }
        }
    }
    free_state_list(statelist);
    return H;
}

StateComplexe* convert_state_to_complex(State* s) {
    StateComplexe* sc = malloc(sizeof(StateComplexe));
    sc->size = s->size;
    sc->vector = malloc(sc->size * sizeof(cplx));
    for (int i = 0; i < sc->size; i++)
        sc->vector[i] = (cplx)(s->occupancy[i]) + 0.0 * I;
    return sc;
}

void free_state_complexe(StateComplexe* sc) {
    if (sc) {
        free(sc->vector);
        free(sc);
    }
}

void print_state_list_complexe(StateListComplexe* list) {
    if (!list) {
        printf("Liste vide (NULL).\n");
        return;
    }
    for (long long i = 0; i < list->count; i++) {
        StateComplexe state = list->complexe_state[i];
        printf("État #%lld : [", i);
        for (int j = 0; j < state.size; j++) {
            double complex z = state.vector[j];
            printf(" %.3f%+.3fi ", creal(z), cimag(z));
            if (j < state.size - 1) printf("|");
        }
        printf("]\n");
    }
}

void print_state_complexe(StateComplexe* state) {
    if (!state || !state->vector) {
        printf("État complexe vide ou NULL.\n");
        return;
    }
    printf("StateComplexe (taille = %d):\n[", state->size);
    for (int i = 0; i < state->size; i++) {
        double complex z = state->vector[i];
        printf(" %.3f%+.3fi ", creal(z), cimag(z));
        if (i < state->size - 1) printf("|");
    }
    printf("]\n");
}

double* transition_probability_over_time(StateComplexe* left_state, StateListComplexe* list) {
    // Alloue le tableau résultat
    double* probabilities = malloc(list->count * sizeof(double));
    if (!probabilities) return NULL;

    for (long long t_idx = 0; t_idx < list->count; t_idx++) {
        cplx inner_product = 0.0 + 0.0 * I;
        StateComplexe* right_state = &list->complexe_state[t_idx];
        if (right_state->size != left_state->size) {
            fprintf(stderr, "Erreur : tailles d'états incompatibles à t = %lld\n", t_idx);
            probabilities[t_idx] = 0.0;
            continue;
        }
        for (int i = 0; i < left_state->size; i++) {
            inner_product += conj(left_state->vector[i]) * right_state->vector[i];
        }
        probabilities[t_idx] = pow(cabs(inner_product), 2);
    }
    return probabilities;
}

void print_double_array(const double* arr, long long size) {
    printf("[");
    for (long long i = 0; i < size; i++) {
        printf("%.6f", arr[i]);
        if (i < size - 1) printf(", ");
    }
    printf("]\n");
}

void save_top_hubbard_states_to_csv(
        StateListComplexe* psi_t,
        StateList* state_list_init,
        int top_n,
        const char* filename)
{
    FILE* f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "Impossible d'ouvrir %s\n", filename);
        return;
    }

    // Header
    fprintf(f, "t_idx,idx,proba,real,imag,occupation\n");

    for (int t = 0; t < psi_t->count; t++) {
        // Calculer la norme totale à ce t pour normalisation éventuelle
        double norm = 0.0;
        for (int j = 0; j < psi_t->complexe_state[t].size; j++)
            norm += pow(cabs(psi_t->complexe_state[t].vector[j]), 2);

        for (int k = 0; k < top_n; k++) {
            // Cherche le top-n en amplitude (naïf, à améliorer si besoin)
            int idx = k;
            double amp = pow(cabs(psi_t->complexe_state[t].vector[idx]), 2);
            double re  = creal(psi_t->complexe_state[t].vector[idx]);
            double im  = cimag(psi_t->complexe_state[t].vector[idx]);

            // État d'occupation
            fprintf(f, "%d,%d,%.12f,%.12f,%.12f,[", t, idx, amp, re, im);
            for (int s = 0; s < state_list_init->states[idx].size; s++) {
                fprintf(f, "%lld", state_list_init->states[idx].occupancy[s]);
                if (s < state_list_init->states[idx].size - 1)
                    fprintf(f, " ");
            }
            fprintf(f, "]\n");
        }
    }
    fclose(f);
}

void top_hubbard_states_interface(int N, double U, double T_final, int nbr_pts, int top_n, State* init_state, const char* filename) {
    // 1. Génére la base
    StateList* state_list = get_hubbard_states(N);

    // 2. Hamiltonien
    Matrix *tmat = create_tridiagonal_matrix(N);
    Matrix *H    = hubbard_hamiltonian_matrix(N, tmat, U);

    // 3. État initial : choisis-le ou fais-le passer en argument si besoin
    //State* state_init = &state_list->states[1]; // ou let user choose

    // 4. Évolution temporelle
    double* T_array = malloc(nbr_pts * sizeof(double));
    generate_time_array(T_final, nbr_pts, 1.0, T_array);

    StateListComplexe* psi_t = time_evol_state(H, T_array, nbr_pts, convert_state_to_complex(init_state));

    // 5. Appelle l’export (ex : .csv ou .txt, ou écris un format numpy lisible, sinon tu ajoutes l’export .npz)
    save_top_hubbard_states_to_csv(psi_t, state_list, top_n, filename);

    // 6. Free
    free(T_array);
    free_memory_matrix(tmat);
    free_memory_matrix(H);
    free_state_list(state_list);
    // etc.
}

int main(void) {
    int N = 4;
    double U = 2.0;
    int nbr_pts = 100;
    double T_final = 1e-20;

    StateList* state_list = get_hubbard_states(N);
    int dim = state_list->count;

    // --- Recherche de l'indice correspondant à [down up] ---
    int init_idx = 40;

    StateComplexe* v0 = basis_vector(dim, init_idx);

    Matrix* tmat = create_tridiagonal_matrix(N);
    Matrix* H = hubbard_hamiltonian_matrix(N, tmat, U);

    double* T_array = malloc(nbr_pts * sizeof(double));
    generate_time_array(T_final, nbr_pts, 1.0, T_array);

    StateListComplexe* psi_t = time_evol_state(H, T_array, nbr_pts, v0);

    double* proba = transition_probability_over_time(v0, psi_t);

    printf("[");
    for (int i = 0; i < nbr_pts; i++) {
        printf("%f", proba[i]);
        if (i < nbr_pts - 1) printf(", ");
    }
    printf("]\n");


    int top_n = 4; // Ou ce que tu veux
    save_top_hubbard_states_to_csv(psi_t, state_list, top_n, "top_hubbard_states.csv");

    // --- Free mémoire ---
    free_state_list(state_list);
    free_state_complexe(v0);
    free(T_array);
    free_memory_matrix(tmat);
    free_memory_matrix(H);
    for (int i = 0; i < psi_t->count; i++)
        free(psi_t->complexe_state[i].vector);
    free(psi_t->complexe_state);
    free(psi_t);
    free(proba);

    return 0;
}