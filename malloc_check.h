#ifndef MALLOC_CHECK_H
#define MALLOC_CHECK_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// =============================================================================
// OPTION 1: MACROS SIMPLES POUR VÉRIFICATION MALLOC
// =============================================================================

// Macro pour malloc avec vérification
#define SAFE_MALLOC(ptr, size, type) do { \
    ptr = (type*)malloc(size); \
    if (ptr == NULL) { \
        fprintf(stderr, "ERREUR MALLOC: Échec allocation de %zu octets à la ligne %d dans %s\n", \
                (size_t)(size), __LINE__, __FILE__); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

// Macro pour calloc avec vérification
#define SAFE_CALLOC(ptr, count, size, type) do { \
    ptr = (type*)calloc(count, size); \
    if (ptr == NULL) { \
        fprintf(stderr, "ERREUR CALLOC: Échec allocation de %zu éléments de %zu octets à la ligne %d dans %s\n", \
                (size_t)(count), (size_t)(size), __LINE__, __FILE__); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

// Macro pour realloc avec vérification
#define SAFE_REALLOC(ptr, size, type) do { \
    type* temp = (type*)realloc(ptr, size); \
    if (temp == NULL && size != 0) { \
        fprintf(stderr, "ERREUR REALLOC: Échec réallocation de %zu octets à la ligne %d dans %s\n", \
                (size_t)(size), __LINE__, __FILE__); \
        free(ptr); \
        exit(EXIT_FAILURE); \
    } \
    ptr = temp; \
} while(0)

// =============================================================================
// OPTION 2: FONCTIONS AVEC GESTION D'ERREUR PERSONNALISÉE
// =============================================================================

// Fonction malloc sécurisée avec message d'erreur personnalisé
static inline void* safe_malloc_func(size_t size, const char* var_name, const char* file, int line) {
    void* ptr = malloc(size);
    if (ptr == NULL) {
        fprintf(stderr, "ERREUR MALLOC: Impossible d'allouer %zu octets pour '%s' dans %s:%d\n", 
                size, var_name, file, line);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

// Fonction calloc sécurisée
static inline void* safe_calloc_func(size_t count, size_t size, const char* var_name, const char* file, int line) {
    void* ptr = calloc(count, size);
    if (ptr == NULL) {
        fprintf(stderr, "ERREUR CALLOC: Impossible d'allouer %zu éléments de %zu octets pour '%s' dans %s:%d\n", 
                count, size, var_name, file, line);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

// Macros utilisant les fonctions ci-dessus
#define MALLOC_CHECK(var, size) \
    var = safe_malloc_func(size, #var, __FILE__, __LINE__)

#define CALLOC_CHECK(var, count, size) \
    var = safe_calloc_func(count, size, #var, __FILE__, __LINE__)

// =============================================================================
// OPTION 3: SYSTÈME DE NETTOYAGE AUTOMATIQUE EN CAS D'ERREUR
// =============================================================================

// Structure pour gérer une liste de pointeurs à libérer en cas d'erreur
typedef struct {
    void** ptrs;
    size_t count;
    size_t capacity;
} CleanupList;

// Initialiser la liste de nettoyage
static inline CleanupList* init_cleanup_list(void) {
    CleanupList* list = malloc(sizeof(CleanupList));
    if (!list) return NULL;
    
    list->ptrs = malloc(10 * sizeof(void*));
    if (!list->ptrs) {
        free(list);
        return NULL;
    }
    
    list->count = 0;
    list->capacity = 10;
    return list;
}

// Ajouter un pointeur à la liste de nettoyage
static inline int add_to_cleanup(CleanupList* list, void* ptr) {
    if (!list || !ptr) return -1;
    
    if (list->count >= list->capacity) {
        list->capacity *= 2;
        void** temp = realloc(list->ptrs, list->capacity * sizeof(void*));
        if (!temp) return -1;
        list->ptrs = temp;
    }
    
    list->ptrs[list->count++] = ptr;
    return 0;
}

// Nettoyer tous les pointeurs et libérer la liste
static inline void cleanup_all(CleanupList* list) {
    if (!list) return;
    
    for (size_t i = 0; i < list->count; i++) {
        if (list->ptrs[i]) {
            free(list->ptrs[i]);
            list->ptrs[i] = NULL;
        }
    }
    
    free(list->ptrs);
    free(list);
}

// Malloc avec ajout automatique à la liste de nettoyage
static inline void* malloc_with_cleanup(size_t size, CleanupList* cleanup, const char* var_name, const char* file, int line) {
    void* ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "ERREUR MALLOC: Impossible d'allouer %zu octets pour '%s' dans %s:%d\n", 
                size, var_name, file, line);
        cleanup_all(cleanup);
        exit(EXIT_FAILURE);
    }
    
    if (add_to_cleanup(cleanup, ptr) != 0) {
        fprintf(stderr, "ERREUR: Impossible d'ajouter le pointeur à la liste de nettoyage\n");
        free(ptr);
        cleanup_all(cleanup);
        exit(EXIT_FAILURE);
    }
    
    return ptr;
}

#define MALLOC_CLEANUP(var, size, cleanup_list) \
    var = malloc_with_cleanup(size, cleanup_list, #var, __FILE__, __LINE__)

// =============================================================================
// EXEMPLES D'UTILISATION
// =============================================================================

// Exemple 1: Utilisation des macros simples
void example_simple_macros() {
    int* arr;
    double* matrix;
    
    // Au lieu de: arr = malloc(100 * sizeof(int));
    SAFE_MALLOC(arr, 100 * sizeof(int), int);
    
    // Au lieu de: matrix = calloc(50, sizeof(double));
    SAFE_CALLOC(matrix, 50, sizeof(double), double);
    
    // Utilisation normale
    for (int i = 0; i < 100; i++) {
        arr[i] = i;
    }
    
    // Réallocation si nécessaire
    SAFE_REALLOC(arr, 200 * sizeof(int), int);
    
    free(arr);
    free(matrix);
}

// Exemple 2: Utilisation avec fonctions personnalisées
void example_custom_functions() {
    char* buffer;
    int* numbers;
    
    MALLOC_CHECK(buffer, 1024);
    CALLOC_CHECK(numbers, 100, sizeof(int));
    
    strcpy(buffer, "Hello, World!");
    
    for (int i = 0; i < 100; i++) {
        numbers[i] = i * i;
    }
    
    free(buffer);
    free(numbers);
}

// Exemple 3: Utilisation avec nettoyage automatique
void example_with_cleanup() {
    CleanupList* cleanup = init_cleanup_list();
    if (!cleanup) {
        fprintf(stderr, "Impossible d'initialiser la liste de nettoyage\n");
        return;
    }
    
    int* arr1;
    double* arr2;
    char* buffer;
    
    // Ces allocations seront automatiquement libérées en cas d'erreur
    MALLOC_CLEANUP(arr1, 1000 * sizeof(int), cleanup);
    MALLOC_CLEANUP(arr2, 500 * sizeof(double), cleanup);
    MALLOC_CLEANUP(buffer, 256, cleanup);
    
    // Utilisation normale...
    for (int i = 0; i < 1000; i++) {
        arr1[i] = i;
    }
    
    // Si tout va bien, nettoyer manuellement
    cleanup_all(cleanup);
    // Les pointeurs sont automatiquement mis à NULL
}

// =============================================================================
// ADAPTATION POUR VOS STRUCTURES SPÉCIFIQUES
// =============================================================================

// Macro spécialisée pour allouer une Matrix
#define ALLOC_MATRIX(matrix_ptr, dim) do { \
    SAFE_MALLOC(matrix_ptr, sizeof(Matrix), Matrix); \
    matrix_ptr->dim = dim; \
    SAFE_MALLOC(matrix_ptr->matrix, dim * sizeof(double*), double*); \
    for (int i = 0; i < dim; i++) { \
        SAFE_MALLOC(matrix_ptr->matrix[i], dim * sizeof(double), double); \
    } \
} while(0)

// Macro spécialisée pour allouer un State
#define ALLOC_STATE(state_ptr, size) do { \
    SAFE_MALLOC(state_ptr, sizeof(State), State); \
    state_ptr->size = size; \
    SAFE_CALLOC(state_ptr->occupancy, size, sizeof(long long), long long); \
} while(0)

// Macro spécialisée pour allouer une StateList
#define ALLOC_STATE_LIST(list_ptr, count) do { \
    SAFE_MALLOC(list_ptr, sizeof(StateList), StateList); \
    list_ptr->count = count; \
    SAFE_MALLOC(list_ptr->states, count * sizeof(State), State); \
} while(0)

#endif // MALLOC_CHECK_H