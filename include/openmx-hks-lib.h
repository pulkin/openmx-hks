#define VERSION "0.1"
#define DEBUG 1

#define SUCCESS 0
#define ERR_FILE_IO 1
#define ERR_VERSION 2
#define ERR_FILE_STRUCTURE 3

#define SPIN_MODE_NONE 0
#define SPIN_MODE_FULL 3

#define SPINORB_FIRST_DIM 3

#define SPIN_SIZE(a) ((a)->spin_mode+1)
#define SPIN_BASIS_SIZE(a) ((a)->spin_mode == 0 ? 1 : ((a)->spin_mode == 3 ? 2 : -1))
#define TOTAL_GRID_SIZE(a) ((a)->grid_size[0] * (a)->grid_size[1] * (a)->grid_size[2])
#define INDEX_OF(array, element) (int)(element-array)

struct specimen {

    int id;
    int basis_size;
    int contracted_basis_size;

};

struct atom_replica {

    struct atom *atom;
    struct atom *source;
    struct cell_replica *cell;
    double ***overlap;
    double ***hamiltonian;
    double ***spinorb;
    double ***density;

};

struct atom {

    int id;
    struct specimen *specimen;
    double coordinates[3];
    int neighbours_number;
    struct atom_replica *neighbours;

};

struct cell_replica {
    
    int id;
    int index[3];
    
};

struct hks_data {

    // Renamed
    int version_major;
    int version_minor;
    int spin_mode;
    int atoms_number;
    int species_number;
    int first_neighbors_max;
    int grid_size[3];
    double unit_cell_vectors[3][3];
    double grid_vectors[3][3];
    double grid_origin[3];
    double fermi;
    double *charge_density;
    double *hartree_potential;
    
    struct specimen *species;
    struct atom *atoms;
    struct cell_replica *cell_replicas;
    int cell_replica_number;
    
    // Kept
    int Matomnum;
    int MatomnumF;
    int MatomnumS;
    int Num_Cells0;
    double ScaleSize;
    
    // Not read
    unsigned int rest_length;
    char *rest;

};

struct plain_index {
    int spin;
    int atom;
    int orbital;
};

struct basis_description {
    struct hks_data *data;
    int ***s2r;
    struct plain_index *r2s;
    int size;
};

struct F_complex {double r,i;};

int read_hks(FILE *f, struct hks_data *data);
int write_hks(FILE *f, struct hks_data *data);
void make_basis(struct hks_data *data, struct basis_description *basis);
void dispose_basis(struct basis_description *basis);
int calculate_block(struct basis_description *basis, int x, int y, int z, struct F_complex *H, struct F_complex *S);
void dense2csr(struct F_complex *data, int w, int h, int offset, struct F_complex *out_data, int *out_indices, int *out_indptr);
void dispose_hks(struct hks_data *data);
void slice_basis(struct basis_description *basis, int *slice);
int block_number(struct hks_data *data, char* blocks);
void compress(struct hks_data *data);
