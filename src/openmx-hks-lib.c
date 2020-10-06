#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <argp.h>
#include "openmx-hks-lib.h"

#define SAFE_READ(destination, size, number, file) do {\
    int x;\
    if ((x=fread(destination, size, number, file))!=number) {\
        if (DEBUG) printf("SAFE_READ failure: %d vs %d\n", x, number);\
        return ERR_FILE_STRUCTURE;\
    }\
} while(0)
#define FREE_3D(array, n, m) do {if (array) {\
    int _i;\
    for (_i=0; _i<n; _i++) if (array[_i]) {\
        int _j;\
        for (_j=0; _j<m; _j++) if (array[_i][_j]) free(array[_i][_j]);\
        free(array[_i]);\
    }\
    free(array);\
}} while(0)

int read_hks(FILE *f, struct hks_data *data) {
/* Reads hks data from file with the name provided as an argument.
 * 
 *     f: file to read from
 *     data: hks data
 * 
 * Returns: status of operation. */
 
    int i,j,k,l;
    long current_position, eof_position;
    
    SAFE_READ(&data->version_major, sizeof(int), 1, f);
    SAFE_READ(&data->version_minor, sizeof(int), 1, f);
    
    if (data->version_major == 1 && data->version_minor == 0) {
        
        SAFE_READ(&data->spin_mode, sizeof(int), 1, f);
        SAFE_READ(&data->atoms_number, sizeof(int), 1, f);
        SAFE_READ(&data->species_number, sizeof(int), 1, f);
        SAFE_READ(&data->first_neighbors_max, sizeof(int), 1, f);
        SAFE_READ(&data->cell_replica_number, sizeof(int), 1, f);
        data->cell_replica_number++;
        SAFE_READ(&data->Matomnum, sizeof(int), 1, f);
        SAFE_READ(&data->MatomnumF, sizeof(int), 1, f);
        SAFE_READ(&data->MatomnumS, sizeof(int), 1, f);
        SAFE_READ(data->grid_size, sizeof(int), 3, f);
        SAFE_READ(&data->Num_Cells0, sizeof(int), 1, f);
        SAFE_READ(&data->ScaleSize, sizeof(double), 1, f);
        SAFE_READ(data->unit_cell_vectors, sizeof(double), 9, f);
        SAFE_READ(data->grid_vectors, sizeof(double), 9, f);
        SAFE_READ(data->grid_origin, sizeof(double), 3, f);

        // Atomic coordinates
        data->atoms = (struct atom*)malloc(sizeof(struct atom)*(data->atoms_number));
        for (i=0; i<data->atoms_number; i++) {
            SAFE_READ(data->atoms[i].coordinates, sizeof(double), 3, f);
            data->atoms[i].id = i+1;
        }
        
        // Fermi
        SAFE_READ(&data->fermi, sizeof(double), 1, f);
        
        // Skip zeroth specimen
        fseek(f, sizeof(int), SEEK_CUR);
        
        // Read specimens
        // Initialize specimen array first
        data->species = (struct specimen*)malloc(sizeof(struct specimen)*(data->species_number));
        // Then read the data
        for (i=0; i<data->atoms_number; i++) {
            SAFE_READ(&k, sizeof(int), 1, f);
            data->atoms[i].specimen = &data->species[k];
        }
        
        // Specimen specifications
        for (i=0; i<data->species_number; i++) {
            SAFE_READ(&data->species[i].contracted_basis_size, sizeof(int), 1, f);
            data->species[i].id = i;
        }
        for (i=0; i<data->species_number; i++)
            SAFE_READ(&data->species[i].basis_size, sizeof(int), 1, f);

        // Structure: nearest neighbours numbers (including self)
        fseek(f, sizeof(int), SEEK_CUR);
        for (i=0; i<data->atoms_number; i++) {
            SAFE_READ(&data->atoms[i].neighbours_number, sizeof(int), 1, f);
            data->atoms[i].neighbours_number++;
        }

        // Structure: nearest neighbours id
        int size1 = (int)data->first_neighbors_max*data->ScaleSize+1;
        fseek(f, sizeof(int)*size1, SEEK_CUR);
        for (i=0; i< data->atoms_number; i++) {
            data->atoms[i].neighbours = (struct atom_replica*)malloc(sizeof(struct atom_replica)*data->atoms[i].neighbours_number);
            for (j=0; j<data->atoms[i].neighbours_number; j++) {
                SAFE_READ(&k,sizeof(int),1,f);
                data->atoms[i].neighbours[j].atom = &data->atoms[k-1];
                data->atoms[i].neighbours[j].source = &data->atoms[i];
            }
            fseek(f, sizeof(int)*(size1-data->atoms[i].neighbours_number), SEEK_CUR);
        }

        // Structure: nearest neighbours cell id
        // Initialize cell replica array first
        data->cell_replicas = (struct cell_replica*)malloc(sizeof(struct cell_replica)*(data->cell_replica_number));
        fseek(f, sizeof(int)*size1, SEEK_CUR);
        for (i=0; i< data->atoms_number; i++) {
            for (j=0; j<data->atoms[i].neighbours_number; j++) {
                SAFE_READ(&k,sizeof(int),1,f);
                data->atoms[i].neighbours[j].cell = &data->cell_replicas[k];
            }
            fseek(f, sizeof(int)*(size1-data->atoms[i].neighbours_number), SEEK_CUR);
        }

        // Structure: neighbour cell
        for (i=0; i< data->cell_replica_number; i++) {
            fseek(f, sizeof(int), SEEK_CUR);
            SAFE_READ(data->cell_replicas[i].index,sizeof(int),3,f);
            data->cell_replicas[i].id = i;
        }

        // Overlap
        for (i=0; i<4; i++) {
            for (j=0; j<data->atoms_number; j++) {
                int basis_1 = data->atoms[j].specimen->basis_size;
                for (k=0; k<data->atoms[j].neighbours_number; k++) {
                    struct atom_replica *rpl = &data->atoms[j].neighbours[k];
                    int basis_2 = rpl->atom->specimen->basis_size;
                    if (i==0)
                        rpl->overlap = (double***)malloc(sizeof(double**)*4);
                    rpl->overlap[i] = (double**)malloc(sizeof(double*)*basis_1);
                    for (l=0; l<basis_1; l++) {
                        rpl->overlap[i][l] = (double*)malloc(sizeof(double)*basis_2);
                        SAFE_READ(rpl->overlap[i][l], sizeof(double), basis_2, f);
                    }
                }
            }
        }

        // Hamiltonian
        for (i=0; i<SPIN_SIZE(data); i++) {
            for (j=0; j<data->atoms_number; j++) {
                int basis_1 = data->atoms[j].specimen->basis_size;
                for (k=0; k<data->atoms[j].neighbours_number; k++) {
                    struct atom_replica *rpl = &data->atoms[j].neighbours[k];
                    int basis_2 = rpl->atom->specimen->basis_size;
                    if (i==0)
                        rpl->hamiltonian = (double***)malloc(sizeof(double**)*SPIN_SIZE(data));
                    rpl->hamiltonian[i] = (double**)malloc(sizeof(double*)*basis_1);
                    for (l=0; l<basis_1; l++) {
                        rpl->hamiltonian[i][l] = (double*)malloc(sizeof(double)*basis_2);
                        SAFE_READ(rpl->hamiltonian[i][l], sizeof(double), basis_2, f);
                    }
                }
            }
        }
        
        // Hamiltonian spinorb
        if (data->spin_mode == SPIN_MODE_FULL) {
            for (i=0; i<SPINORB_FIRST_DIM; i++) {
                for (j=0; j<data->atoms_number; j++) {
                    int basis_1 = data->atoms[j].specimen->basis_size;
                    for (k=0; k<data->atoms[j].neighbours_number; k++) {
                        struct atom_replica *rpl = &data->atoms[j].neighbours[k];
                        int basis_2 = rpl->atom->specimen->basis_size;
                        if (i==0)
                            rpl->spinorb = (double***)malloc(sizeof(double**)*SPINORB_FIRST_DIM);
                        rpl->spinorb[i] = (double**)malloc(sizeof(double*)*basis_1);
                        for (l=0; l<basis_1; l++) {
                            rpl->spinorb[i][l] = (double*)malloc(sizeof(double)*basis_2);
                            SAFE_READ(rpl->spinorb[i][l], sizeof(double), basis_2, f);
                        }
                    }
                }
            }
        }

        // Density matrix
        for (i=0; i<SPIN_SIZE(data); i++) {
            for (j=0; j<data->atoms_number; j++) {
                int basis_1 = data->atoms[j].specimen->basis_size;
                for (k=0; k<data->atoms[j].neighbours_number; k++) {
                    struct atom_replica *rpl = &data->atoms[j].neighbours[k];
                    int basis_2 = rpl->atom->specimen->basis_size;
                    if (i==0)
                        rpl->density = (double***)malloc(sizeof(double**)*SPIN_SIZE(data));
                    rpl->density[i] = (double**)malloc(sizeof(double*)*basis_1);
                    for (l=0; l<basis_1; l++) {
                        rpl->density[i][l] = (double*)malloc(sizeof(double)*basis_2);
                        SAFE_READ(rpl->density[i][l], sizeof(double), basis_2, f);
                    }
                }
            }
        }

        // Density in real space
        data->charge_density = (double*)malloc(sizeof(double)*TOTAL_GRID_SIZE(data));
        SAFE_READ(data->charge_density, sizeof(double), TOTAL_GRID_SIZE(data), f);

        // V Hartree in real space
        data->hartree_potential = (double*)malloc(sizeof(double)*TOTAL_GRID_SIZE(data));
        SAFE_READ(data->hartree_potential, sizeof(double), TOTAL_GRID_SIZE(data), f);

        // Rest of the file
        current_position = ftell(f);
        fseek(f, 0L, SEEK_END);
        eof_position = ftell(f);
        data->rest_length = eof_position-current_position;
        data->rest = (char*)malloc(sizeof(char)*(data->rest_length));
        fseek(f, current_position, SEEK_SET);
        SAFE_READ(data->rest, sizeof(char), data->rest_length, f);
        
    } else {
        return ERR_VERSION;
    }
    
    return SUCCESS;
}

void write_zeros(FILE *f, int n) {
/* Writes zeros into a file.
 * 
 *     f: a file
 *     n: number of zeros */
    char j[4096] = {0};
    while (n>4096) {
        fwrite(j, sizeof(char), 4096, f);
        n -= 4096;
    }
    fwrite(j, sizeof(char), n, f);
}

int write_hks(FILE *f, struct hks_data *data) {
/* Writes hks data to a file with the name provided as an argument.
 * 
 *     f: file to write to
 *     data: hks data
 * 
 * Returns: status of operation. */
    int i,j,k,l;
    
    if (data->version_major == 1 && data->version_minor == 0) {
        
        fwrite(&data->version_major, sizeof(int), 1, f);
        fwrite(&data->version_minor, sizeof(int), 1, f);
        fwrite(&data->spin_mode, sizeof(int), 1, f);
        fwrite(&data->atoms_number, sizeof(int), 1, f);
        fwrite(&data->species_number, sizeof(int), 1, f);
        fwrite(&data->first_neighbors_max, sizeof(int), 1, f);
        data->cell_replica_number --;
        fwrite(&data->cell_replica_number, sizeof(int), 1, f);
        data->cell_replica_number ++;
        fwrite(&data->Matomnum, sizeof(int), 1, f);
        fwrite(&data->MatomnumF, sizeof(int), 1, f);
        fwrite(&data->MatomnumS, sizeof(int), 1, f);
        fwrite(data->grid_size, sizeof(int), 3, f);
        fwrite(&data->Num_Cells0, sizeof(int), 1, f);
        fwrite(&data->ScaleSize, sizeof(double), 1, f);
        fwrite(data->unit_cell_vectors, sizeof(double), 9, f);
        fwrite(data->grid_vectors, sizeof(double), 9, f);
        fwrite(data->grid_origin, sizeof(double), 3, f);
            
        // Atomic coordinates
        for (i=0; i<data->atoms_number; i++) {
            fwrite(data->atoms[i].coordinates, sizeof(double), 3, f);
        }
        
        // Fermi
        fwrite(&data->fermi, sizeof(double), 1, f);
        
        // Skip zeroth specimen
        write_zeros(f, sizeof(int));
        
        // Write specimens
        for (i=0; i<data->atoms_number; i++) {
            fwrite(&data->atoms[i].specimen->id, sizeof(int), 1, f);
        }

        // Specimen specifications
        for (i=0; i<data->species_number; i++) {
            fwrite(&data->species[i].contracted_basis_size, sizeof(int), 1, f);
        }
        for (i=0; i<data->species_number; i++)
            fwrite(&data->species[i].basis_size, sizeof(int), 1, f);

        // Structure: nearest neighbours numbers (including self)
        write_zeros(f, sizeof(int));
        for (i=0; i<data->atoms_number; i++) {
            data->atoms[i].neighbours_number--;
            fwrite(&data->atoms[i].neighbours_number, sizeof(int), 1, f);
            data->atoms[i].neighbours_number++;
        }

        // Structure: nearest neighbours id
        int size1 = (int)data->first_neighbors_max*data->ScaleSize+1;
        write_zeros(f, sizeof(int)*size1);
        for (i=0; i< data->atoms_number; i++) {
            for (j=0; j<data->atoms[i].neighbours_number; j++) {
                fwrite(&data->atoms[i].neighbours[j].atom->id,sizeof(int),1,f);
            }
            write_zeros(f, sizeof(int)*(size1 - data->atoms[i].neighbours_number));
        }

        // Structure: nearest neighbours cell id
        write_zeros(f, sizeof(int)*size1);
        for (i=0; i< data->atoms_number; i++) {
            for (j=0; j<data->atoms[i].neighbours_number; j++) {
                fwrite(&data->atoms[i].neighbours[j].cell->id,sizeof(int),1,f);
            }
            write_zeros(f, sizeof(int)*(size1 - data->atoms[i].neighbours_number));
        }

        // Structure: neighbour cell
        for (i=0; i< data->cell_replica_number; i++) {
            write_zeros(f, sizeof(int));
            fwrite(data->cell_replicas[i].index,sizeof(int),3,f);
        }

        // Overlap
        for (i=0; i<4; i++) {
            for (j=0; j<data->atoms_number; j++) {
                int basis_1 = data->atoms[j].specimen->basis_size;
                for (k=0; k<data->atoms[j].neighbours_number; k++) {
                    struct atom_replica *rpl = &data->atoms[j].neighbours[k];
                    int basis_2 = rpl->atom->specimen->basis_size;
                    for (l=0; l<basis_1; l++) {
                        fwrite(rpl->overlap[i][l], sizeof(double), basis_2, f);
                    }
                }
            }
        }
        
        // Hamiltonian
        for (i=0; i<SPIN_SIZE(data); i++) {
            for (j=0; j<data->atoms_number; j++) {
                int basis_1 = data->atoms[j].specimen->basis_size;
                for (k=0; k<data->atoms[j].neighbours_number; k++) {
                    struct atom_replica *rpl = &data->atoms[j].neighbours[k];
                    int basis_2 = rpl->atom->specimen->basis_size;
                    for (l=0; l<basis_1; l++) {
                        fwrite(rpl->hamiltonian[i][l], sizeof(double), basis_2, f);
                    }
                }
            }
        }
        
        // Hamiltonian spinorb
        if (data->spin_mode == SPIN_MODE_FULL) {
            for (i=0; i<SPINORB_FIRST_DIM; i++) {
                for (j=0; j<data->atoms_number; j++) {
                    int basis_1 = data->atoms[j].specimen->basis_size;
                    for (k=0; k<data->atoms[j].neighbours_number; k++) {
                        struct atom_replica *rpl = &data->atoms[j].neighbours[k];
                        int basis_2 = rpl->atom->specimen->basis_size;
                        for (l=0; l<basis_1; l++) {
                            fwrite(rpl->spinorb[i][l], sizeof(double), basis_2, f);
                        }
                    }
                }
            }
        }
        
        // Density matrix
        for (i=0; i<SPIN_SIZE(data); i++) {
            for (j=0; j<data->atoms_number; j++) {
                int basis_1 = data->atoms[j].specimen->basis_size;
                for (k=0; k<data->atoms[j].neighbours_number; k++) {
                    struct atom_replica *rpl = &data->atoms[j].neighbours[k];
                    int basis_2 = rpl->atom->specimen->basis_size;
                    for (l=0; l<basis_1; l++) {
                        fwrite(rpl->density[i][l], sizeof(double), basis_2, f);
                    }
                }
            }
        }
        
        // Density in real space
        fwrite(data->charge_density, sizeof(double), TOTAL_GRID_SIZE(data), f);
        // V Hartree in real space
        fwrite(data->hartree_potential, sizeof(double), TOTAL_GRID_SIZE(data), f);
        
        // Rest of the file
        fwrite(data->rest, sizeof(char), data->rest_length, f);
        
    } else {
        return ERR_VERSION;
    }
    
    return SUCCESS;
}

void make_basis(struct hks_data *data, struct basis_description *basis) {
/* Writes a lookup table for a sparse index into plain index.
 * 
 *     data: hks data
 * 
 *     *basis: a pointer to the memory where to write the basis struct */
    
    int i,j,k;
    
    basis->data = data;
    basis->s2r = (int***)malloc(sizeof(int**)*SPIN_BASIS_SIZE(data));
    basis->size = 0;
    
    // Forward "sparse" to "real"
    for (i=0; i<SPIN_BASIS_SIZE(data); i++) {
        basis->s2r[i] = (int**)malloc(sizeof(int*)*data->atoms_number);
        for (j=0; j<data->atoms_number; j++) {
            basis->s2r[i][j] = (int*)malloc(sizeof(int)*data->atoms[j].specimen->basis_size);
            for (k=0; k<data->atoms[j].specimen->basis_size; k++) {
                basis->s2r[i][j][k] = basis->size;
                basis->size ++;
            }
        }
    }

    basis->r2s = (struct plain_index*)malloc(sizeof(struct plain_index)*basis->size);
    int basis_size = 0;
    
    // Backward "real" to "sparse"
    for (i=0; i<SPIN_BASIS_SIZE(data); i++) {
        for (j=0; j<data->atoms_number; j++) {
            for (k=0; k<data->atoms[j].specimen->basis_size; k++) {
                basis->r2s[basis_size].spin = i;
                basis->r2s[basis_size].atom = j;
                basis->r2s[basis_size].orbital = k;
                basis_size++;
            }
        }
    }
}

void slice_basis(struct basis_description *basis, int *slice) {
/* Reduces the basis set by slicing. */

    int i, spin, at, orb;
    int counter = 0;
    
    for (i=0; i<basis->size; i++) {
        spin = basis->r2s[i].spin;
        at = basis->r2s[i].atom;
        orb = basis->r2s[i].orbital;
        
        if (slice[i]) {
            basis->s2r[spin][at][orb] = counter;
            basis->r2s[counter].spin = spin;
            basis->r2s[counter].atom = at;
            basis->r2s[counter].orbital = orb;
            counter++;
        } else {
            basis->s2r[spin][at][orb] = -1;
        }
    }
    basis->size = counter;
}

void dispose_basis(struct basis_description *basis) {
/* Disposes sparse-to-real basis */

    int i,j;

    struct hks_data *data = basis->data;
    for (i=0; i<SPIN_BASIS_SIZE(data); i++) {
        for (j=0; j<data->atoms_number; j++) {
            free(basis->s2r[i][j]);
        }
        free(basis->s2r[i]);
    }
    free(basis->s2r);
    free(basis->r2s);
    
}    

int calculate_block(struct basis_description *basis, int x, int y, int z, struct F_complex *H, struct F_complex *S) {
/* Calculates a tight binding block.
 * 
 *     basis: hks basis
 *     x,y,z: block indexes
 *     H,S: pointers to write Hamiltonian and overlap
 * 
 * Returns: 1 if block is non-zero, 0 otherwise. */
 
    struct hks_data *data = basis->data;
    
    int j,j2,k,k2,row,column,index,result = 0;
    
    if (H) memset(H, 0, sizeof(struct F_complex)*basis->size*basis->size);
    if (S) memset(S, 0, sizeof(struct F_complex)*basis->size*basis->size);
    
    // Iterate over elements
    for (j=0; j<data->atoms_number; j++) {
        struct atom *a = data->atoms + j;
        for (j2=0; j2<a->neighbours_number; j2++) {
            
            struct atom_replica r = a->neighbours[j2];
            struct atom *a2 = r.atom;
            
            if (r.cell->index[0] == x && r.cell->index[1] == y && r.cell->index[2] == z) {
                
                result = 1;
                
                if (!H || !S) return result;
                
                for (k=0; k<a->specimen->basis_size; k++) {
                    for (k2=0; k2<a2->specimen->basis_size; k2++) {
                        if (data->spin_mode == SPIN_MODE_FULL) {
                            
                            // Spin 0 0
                            row = basis->s2r[0][j][k];
                            column = basis->s2r[0][INDEX_OF(data->atoms,a2)][k2];
                            if (row>=0 && column>=0) {
                                index = row*basis->size + column;
                                H[index].r += r.hamiltonian[0][k][k2];
                                H[index].i += r.spinorb[0][k][k2];
                                S[index].r += r.overlap[0][k][k2];
                            }
                            
                            // Spin 1 1
                            row = basis->s2r[1][j][k];
                            column = basis->s2r[1][INDEX_OF(data->atoms,a2)][k2];
                            if (row>=0 && column>=0) {
                                index = row*basis->size + column;
                                H[index].r += r.hamiltonian[1][k][k2];
                                H[index].i += r.spinorb[1][k][k2];
                                S[index].r += r.overlap[0][k][k2];
                            }
                            
                            // Spin 0 1
                            row = basis->s2r[0][j][k];
                            column = basis->s2r[1][INDEX_OF(data->atoms,a2)][k2];
                            if (row>=0 && column>=0) {
                                index = row*basis->size + column;
                                H[index].r += r.hamiltonian[2][k][k2];
                                H[index].i += r.hamiltonian[3][k][k2] + r.spinorb[2][k][k2];
                            }
                            
                        } else if (data->spin_mode == SPIN_MODE_NONE) {
                            
                            row = basis->s2r[0][j][k];
                            column = basis->s2r[0][INDEX_OF(data->atoms,a2)][k2];
                            if (row>=0 && column>=0) {
                                index = row*basis->size + column;
                                H[index].r += r.hamiltonian[0][k][k2];
                                S[index].r += r.overlap[0][k][k2];
                            }
                            
                        }
                    }
                }
                
            }
            
            if (r.cell->index[0] == -x && r.cell->index[1] == -y && r.cell->index[2] == -z) {

                result = 1;
                
                if (!H || !S) return result;
                
                for (k=0; k<a->specimen->basis_size; k++) {
                    for (k2=0; k2<a2->specimen->basis_size; k2++) {
                        if (data->spin_mode == SPIN_MODE_FULL) {
    
                            // Spin 1 0 (which is a spin 0 1 in a conjugate block of a matrix)
                            row = basis->s2r[0][j][k];
                            column = basis->s2r[1][INDEX_OF(data->atoms,a2)][k2];
                            if (row>=0 && column>=0) {
                                index = column*basis->size + row;
                                H[index].r += r.hamiltonian[2][k][k2];
                                H[index].i -= r.hamiltonian[3][k][k2] + r.spinorb[2][k][k2];
                            }
                            
                        }
                    }
                }
            }
        }
    }
    return result;
}

void dense2csr(struct F_complex *data, int w, int h, int offset, struct F_complex *out_data, int *out_indices, int *out_indptr) {
/* Converts dense complex dense matrix into compressed sparse row representation.
 *
 *      data: dense 2D matrix;
 *      w: matrix width;
 *      h: matrix height;
 *      offset: offset of *data in a larger array, if any;
 *      out_data: non-zero entries;
 *      out_indices: non-zero element second index;
 *      out_indptr: boundaries of non-zero blocks; */

    int i, j;
    struct F_complex *data_column = data;

    for (i=0; i<w; i++) {
        out_indptr[i] = offset;
        for (j=0; j<h; j++) if (data_column[j].r || data_column[j].i) {
            out_data[offset].r = data_column[j].r;
            out_data[offset].i = data_column[j].i;
            out_indices[offset] = j;
            offset++;
        }
        data_column += h;
    }
    out_indptr[w] = offset;
}

int block_number(struct hks_data *data, char* blocks) {
/* Calculates the number of non-zero blocks. Optionally, stores 1 in
 * the array supplied if the block is non-zero otherwise stores 0.
 * 
 *      data: hks data;
 * 
 *      blocks: array of chars of the size data->cell_replica_number.
 *      NULL is accepted;
 * 
 * Returns: the number of non-zero blocks. */
    
    char do_free = 0;
    
    if (!blocks) {
        blocks = malloc(data->cell_replica_number);
        do_free = 1;
    }
    memset(blocks, 0, data->cell_replica_number);
    
    int j, j2;
    
    for (j=0; j<data->atoms_number; j++) {
        struct atom *a = data->atoms + j;
        for (j2=0; j2<a->neighbours_number; j2++) blocks[a->neighbours[j2].cell->id] = 1;
    }
    
    int result = 0;
    for (j=0; j<data->cell_replica_number; j++) result += blocks[j];
    
    if (do_free) free(blocks);
    return result;
}

void compress(struct hks_data *data) {
/* Compresses HKS data and removes zero blocks.
 * 
 *      data: hks data; */
 
    char nonzero[data->cell_replica_number];
    int N = block_number(data, nonzero);
    int id[data->cell_replica_number];
    id[0] = 0;
    
    int j,j2;
    for (j=1; j<data->cell_replica_number; j++) id[j] = id[j-1] + nonzero[j-1];
    
    // Fix pointer data->atoms.neighbours->cell_replica
    for (j=0; j<data->atoms_number; j++) {
        struct atom *a = data->atoms + j;
        for (j2=0; j2<a->neighbours_number; j2++) {
            int old_id = a->neighbours[j2].cell->id;
            a->neighbours[j2].cell = data->cell_replicas + id[old_id];
        }
    }
    
    // Fix array data->cell_replicas
    j2 = 0;
    for (j=0; j<data->cell_replica_number; j++) if (nonzero[j]) {
        data->cell_replicas[j2] = data->cell_replicas[j];
        j2++;
    }
    
    // Fix number of replicas
    data->cell_replica_number = N;
}

void dispose_hks(struct hks_data *data) {
/* Disposes hks.
 * 
 *     data - hks structure */
    int i;
    
    if (data->atoms) { // data->atoms
        for (i=0; i<data->atoms_number; i++) {
            struct atom a = data->atoms[i];
            if (a.neighbours) { // data->atoms[i].neighbours
                int j;
                for (j=0; j<a.neighbours_number; j++) {
                    struct atom_replica r = a.neighbours[j];
                    FREE_3D(r.overlap,4,a.specimen->basis_size);
                    FREE_3D(r.hamiltonian,SPIN_SIZE(data),a.specimen->basis_size);
                    if (data->spin_mode == SPIN_MODE_FULL) {
                        FREE_3D(r.spinorb,SPINORB_FIRST_DIM,a.specimen->basis_size);
                    }
                    FREE_3D(r.density,SPIN_SIZE(data),a.specimen->basis_size);
                }
                free(a.neighbours);
            }
        }
        free(data->atoms);
    }
    
    if (data->charge_density) free(data->charge_density);
    if (data->hartree_potential) free(data->hartree_potential);
    if (data->cell_replicas) free(data->cell_replicas);
    if (data->species) free(data->species);
    if (data->rest) free(data->rest);
}
