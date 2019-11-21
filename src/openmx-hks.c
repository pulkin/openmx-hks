#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <argp.h>
#include "openmx-hks-lib.h"
#include "simplemat.h"
#include "simplejson.h"
#include "simpleh5.h"

#define ACTION_SET_FERMI 0
#define ACTION_COPY_FERMI 1
#define ACTION_DISPLAY 2
#define ACTION_EXTRACT_HAMILTONIAN 3
#define ACTION_SHIFT_HAMILTONIAN 4
#define ACTION_EXTRACT_STRUCTURE 5

#define BohrR 0.529177249
#define Hartree 27.2113845

const char *argp_program_version = "openmx-hks " VERSION;
const char *argp_program_bug_address = "<gpulkin@gmail.com>";
static char doc[] = "openmx-hks: performs various operations on OpenMX HKS files";
static char args_doc[] = "ACTION (display, copy-fermi, set-fermi, shift-hamiltonian, extract-hamiltonian, extract-structure) FILE [ARGS]";
static struct argp_option options[] = {
    {"verbose", 'v', 0, 0, "Verbose output" },
    {"float", 'f', "FORMAT", 0, "Float format" },
    { 0 }

};
static char size_suffixes[5] = "bKMGT";

struct arguments {
    int action;
    int verbose;
    char *float_format;
    
    char *input;
    char *output;
    
    double fermi;
    double shift;
    char *atom_names;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = state->input;
    
    switch (key) {
        case 'v':
            arguments->verbose += 1;
            break;
            
        case 'f':
            arguments->float_format = arg;
            break;
            
        case ARGP_KEY_ARG:
            if (state->arg_num >= 4) argp_usage(state);
            if (state->arg_num == 0) {
                if (strcmp(arg, "set-fermi") == 0)
                    arguments->action = ACTION_SET_FERMI;
                else if (strcmp(arg, "copy-fermi") == 0)
                    arguments->action = ACTION_COPY_FERMI;
                else if (strcmp(arg, "display") == 0)
                    arguments->action = ACTION_DISPLAY;
                else if (strcmp(arg, "extract-hamiltonian") == 0)
                    arguments->action = ACTION_EXTRACT_HAMILTONIAN;
                else if (strcmp(arg, "shift-hamiltonian") == 0)
                    arguments->action = ACTION_SHIFT_HAMILTONIAN;
                else if (strcmp(arg, "extract-structure") == 0)
                    arguments->action = ACTION_EXTRACT_STRUCTURE;
                else argp_usage(state);
            }
            if (state->arg_num == 1) {
                switch (arguments->action) {
                    case ACTION_SET_FERMI:
                        arguments->output = arg;
                        break;
                    case ACTION_COPY_FERMI:
                        arguments->input = arg;
                        break;
                    case ACTION_DISPLAY:
                        arguments->input = arg;
                        break;
                    case ACTION_EXTRACT_HAMILTONIAN:
                        arguments->input = arg;
                        break;
                    case ACTION_SHIFT_HAMILTONIAN:
                        arguments->output = arg;
                        break;
                    case ACTION_EXTRACT_STRUCTURE:
                        arguments->input = arg;
                        break;
                    default:
                        argp_usage(state);
                        break;
                }
            }
            if (state->arg_num == 2) {
                switch (arguments->action) {
                    case ACTION_SET_FERMI:
                        if (arg[0] == '_') arg[0] = '-';
                        arguments->fermi = atof(arg);
                        break;
                    case ACTION_COPY_FERMI:
                        arguments->output = arg;
                        break;
                    case ACTION_EXTRACT_HAMILTONIAN:
                        arguments->output = arg;
                        break;
                    case ACTION_SHIFT_HAMILTONIAN:
                        if (arg[0] == '_') arg[0] = '-';
                        arguments->shift = atof(arg);
                        break;
                    case ACTION_EXTRACT_STRUCTURE:
                        arguments->output = arg;
                        break;
                    default:
                        argp_usage(state);
                        break;
                }
            }
            if (state->arg_num == 3) {
                switch (arguments->action) {
                    case ACTION_EXTRACT_STRUCTURE:
                        arguments->atom_names = arg;
                        break;
                    default:
                        argp_usage(state);
                        break;
                }
            }
            break;
            
        case ARGP_KEY_END:
            if (state->arg_num < 2) argp_usage(state);
            switch (arguments->action) {
                case ACTION_SET_FERMI:
                    if (state->arg_num != 3) argp_usage(state);
                    break;
                case ACTION_COPY_FERMI:
                    if (state->arg_num != 3) argp_usage(state);
                    break;
                case ACTION_DISPLAY:
                    if (state->arg_num != 2) argp_usage(state);
                    break;
                case ACTION_EXTRACT_HAMILTONIAN:
                    if (state->arg_num != 3) argp_usage(state);
                    break;
                case ACTION_SHIFT_HAMILTONIAN:
                    if (state->arg_num != 3) argp_usage(state);
                    break;
                case ACTION_EXTRACT_STRUCTURE:
                    if (state->arg_num != 4) argp_usage(state);
                    break;
                default:
                    argp_usage(state);
                    break;
            }
            break;
            
        default:
            return ARGP_ERR_UNKNOWN;
    }
    
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

void print_bytes(const void *object, size_t size)
{
  size_t i;

  for(i = 0; i < size; i++) {
    printf("%02X", ((const unsigned char *) object)[i] & 0xff);
    if (i<size-1) printf(" ");
  }
}

void print_size(double s) {
    int i = 0;
    while (s>512 && i < sizeof(size_suffixes)-1) {
        s = s/1024;
        i++;
    }
    if (s<1) printf("%.2f %c", s, size_suffixes[i]);
    else printf("%.0f %c", s, size_suffixes[i]);
}

void print_hks(struct hks_data *data, int verbosity) {
    
    if (verbosity>-1) {
        
        printf("[INFO] File version\t\t\t%i.%i\n", data->version_major, data->version_minor);
        printf("[INFO] Unrecognized data:\t\t%lu bytes\n", data->rest_length);
        
        printf("[INFO] Spin treatment:\t\t\t");
        switch (data->spin_mode) {
            case SPIN_MODE_NONE:
            printf("none");
            break;
            case SPIN_MODE_FULL:
            printf("fully-realtivistic");
            break;
            default:
            printf("unknown (int: %d, hex:",data->spin_mode);
                    print_bytes(&data->spin_mode, sizeof(data->spin_mode));
                    printf(")");
        }
        printf("\n");
        printf("[INFO] Fermi level (Hartree, data):\t%2.10f\n", data->fermi);
        printf("[INFO] Fermi level (eV, calculated):\t%2.10f\n", data->fermi*Hartree);
        printf("[INFO] Species:\t\t\t\t%d\n", data->species_number);
        printf("[INFO] Atoms per unit cell:\t\t%d\n", data->atoms_number);
        printf("[INFO] TB neighbours:\t\t\t%d\n", data->cell_replica_number);
        
    }
    
    if (verbosity > 0) {
        
        int i,j;
        
        printf("[INFO] Unit cell matrix (Bohr radius units, data):\n");
        for (i=0; i<3; i++) {
            printf("[INFO]");
            for (j=0; j<3; j++)
                printf(" %17.13f", data->unit_cell_vectors[i][j]);
            printf("\n");
        }
        
        printf("[INFO] Unit cell matrix (angstroms, calculated):\n");
        for (i=0; i<3; i++) {
            printf("[INFO]");
            for (j=0; j<3; j++)
                printf(" %17.13f", data->unit_cell_vectors[i][j]*BohrR);
            printf("\n");
        }
        
        printf("[INFO] Spherical basis functions per specimen:");
        for (i=0; i<data->species_number; i++) {
            printf(" #%d: %d",data->species[i].id,data->species[i].basis_size);
            if (i<data->species_number-1) printf(",");
        }
        printf("\n");

        printf("[INFO] Atomic cartesian coordinates in the unit cell (Bohr radius units, data):\n");
        for (i=0; i<data->atoms_number; i++) {
            struct atom a = data->atoms[i];
            printf("[INFO]  #%d",a.specimen->id);
            for (j=0; j<3; j++)
                printf(" %17.13f", a.coordinates[j]);
            printf("\n");
        }
        
        printf("[INFO] Atomic cartesian coordinates in the unit cell (angstroms, calculated):\n");
        for (i=0; i<data->atoms_number; i++) {
            struct atom a = data->atoms[i];
            printf("[INFO]  #%d",a.specimen->id);
            for (j=0; j<3; j++)
                printf(" %17.13f", a.coordinates[j]*BohrR);
            printf("\n");
        }
        
    }
    
}

void shift_hamiltonian(struct hks_data *data, double shift) {
    int i,j,k,l;
    for (i=0; i<data->atoms_number; i++) {
        struct atom a = data->atoms[i];
        for (j=0; j<a.neighbours_number; j++) {
            struct atom_replica ar = a.neighbours[j];
            for (k=0; k<a.specimen->basis_size; k++) {
                for (l=0; l<ar.atom->specimen->basis_size; l++) {
                    if (data->spin_mode == SPIN_MODE_NONE) {
                        ar.hamiltonian[0][k][l] += ar.overlap[0][k][l]*shift;
                    } else if (data->spin_mode == SPIN_MODE_FULL) {
                        ar.hamiltonian[0][k][l] += ar.overlap[0][k][l]*shift;
                        ar.hamiltonian[1][k][l] += ar.overlap[0][k][l]*shift;
                    }
                }
            }
        }
    }
}

struct hks_data read_and_print_hks(char *name, int verbosity) {

    FILE *f;
    
    if (!(f = fopen(name, "r"))) {
        printf("[ERRO] Could not open file '%s' for reading\n", name);
        exit(1);
    }
    
    struct hks_data data;

    switch (read_hks(f, &data)) {
        case SUCCESS:
            if (verbosity>-1) {
                printf("[INFO] Read file '%s'\n", name);
                print_hks(&data, verbosity);
            }
            break;
        case ERR_VERSION:
            printf("[ERRO] Cannot recognize version %i.%i\n", data.version_major, data.version_minor);
            exit(1);
    }
    
    fclose(f);

    return data;

}

void write_and_print_hks(char *name, struct hks_data *data, int verbosity) {

    FILE *f;
    
    if (!(f = fopen(name, "w"))) {
        printf("[ERRO] Could not open file '%s' for writing\n", name);
        exit(1);
    }

    switch (write_hks(f, data)) {
        case SUCCESS:
            if (verbosity>-1) {
                printf("[INFO] Write to HKS file '%s'\n", name);
                print_hks(data,verbosity);
            }
            break;
        case ERR_VERSION:
            printf("[ERRO] Cannot recognize version %i.%i\n", data->version_major, data->version_minor);
            break;
    }
    
    fclose(f);

}

void write_and_print_xsf(char *name, struct hks_data *data, char* atom_names, int verbosity) {

    FILE *f;
    
    if (!(f = fopen(name, "w"))) {
        printf("[ERRO] Could not open file '%s' for writing\n", name);
        exit(1);
    }
    
    if (verbosity>-1) printf("[INFO] Writing vectors ...\n");
    
    fprintf(f, "CRYSTAL\n");
    fprintf(f, "PRIMVEC\n");
    int i, j;
    for (i=0; i<3; i++)
        fprintf(f, "%.14f %.14f %.14f\n", data->unit_cell_vectors[i][0]*BohrR, data->unit_cell_vectors[i][1]*BohrR, data->unit_cell_vectors[i][2]*BohrR);
    fprintf(f, "CONVVEC\n");
    for (i=0; i<3; i++)
        fprintf(f, "%.14f %.14f %.14f\n", data->unit_cell_vectors[i][0]*BohrR, data->unit_cell_vectors[i][1]*BohrR, data->unit_cell_vectors[i][2]*BohrR);

    if (verbosity>-1) printf("[INFO] Writing coordinates ...\n");
    
    fprintf(f, "PRIMCOORD\n%d 1\n", data->atoms_number);
    for (i=0; i<data->atoms_number; i++) {
        int c=0;
        for (j=0; j<data->atoms[i].specimen->id; j++) {
            while (atom_names[c] != ',') {
                c++;
                if (atom_names[c] == 0) {
                    printf("[ERRO] Not enough specimen specified in '%s': needed %d, found %d\n", atom_names, data->atoms[i].specimen->id+1, j+1);
                    exit(1);
                }
            }
            c++;
        }
        while ((atom_names[c] != ',') && (atom_names[c] != 0)) {
            fprintf(f, "%c", atom_names[c]);
            c++;
        }
        fprintf(f, " %.14f %.14f %.14f\n", data->atoms[i].coordinates[0]*BohrR, data->atoms[i].coordinates[1]*BohrR, data->atoms[i].coordinates[2]*BohrR);
    }
    fclose(f);

}

void write_and_print_blocks(char *name, struct hks_data *data, int verbosity) {
    
    int i;
    int last_dot = -1;
    for (i=0; name[i] != 0; i++) if (name[i] == '.') last_dot = i;
    if (last_dot == -1) {
        printf("[ERRO] Please provide a proper extension (\".mat\", \".json\", \".h5\") for the output file\n");
        exit(1);
    }
    
    void* (*open_)(char*);
    void (*close_)(void*);
    void (*write_header)(void*);
    void (*write_double_scalar)(void*, char*, double*);
    void (*write_int_1D_array)(void*, char*, int*, int, int);
    void (*write_int_2D_array)(void*, char*, int*, int, int);
    void (*write_complex_3D_array)(void*, char*, double*, int, int, int);
    void (*write_footer)(void*);

    char *ext = name+last_dot+1;
    if (strcmp(ext,"mat") == 0) {
        if (verbosity>-1) printf("[INFO] Extracting to Matlab MAT file '%s'\n", name);
        open_ = &open_mat;
        close_ = &close_mat;
        write_header = &write_mat_header;
        write_double_scalar = &write_mat_double_scalar;
        write_int_1D_array = &write_mat_int_1D_array;
        write_int_2D_array = &write_mat_int_2D_array;
        write_complex_3D_array = &write_mat_complex_3D_array;
        write_footer = &write_mat_footer;
    } else if (strcmp(ext,"json") == 0) {
        if (verbosity>-1) printf("[INFO] Extracting to JSON file '%s'\n", name);
        open_ = &open_json;
        close_ = &close_json;
        write_header = &write_json_header;
        write_double_scalar = &write_json_double_scalar;
        write_int_1D_array = &write_json_int_1D_array;
        write_int_2D_array = &write_json_int_2D_array;
        write_complex_3D_array = &write_json_complex_3D_array;
        write_footer = &write_json_footer;
    } else if (strcmp(ext,"h5") == 0) {
        if (verbosity>-1) printf("[INFO] Extracting to H5 file '%s'\n", name);
        open_ = &open_h5;
        close_ = &close_h5;
        write_header = &write_h5_header;
        write_double_scalar = &write_h5_double_scalar;
        write_int_1D_array = &write_h5_int_1D_array;
        write_int_2D_array = &write_h5_int_2D_array;
        write_complex_3D_array = &write_h5_complex_3D_array;
        write_footer = &write_h5_footer;
    } else {
        printf("[ERRO] Please provide a proper extension (\".mat\", \".json\") for the output file\n");
        exit(1);
    }

    void *f;
    if (!(f = open_(name))) {
        printf("[ERRO] Could not open file '%s' for storage\n", name);
        exit(1);
    }
    
        
    if (verbosity>-1) {
        printf("[INFO] Compressing HKS ... ");
    }
    compress(data);
    if (verbosity>-1) {
        printf("Done\n");
        printf("[INFO] Resulting TB neighbours:\t%d\n",data->cell_replica_number);
    }
    
    struct basis_description basis;
    make_basis(data, &basis);
    if (verbosity>-1) {
        printf("[INFO] Resulting TB block size:\t%d\n",basis.size);
        printf("[INFO] Resulting TB parameters:\t%ld x2 = %ld\n",basis.size*basis.size*data->cell_replica_number,2*basis.size*basis.size*data->cell_replica_number);
    }
    if (verbosity>0) {
        long int defined = 0;
        int j, k, k2;
        for (i=0; i<data->atoms_number; i++) {
            struct atom *a = data->atoms + i;
            for (j=0; j<a->neighbours_number; j++) {
                struct atom_replica r = a->neighbours[j];
                struct atom *a2 = r.atom;
                defined += a->specimen->basis_size * a2->specimen->basis_size * SPIN_SIZE(data);
            }
        }
        printf("[INFO] Defined TB parameters:\t%ld x2 = %ld\n",defined,2*defined);
    }
    
    (*write_header)(f);
    (*write_double_scalar)(f, "fermi", &(data->fermi));
    (*write_int_1D_array)(f,"basis_spin",(int*)basis.r2s,basis.size,3);
    (*write_int_1D_array)(f,"basis_atom",(int*)basis.r2s+1,basis.size,3);
    (*write_int_1D_array)(f,"basis_orbital",(int*)basis.r2s+2,basis.size,3);
    
    char non_zero[data->cell_replica_number];
    int nv[data->cell_replica_number*3];
    if (verbosity>-1) {
        printf("[INFO] Allocating ");
        print_size(2*basis.size*basis.size*data->cell_replica_number*sizeof(struct F_complex));
        printf(" of memory\n");
    }
    struct F_complex *H = malloc(basis.size*basis.size*data->cell_replica_number*sizeof(struct F_complex));
    struct F_complex *S = malloc(basis.size*basis.size*data->cell_replica_number*sizeof(struct F_complex));
    if (!H || !S) {
        printf("[ERRO] Could not allocate memory\n");
        exit(1);
    }
    
    if (verbosity>-1) {
        printf("[INFO] Writing data ... ");
    }
    
    for (i=0; i<data->cell_replica_number; i++) {
        
        int *ind = data->cell_replicas[i].index;
        memcpy(nv+3*i,ind,sizeof(int)*3);
        
        calculate_block(&basis, ind[0], ind[1], ind[2], H + basis.size*basis.size*i, S + basis.size*basis.size*i);
        
    }
    (*write_int_2D_array)(f,"vectors",nv,data->cell_replica_number,3);
    (*write_complex_3D_array)(f, "H", (double*)H, data->cell_replica_number, basis.size, basis.size);
    (*write_complex_3D_array)(f, "S", (double*)S, data->cell_replica_number, basis.size, basis.size);
    (*write_footer)(f);
    
    if (verbosity>-1) {
        printf("done\n");
    }
    if (verbosity>0) {
        long int j;
        long int nonzero = 0;
        for (j=0; j<basis.size*basis.size*data->cell_replica_number; j++) {
            if (!(H[j].r == 0 && H[j].i == 0)) nonzero += 1;
            if (!(S[j].r == 0 && S[j].i == 0)) nonzero += 1;
        }
        printf("[INFO] Non-zero elements:\t%ld\n", nonzero);
        double sparse = 1.0*nonzero/(2*basis.size*basis.size*data->cell_replica_number);
        printf("[INFO] Sparsity:\t\t%.3f\n", 1.0-sparse);
    }
    
    dispose_basis(&basis);
    free(H);
    free(S);
   
    close_(f);
    
}

int main(int argc, char *argv[]) {
    struct arguments arguments;
    arguments.verbose = -1;
    arguments.float_format = "%.14e";
    
    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    
    switch (arguments.action) {
        case ACTION_SET_FERMI: {
            struct hks_data data = read_and_print_hks(arguments.output, arguments.verbose);
            data.fermi = arguments.fermi;
            write_and_print_hks(arguments.output, &data, arguments.verbose);
            dispose_hks(&data);
        }
        break;
        case ACTION_COPY_FERMI: {
            struct hks_data input = read_and_print_hks(arguments.input, arguments.verbose);
            struct hks_data output = read_and_print_hks(arguments.output, arguments.verbose);
            output.fermi = input.fermi;
            write_and_print_hks(arguments.output, &output, arguments.verbose);
            dispose_hks(&input);
            dispose_hks(&output);
        }
        break;
        case ACTION_DISPLAY: {
            struct hks_data input = read_and_print_hks(arguments.input, 1);
            dispose_hks(&input);
        }
        break;
        case ACTION_EXTRACT_HAMILTONIAN: {
            struct hks_data input = read_and_print_hks(arguments.input, arguments.verbose);
            write_and_print_blocks(arguments.output, &input, arguments.verbose);
            dispose_hks(&input);
        }
        break;
        case ACTION_SHIFT_HAMILTONIAN: {
            struct hks_data output = read_and_print_hks(arguments.output, arguments.verbose);
            shift_hamiltonian(&output, arguments.shift);
            write_and_print_hks(arguments.output, &output, arguments.verbose);
            dispose_hks(&output);
        }
        case ACTION_EXTRACT_STRUCTURE: {
            struct hks_data input = read_and_print_hks(arguments.input, arguments.verbose);
            write_and_print_xsf(arguments.output, &input, arguments.atom_names, arguments.verbose);
            dispose_hks(&input);
        }
        break;
    }

    return 0;
}
