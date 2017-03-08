#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <argp.h>
#include "openmx-hks-lib.h"
#include "simplemat.h"
#include "simplejson.h"

#define ACTION_SET_FERMI 0
#define ACTION_COPY_FERMI 1
#define ACTION_DISPLAY 2
#define ACTION_EXTRACT_HAMILTONIAN 3
#define ACTION_SHIFT_HAMILTONIAN 4

#define BohrR 0.529177249

const char *argp_program_version = "openmx-hks " VERSION;
const char *argp_program_bug_address = "<gpulkin@gmail.com>";
static char doc[] = "openmx-hks: performs various operations on OpenMX HKS files";
static char args_doc[] = "ACTION (display, copy-fermi, set-fermi, shift-hamiltonian, extract-hamiltonian) FILE [ARGS]";
static struct argp_option options[] = {
    {"verbose", 'v', 0, 0, "Verbose output" },
    {"float", 'f', "FORMAT", 0, "Float format" },
    { 0 }
};

struct arguments {
    int action;
    int verbose;
    char *float_format;
    
    char *input;
    char *output;
    
    double fermi;
    double shift;
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
            if (state->arg_num >= 3) argp_usage(state);
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
                    case ACTION_DISPLAY:
                        argp_usage(state);
                        break;
                    case ACTION_EXTRACT_HAMILTONIAN:
                        arguments->output = arg;
                        break;
                    case ACTION_SHIFT_HAMILTONIAN:
                        if (arg[0] == '_') arg[0] = '-';
                        arguments->shift = atof(arg);
                        break;
                }
            }
            break;
            
        case ARGP_KEY_END:
            if (state->arg_num < 2) argp_usage(state);
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

void print_hks(struct hks_data *data, int verbosity) {
    
    int i,j;
    
    printf("Version %i %i\n", data->version_major, data->version_minor);
    printf("Unstructured data: %lu bytes\n", data->rest_length);
    
    printf("Spin Hamiltonian: ");
    switch (data->spin_mode) {
        case SPIN_MODE_NONE:
            printf("none");
            break;
        case SPIN_MODE_FULL:
            printf("fully-realtivistic");
            break;
        default:
            printf("unknown (%d)",data->spin_mode);
    }
    printf(" \tint: %d hex: ", data->spin_mode);
    print_bytes(&data->spin_mode, sizeof(data->spin_mode));
    printf("\n");
    printf("Fermi level: \t%2.10f hex: ", data->fermi);
    print_bytes(&data->fermi, sizeof(data->fermi));
    printf("\n");
    printf("Neighbours: %d\n", data->cell_replica_number);
    printf("\nUnit cell:\n");
    if (verbosity == 0) {
        for (i=0; i<3; i++) {
            for (j=0; j<3; j++)
                printf("\t %16.13f", data->unit_cell_vectors[i][j]);
            printf("\n");
        }
    } else if (verbosity>0) {
        printf("  a.u:\n    ");
        for (i=0; i<3; i++) {
            for (j=0; j<3; j++) 
                printf("\t %16.13f", data->unit_cell_vectors[i][j]);
            printf("\n");
        }
        printf("  angstrom:\n    ");
        for (i=0; i<3; i++) {
            for (j=0; j<3; j++) 
                printf("\t %16.13f", data->unit_cell_vectors[i][j]*BohrR);
            printf("\n");
        }
    }
    if (verbosity>0) {
        printf("\nSpecies: %d\n",data->species_number);
        printf("\tID\tbasis\tCbasis\n");
        for (i=0; i<data->species_number; i++) {
            printf("\t%d\t%d\t%d\n",data->species[i].id,data->species[i].basis_size,data->species[i].contracted_basis_size);
        }
    }
    printf("\nAtoms: %d\n", data->atoms_number);
    for (i=0; i<data->atoms_number; i++) {
        struct atom a = data->atoms[i];
        if (verbosity == 0)
            printf("\t%d\t%16.13f\t%16.13f\t%16.13f\n", a.specimen->id, a.coordinates[0], a.coordinates[1], a.coordinates[2]);
        else if (verbosity>0) {
            printf("  #%d\n",a.id);
            printf("    specimen: \t%d\n", a.specimen->id);
            printf("    coordinates, a.u.    : \t%16.13f\t%16.13f\t%16.13f\n", a.coordinates[0], a.coordinates[1], a.coordinates[2]);
            printf("    coordinates, angstrom: \t%16.13f\t%16.13f\t%16.13f\n", a.coordinates[0]*BohrR, a.coordinates[1]*BohrR, a.coordinates[2]*BohrR);
            printf("    neighbours: %d \t(",a.neighbours_number);
            for (j=0; j<a.neighbours_number; j++) {
                struct atom_replica ar = a.neighbours[j];
                printf("%d[%d,%d,%d]",ar.atom->id,ar.cell->index[0],ar.cell->index[1],ar.cell->index[2]);
                if (j<a.neighbours_number-1) printf(",");
            }
            printf(")\n");
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
        printf("Could not open file %s for reading\n", name);
        exit(1);
    }
    
    struct hks_data data;

    switch (read_hks(f, &data)) {
        case SUCCESS:
            if (verbosity>-1) {
                printf("----------\nFile %s\n----------\n", name);
                print_hks(&data, verbosity);
                printf("\n");
            }
            break;
        case ERR_VERSION:
            printf("Unsupported version %i %i\n", data.version_major, data.version_minor);
            exit(1);
    }
    
    fclose(f);

    return data;

}

void write_and_print_hks(char *name, struct hks_data *data, int verbosity) {

    FILE *f;
    
    if (!(f = fopen(name, "w"))) {
        printf("Could not open file %s for writing\n", name);
        exit(1);
    }

    switch (write_hks(f, data)) {
        case SUCCESS:
            if (verbosity>-1) {
                printf("----------\nOutput to %s\n----------\n", name);
                print_hks(data,verbosity);
                printf("\n");
            }
            break;
        case ERR_VERSION:
            printf("Unsupported version %i %i\n", data->version_major, data->version_minor);
            break;
    }
    
    fclose(f);

}

void write_and_print_blocks(char *name, struct hks_data *data, int verbosity) {
    
    int i;
    int last_dot = -1;
    for (i=0; name[i] != 0; i++) if (name[i] == '.') last_dot = i;
    if (last_dot == -1) {
        printf("Please provide a proper extension (\".mat\", \".json\") for the output file\n");
        exit(1);
    }
    
    void (*write_header)(FILE*);
    void (*write_double_scalar)(FILE*, char*, double*);
    void (*write_int_1D_array)(FILE*, char*, int*, int, int);
    void (*write_int_2D_array)(FILE*, char*, int*, int, int);
    void (*write_complex_3D_array)(FILE*, char*, double*, int, int, int);
    void (*write_footer)(FILE*);

    char *ext = name+last_dot+1;
    if (strcmp(ext,"mat") == 0) {
        if (verbosity>-1) printf("Writing Matlab MAT file ...\n");
        write_header = &write_mat_header;
        write_double_scalar = &write_mat_double_scalar;
        write_int_1D_array = &write_mat_int_1D_array;
        write_int_2D_array = &write_mat_int_2D_array;
        write_complex_3D_array = &write_mat_complex_3D_array;
        write_footer = &write_mat_footer;
    } else if (strcmp(ext,"json") == 0) {
        if (verbosity>-1) printf("Writing JSON file ...\n");
        write_header = &write_json_header;
        write_double_scalar = &write_json_double_scalar;
        write_int_1D_array = &write_json_int_1D_array;
        write_int_2D_array = &write_json_int_2D_array;
        write_complex_3D_array = &write_json_complex_3D_array;
        write_footer = &write_json_footer;
    } else {
        printf("Please provide a proper extension (\".mat\", \".json\") for the output file\n");
        exit(1);
    }

    FILE *f;
    if (!(f = fopen(name, "w"))) {
        printf("Could not open file %s for writing\n", name);
        exit(1);
    }
    
        
    if (verbosity>-1) {
        printf("Compressing HKS ... ");
    }
    compress(data);
    if (verbosity>-1) {
        printf("Done\n");
    }

    if (verbosity>-1) {
        printf("----------\nOutput to %s\n----------\n", name);
        print_hks(data,verbosity);
        printf("\n");
    }
    struct basis_description basis;
    make_basis(data, &basis);
    
    (*write_header)(f);
    (*write_double_scalar)(f, "fermi", &(data->fermi));
    (*write_int_1D_array)(f,"basis_spin",(int*)basis.r2s,basis.size,3);
    (*write_int_1D_array)(f,"basis_atom",(int*)basis.r2s+1,basis.size,3);
    (*write_int_1D_array)(f,"basis_orbital",(int*)basis.r2s+2,basis.size,3);
    
    char non_zero[data->cell_replica_number];
    int nv[data->cell_replica_number*3];
    struct F_complex *H = malloc(basis.size*basis.size*data->cell_replica_number*sizeof(struct F_complex));
    struct F_complex *S = malloc(basis.size*basis.size*data->cell_replica_number*sizeof(struct F_complex));
    if (!H || !S) {
        printf("Could not allocate memory, exiting");
        exit(1);
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
    
    dispose_basis(&basis);
    free(H);
    free(S);
   
    fclose(f);
    
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
    }

    return 0;
}
