#define miMATRIX 14
#define miUINT32 6
#define mxDOUBLE_CLASS 6
#define mxINT32_CLASS 12
#define miINT32 5
#define miINT8 1
#define miDOUBLE 9

void write_mat_header(FILE *f);
int calculate_size(int dims, int name_size, int array_size, int blocks);
void write_mat_array_header(FILE *f, int total_size, int kind);
void write_mat_array_shape(FILE *f, int rank, int *dims);
void write_mat_name(FILE *f, char* name);
void write_mat_plain_double(FILE *f, double *data, int size, int step);
void write_mat_plain_int(FILE *f, int *data, int size, int step);

void write_mat_double_2D_array(FILE *f, char* name, double* data, int n, int m);
void write_mat_complex_3D_array(FILE *f, char* name, double* data, int n, int m, int k);
void write_mat_int_2D_array(FILE *f, char* name, int* data, int n, int m);
