void write_mat_header(FILE *f);
void write_mat_double_scalar(FILE *f, char* name, double* data);
void write_mat_int_1D_array(FILE *f, char* name, int* data, int n, int step);
void write_mat_int_2D_array(FILE *f, char* name, int* data, int n, int m);
void write_mat_complex_3D_array(FILE *f, char* name, double* data, int n, int m, int k);
void write_mat_footer(FILE *f);

