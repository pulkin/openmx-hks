void *open_h5(char* fname);
void close_h5(void *f);
void write_h5_header(void *f);
void write_h5_double_scalar(void *f, char* name, double* data);
void write_h5_int_1D_array(void *f, char* name, int* data, int n, int step);
void write_h5_int_2D_array(void *f, char* name, int* data, int n, int m);
void write_h5_complex_3D_array(void *f, char* name, double* data, int n, int m, int k);
void write_h5_footer(void *f);
