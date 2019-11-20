void write_json_header(void *f);
void write_json_double_scalar(void *f, char* name, double* data);
void write_json_int_1D_array(void *f, char* name, int* data, int n, int step);
void write_json_int_2D_array(void *f, char* name, int* data, int n, int m);
void write_json_complex_3D_array(void *f, char* name, double* data, int n, int m, int k);
void write_json_footer(void *f);
