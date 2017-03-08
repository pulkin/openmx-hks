void write_json_header(FILE *f);
void write_json_double_scalar(FILE *f, char* name, double* data);
void write_json_int_1D_array(FILE *f, char* name, int* data, int n, int step);
void write_json_int_2D_array(FILE *f, char* name, int* data, int n, int m);
void write_json_complex_3D_array(FILE *f, char* name, double* data, int n, int m, int k);
void write_json_footer(FILE *f);
