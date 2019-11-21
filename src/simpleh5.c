#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "simpleh5.h"
#include "hdf5.h"

void *open_h5(char* fname) {
    hid_t f = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (f == 0) return NULL;
    hid_t* r = malloc(sizeof(hid_t));
    r[0] = f;
    return r;
}
void close_h5(void *f) {H5Fclose(((hid_t*) f)[0]); free(f);}
void write_h5_header(void *f) {}
void write_h5_footer(void *f) {}
void write_h5_double_scalar(void *f, char* name, double* data) {
    hid_t file = ((hid_t*) f)[0];
    hid_t space = H5Screate(H5S_SCALAR);
    // hid_t attr = H5Acreate2(fid, name, H5T_NATIVE_DOUBLE, aid, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset = H5Dcreate(file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
}
void write_h5_int_1D_array(void *f, char* name, int* data, int n, int stride) {
    hid_t file = ((hid_t*) f)[0];
    hsize_t size[] = {n};
    hid_t space = H5Screate_simple(1, size, NULL);
    hid_t dset = H5Dcreate(file, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int *_data = malloc(sizeof(int) * n);
    for (int i=0; i< n; i++) {
        _data[i] = data[i*stride];
    }
    H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _data);
    free(_data);
}
void write_h5_int_2D_array(void *f, char* name, int* data, int n, int m) {
    hid_t file = ((hid_t*) f)[0];
    hsize_t size[] = {n, m};
    hid_t space = H5Screate_simple(2, size, NULL);
    hid_t dset = H5Dcreate(file, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
}
void write_h5_complex_3D_array(void *f, char* name, double* data, int n, int m, int k) {
    hid_t file = ((hid_t*) f)[0];
    hsize_t size[] = {n, m, k, 2};
    hid_t space = H5Screate_simple(4, size, NULL);
    hid_t dset = H5Dcreate(file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
}

