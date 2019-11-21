#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "simplemat.h"

#define miMATRIX 14
#define miINT8 1
#define miDOUBLE 9
#define miINT32 5
#define miUINT32 6
#define miINT64 12
#define mxINT64_CLASS 14
#define mxDOUBLE_CLASS 6

#define PADDED(X) ((((X)-1)/8+1)*8)

struct mat_header {
    char text[116];
    char subsys[8];
    short version;
    short endian;
};

struct mat_data_header {
    int type;
    int size;
};

void *open_mat(char* fname) {return fopen(fname, "w");}
void close_mat(void *f) {fclose(f);}

void write_mat_header(void *f) {
    
    struct mat_header header;
    
    memset(header.text,0,sizeof(header.text));
    strcpy(header.text,"MATLAB 5.0 MAT-file, created by openmx-hks");
    memset(header.subsys,0,sizeof(header.subsys));
    header.version = 0x0100;
    header.endian = 0x4d49;
    fwrite(&header, 1, sizeof(header), f);
    
}

void write_mat_footer(void *f) {}

int calculate_size(int dims, int name_size, int array_size, int blocks) {
    int header = 16;
    int shape = PADDED(dims*4)+8;
    int name = PADDED(name_size)+8;
    int data = PADDED(array_size)+8;
    return  header + shape + name + data*blocks;
}

void write_mat_array_header(void *f, int total_size, int kind) {
    
    char d[24];
    
    memset(d,0,sizeof(d));
    
    *((int *)d) = miMATRIX;
    *((int *)(d+4)) = total_size;
    *((int *)(d+8)) = miUINT32;
    *((int *)(d+12)) = 8;
    *((int *)(d+16)) = kind;
    
    fwrite(d, sizeof(d), sizeof(char), f);
    
}


void write_mat_array_shape(void *f, int rank, int *dims) {
    
    int l = rank*4;
    int padded_l = PADDED(l);
    char d[padded_l+8];
    
    memset(d,0,sizeof(d));
    *((int *)(d)) = miINT32;
    *((int *)(d+4)) = l;
    memcpy(d+8,dims,l);
    
    fwrite(d, sizeof(d), sizeof(char), f);

}


void write_mat_name(void *f, char* name) {
    
    int l = strlen(name);
    int padded_l = PADDED(l);
    char d[padded_l+8];
    
    memset(d,0,sizeof(d));
    *((int *)(d)) = miINT8;
    *((int *)(d+4)) = l;
    memcpy(d+8,name,l);
    
    fwrite(d, sizeof(d), sizeof(char), f);
    
}

void _write_mat_int(void *f, int *data, u_int8_t ndims, u_int32_t *dims, u_int32_t stride) {
    if (ndims == 1) {
        for (int i=0; i<dims[0]; i++) {
            int64_t x = data[i*stride];
            fwrite(&x, 1, sizeof(int64_t), f);
        }
    } else {
        for (int i=0; i<dims[ndims-1]; i++) {
            _write_mat_int(f, data+stride*i, ndims-1, dims, stride*dims[ndims-1]);
        }
    }
}

void write_mat_plain_int(void *f, int *data, u_int8_t ndims, u_int32_t *dims, u_int32_t stride) {
    
    u_int32_t l = 1;
    for (int i=0; i< ndims; i++) {
        l = l * dims[i];
    }
    char d[8];
    
    *((int *)(d)) = miINT64;
    *((int *)(d+4)) = l * 8;
    fwrite(d, sizeof(d), sizeof(char), f);
    
    _write_mat_int(f, data, ndims, dims, stride);
}

void _write_mat_double(void *f, double *data, u_int8_t ndims, u_int32_t *dims, u_int32_t stride) {
    if (ndims == 1) {
        for (int i=0; i<dims[0]; i++) {
            fwrite(data + i*stride, 1, sizeof(double), f);
        }
    } else {
        for (int i=0; i<dims[ndims-1]; i++) {
            _write_mat_double(f, data+stride*i, ndims-1, dims, stride*dims[ndims-1]);
        }
    }
}

void write_mat_plain_double(void *f, double *data, u_int8_t ndims, u_int32_t *dims, u_int32_t stride) {
    
    u_int32_t l = 1;
    for (int i=0; i< ndims; i++) {
        l = l * dims[i];
    }
    char d[8];
    
    *((int *)(d)) = miDOUBLE;
    *((int *)(d+4)) = l * 8;
    fwrite(d, sizeof(d), sizeof(char), f);
    
    _write_mat_double(f, data, ndims, dims, stride);
}

void write_mat_double_2D_array(void *f, char* name, double* data, int m, int n) {
    
    write_mat_array_header(f, calculate_size(2,strlen(name),n*m*sizeof(double),1), mxDOUBLE_CLASS);
    int dims[2] = {n,m};
    write_mat_array_shape(f,2,dims);
    write_mat_name(f,name);
    write_mat_plain_double(f,data,2,dims,1);
    
}

void write_mat_double_scalar(void *f, char* name, double* data) {
    
    write_mat_double_2D_array(f, name, data, 1, 1);
    
}

void write_mat_complex_3D_array(void *f, char* name, double* data, int k, int m, int n) {
    write_mat_array_header(f, calculate_size(3,strlen(name),k*m*n*sizeof(double),2), mxDOUBLE_CLASS+0x0800);
    int dims[3] = {k,m,n};
    write_mat_array_shape(f,3,dims);
    write_mat_name(f,name);
    write_mat_plain_double(f,data,3,dims,2);
    write_mat_plain_double(f,data+1,3,dims,2); 
}

void write_mat_int_2D_array(void *f, char* name, int* data, int m, int n) {
    
    write_mat_array_header(f, calculate_size(2,strlen(name),n*m*sizeof(long int),1), mxINT64_CLASS);
    int dims[2] = {n,m};
    write_mat_array_shape(f,2,dims);
    write_mat_name(f,name);
    write_mat_plain_int(f,data,2,dims,1);
    
}

void write_mat_int_1D_array(void *f, char* name, int* data, int n, int step) {
    
    write_mat_array_header(f, calculate_size(2,strlen(name),n*sizeof(long int),1), mxINT64_CLASS);
    int dims[2] = {n,1};
    write_mat_array_shape(f,2,dims);
    write_mat_name(f,name);
    write_mat_plain_int(f,data,2,dims,step);
    
}
