#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "simplemat.h"

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

void write_mat_header(FILE *f) {
    
    struct mat_header header;
    
    memset(header.text,0,sizeof(header.text));
    strcpy(header.text,"MATLAB 5.0 MAT-file, created by openmx-hks");
    memset(header.subsys,0,sizeof(header.subsys));
    header.version = 0x0100;
    header.endian = 0x4d49;
    fwrite(&header, 1, sizeof(header), f);
    
}

int calculate_size(int dims, int name_size, int array_size, int blocks) {
    int header = 16;
    int shape = PADDED(dims*4)+8;
    int name = PADDED(name_size)+8;
    int data = PADDED(array_size)+8;
    return  header + shape + name + data*blocks;
}

void write_mat_array_header(FILE *f, int total_size, int kind) {
    
    char d[24];
    
    memset(d,0,sizeof(d));
    
    *((int *)d) = miMATRIX;
    *((int *)(d+4)) = total_size;
    *((int *)(d+8)) = miUINT32;
    *((int *)(d+12)) = 8;
    *((int *)(d+16)) = kind;
    
    fwrite(d, sizeof(d), sizeof(char), f);
    
}


void write_mat_array_shape(FILE *f, int rank, int *dims) {
    
    int l = rank*4;
    int padded_l = PADDED(l);
    char d[padded_l+8];
    
    memset(d,0,sizeof(d));
    *((int *)(d)) = miINT32;
    *((int *)(d+4)) = l;
    memcpy(d+8,dims,l);
    
    fwrite(d, sizeof(d), sizeof(char), f);

}


void write_mat_name(FILE *f, char* name) {
    
    int l = strlen(name);
    int padded_l = PADDED(l);
    char d[padded_l+8];
    
    memset(d,0,sizeof(d));
    *((int *)(d)) = miINT8;
    *((int *)(d+4)) = l;
    memcpy(d+8,name,l);
    
    fwrite(d, sizeof(d), sizeof(char), f);
    
}

void write_mat_plain_double(FILE *f, double *data, int size, int step) {
    
    int l = size*8;
    char d[8];
    
    *((int *)(d)) = miDOUBLE;
    *((int *)(d+4)) = l;
    fwrite(d, sizeof(d), sizeof(char), f);
    
    int i;
    for (i=0; i<size; i++) {
        fwrite(data + i*step, 1, sizeof(double), f);
    }
    
}

void write_mat_plain_int(FILE *f, int *data, int size, int step) {
    
    int l = size*4;
    int padded_l = PADDED(l);
    char d[padded_l+8];
    
    memset(d,0,sizeof(d));
    *((int *)(d)) = miINT32;
    *((int *)(d+4)) = l;
    
    int i;
    for (i=0; i<size; i++) {
        *((int *)(d+4*i+8)) = data[i*step];
    }

    fwrite(d, sizeof(d), sizeof(char), f);
    
}

void write_mat_double_2D_array(FILE *f, char* name, double* data, int n, int m) {
    
    write_mat_array_header(f, calculate_size(2,strlen(name),n*m*sizeof(double),1), mxDOUBLE_CLASS);
    int dims[2] = {n,m};
    write_mat_array_shape(f,2,dims);
    write_mat_name(f,name);
    write_mat_plain_double(f,data,n*m,1);
    
}

void write_mat_double_scalar(FILE *f, char* name, double* data) {
    
    write_mat_double_2D_array(f, name, data, 1, 1);
    
}

void write_mat_complex_3D_array(FILE *f, char* name, double* data, int n, int m, int k) {
    
    write_mat_array_header(f, calculate_size(3,strlen(name),n*m*k*sizeof(double),2), mxDOUBLE_CLASS+0x0800);
    int dims[3] = {n,m,k};
    write_mat_array_shape(f,3,dims);
    write_mat_name(f,name);
    write_mat_plain_double(f,data,n*m*k,2);
    write_mat_plain_double(f,data+1,n*m*k,2);
    
}

void write_mat_int_2D_array(FILE *f, char* name, int* data, int n, int m) {
    
    write_mat_array_header(f, calculate_size(2,strlen(name),n*m*sizeof(int),1), mxINT32_CLASS);
    int dims[2] = {n,m};
    write_mat_array_shape(f,2,dims);
    write_mat_name(f,name);
    write_mat_plain_int(f,data,n*m,1);
    
}

void write_mat_int_1D_array(FILE *f, char* name, int* data, int n, int step) {
    
    write_mat_array_header(f, calculate_size(2,strlen(name),n*sizeof(int),1), mxINT32_CLASS);
    int dims[2] = {n,1};
    write_mat_array_shape(f,2,dims);
    write_mat_name(f,name);
    write_mat_plain_int(f,data,n,step);
    
}
