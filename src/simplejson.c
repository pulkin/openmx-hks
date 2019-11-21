#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "simplejson.h"

#define JSON_NAN -1.0
#define JSON_NEWLINE do {fprintf(f,"\n"); int __i; for (__i=0; __i<indent; __i++) fprintf(f," ");} while(0)
#define JSON_COMMA if (comma) fprintf(f,",")
#define JSON_O_S do {JSON_COMMA; fprintf(f,"{"); indent++; comma = 0;} while(0)
#define JSON_O_E do {fprintf(f,"}"); indent--; comma = 1;} while(0)
#define JSON_A_S do {JSON_COMMA; fprintf(f,"["); indent++; comma = 0;} while(0)
#define JSON_A_E do {fprintf(f,"]"); indent--; comma = 1;} while(0)
#define JSON_S_KEY(k) do {JSON_COMMA; fprintf(f,"\"%s\":",k); comma = 0;} while(0)
#define JSON_S(s) do {JSON_COMMA; fprintf(f,"\"%s\"",s); comma = 1;} while(0)
#define JSON_D(d) do {JSON_COMMA; fprintf(f,"%.16e",d == d ? d : JSON_NAN); comma = 1;} while(0)
#define JSON_I(i) do {JSON_COMMA; fprintf(f,"%d",i); comma = 1;} while(0)

int indent, comma = 0;

void *open_json(char* fname) {return fopen(fname, "w");}
void close_json(void *f) {fclose(f);}
void write_json_header(void *f) {indent = 0; comma = 0; JSON_O_S;}
void write_json_footer(void *f) {JSON_O_E;}
void write_json_double_scalar(void *f, char* name, double* data) {JSON_S_KEY(name); JSON_D(*data);}
void write_json_int_1D_array(void *f, char* name, int* data, int n, int stride) {
    int i=0;
    JSON_S_KEY(name);
    JSON_A_S;
    for (i=0; i<n; i++) {JSON_I(data[stride*i]);}
    JSON_A_E;
}
void write_json_int_2D_array(void *f, char* name, int* data, int n, int m) {
    int i,j=0;
    JSON_S_KEY(name);
    JSON_A_S;
    for (i=0; i<n; i++) {
        JSON_A_S;
        for (j=0; j<m; j++) {
            JSON_I(data[i*m+j]);
        }
        JSON_A_E;
    }
    JSON_A_E;
}
void write_json_complex_3D_array(void *f, char* name, double* data, int n, int m, int k) {
    int i,j,l,c=0;
    JSON_S_KEY(name);

    JSON_A_S;
    for (i=0; i<n; i++) {
        JSON_A_S;
        for (j=0; j<m; j++) {
            JSON_A_S;
            for (l=0; l<k; l++) {
                JSON_A_S;
                for (c=0; c<2; c++) {
                    JSON_D(data[((i*m+j)*k+l)*2+c]);
                }
                JSON_A_E;
            }
            JSON_A_E;
        }
        JSON_A_E;
    }
    JSON_A_E;
}
