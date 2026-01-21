#pragma once

#include <iostream>
#include <fstream>
#include <complex>
#include <iomanip>
#include <time.h>
#include <stdlib.h>
#include <unordered_map>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "constants.h"

using complex = std::complex<double>;

using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using SparseMatrixXd = Eigen::SparseMatrix<double>;

#define do_once  {static int flag_once=0; if (flag_once==0) {flag_once=1;
#define end_do_once ;};}

#define for1(i, N) for (int i=0; i<N; i++) 
#define for2(i, j, N) for (int i=0; i<N; i++) for (int j=0; j<N; j++) 
#define for3(i, j, k, N) for (int i=0; i<N; i++) for (int j=0; j<N; j++) for (int k=0; k<N; k++)


struct LanczosTridiagRep {
	VectorXd a_vec;
	VectorXd b_vec;
	int krylov_dim;
};

int sgn(double x);
double sqr(double x);
complex sqr(const complex & x);

double norm2(const complex & x);

void swing(int & x, int & y);
void swing(double & x, double & y);
void swing(complex & x, complex & y);
void swing(double * & x, double * & y);
void swing(complex * & x, complex * & y);

int min(int & x, int &y);
double min(double x, double y);

int max(int & x, int &y);

int int_rnd();
double rnd ();
int rnd (int k);

complex rnd_gauss2();


int ** new_int2(int n1, int n2);
void delete_int2(int ** &r, int n1, int n2);

double ** new_double2(int n1, int n2);
void delete_double2(double ** &r, int n1, int n2);

double ** new_ldouble2(int n1, int n2);
void delete_ldouble2(double ** &r, int n1, int n2);

complex ** new_complex2(int n1, int n2);
void delete_complex2(complex ** &r, int n1, int n2);

int ***  new_int3(int n1, int n2, int n3);
void delete_int3(int *** &r, int n1, int n2, int n3);

double *** new_double3(int n1, int n2, int n3);
void delete_double3(double *** &r, int n1, int n2, int n3);

complex *** new_complex3(int n1, int n2, int n3);
void delete_complex3(complex *** &r, int n1, int n2, int n3);

int **** new_int4(int n1, int n2, int n3, int n4);
complex **** new_complex4(int n1, int n2, int n3, int n4);
void delete_complex4(complex **** &r, int n1, int n2, int n3, int n4);


void Inverse(complex ** &a, int size); //Gauss with partial pivoting
void div_left(complex ** a, complex *  r, int size); //a=r^{-1}*a, Gauss with partial pivoting;  destroys data in a
bool div_left(double ** a, double *  r, int size); //a=r^{-1}*a, Gauss with partial pivoting;    destroys data in a
bool div_left(double ** a, VectorXd  r, int size);

void EigenJacobi (double ** h, double ** a, double * e, int n, double accuracy=1e-10, int max_sweep=100);
void EigenBlock (double ** h, double ** a, double * e, int n, int * block, int n_block, double accuracy=1e-10, int max_sweep=100);

LanczosTridiagRep LanczosTridiag(SparseMatrixXd& H, VectorXd& v0, int max_Krylov=20, double stop_tol=1e-10);
