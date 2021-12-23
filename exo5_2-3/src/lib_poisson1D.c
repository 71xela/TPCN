/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void print_rowMajor(double *AB, int *lab, int *la)
{
    for(int i = 0; i < *lab; i++)
    {
        for(int j = 0; j < *la; j++)
            printf("%lf ", AB[j + i * *la]);
        printf("\n");
    }
}

void print_colMajor(double *AB, int *lab, int *la)
{
    for(int i = 0; i < *lab; i++)
    {
        for(int j = 0; j < *la; j++)
            printf("%lf ", AB[j + i * *lab]);
    }
    printf("\n");
}

void set_GB_operator_rowMajor_poisson1D(double *AB, int *lab, int *la, int *kv)
{
    for (int j = 0; j < *la; j++)
    {
        if(kv != 0)
            AB[j] = 0.0;
        AB[j + *kv * *la] = -1.0;
        AB[j + (*kv+1) * *la] = 2.0;
        AB[j + (*kv+2) * *la] = -1.0;
    }
    
    AB[*kv * *la] = 0.0;
    AB[*lab * *la - 1] = 0.0;
}

// Initialise matrice de Poisson 1D en colonne 1D
void set_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, int *kv)
{
    int kk;
    for (int j = 0; j < *la; j++){
        kk = j * *lab;
        if(*kv >= 0)
            for (int i = 0; i < *kv; i++)
                AB[kk+i] = 0.0;
        AB[kk + *kv] = -1.0;
        AB[kk + *kv +1] = 2.0;
        AB[kk + *kv +2] = -1.0;
    }
    
    AB[0] = 0.0;
    if(*kv == 1)
        AB[1] = 0;

    AB[*lab * *la - 1] = 0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double *AB, int *lab, int *la, int *kv)
{
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

// Initialisation vecteur f
void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1)
{
    RHS[0] = *BC0; //T0
    RHS[*la - 1]= *BC1; //T1
    for (int i = 1; i < *la - 1; i++)
        RHS[i] = 0.0;
}  

// Initialise solution analytique 
void set_analytical_solution_DBC_1D(double *EX_SOL, double *X, int *la, double *BC0, double *BC1)
{
    double DELTA_T;
    DELTA_T = *BC1 - *BC0;
    for (int i = 0; i < *la; i++)
        EX_SOL[i] = *BC0 + X[i]*DELTA_T;
}  

// Discrétisation 1D
void set_grid_points_1D(double *x, int *la)
{
    double h = 1.0/(1.0 * (*la + 1));
    for (int i = 0; i < (*la); i++)
        x[i] = (i + 1) * h;
}

void write_GB_operator_rowMajor_poisson1D(double *AB, int *lab, int *la, char *filename)
{
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, char *filename)
{
  //TODO
}

// Écrit le vecteur dans un fichier .dat
void write_vec(double *vec, int *la, char *filename)
{
    FILE * file;
    file = fopen(filename, "w");
    // Numbering from 1 to la
    if(file != NULL){
        for (int i = 0; i < *la; i++)
            fprintf(file,"%lf\n", vec[i]);
        fclose(file);
    }
    else
        perror(filename);
}  

void write_xy(double *vec, double *x, int *la, char *filename)
{
    FILE *file;
    file = fopen(filename, "w");
    // Numbering from 1 to la
    if (file != NULL){
        for (int i = 0; i < *la; i++)
            fprintf(file,"%lf\t%lf\n", x[i], vec[i]);
        fclose(file);
    }
    else
        perror(filename);
}  

void eig_poisson1D(double *eigval, int *la)
{
  int ii;
  double scal;
  for (ii=0; ii< *la; ii++){
    scal=(1.0*ii+1.0)*M_PI_2*(1.0/(*la+1));
    eigval[ii]=sin(scal);
    eigval[ii]=4*eigval[ii]*eigval[ii];
  } 
}

double eigmax_poisson1D(int *la)
{
  double eigmax;
  eigmax=sin(*la *M_PI_2*(1.0/(*la+1)));
  eigmax=4*eigmax*eigmax;
  return eigmax;
}

double eigmin_poisson1D(int *la)
{
  double eigmin;
  eigmin=sin(M_PI_2*(1.0/(*la+1)));
  eigmin=4*eigmin*eigmin;
  return eigmin;
}

double richardson_alpha_opt(int *la)
{
  //TODO
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit)
{
  //TODO
}
