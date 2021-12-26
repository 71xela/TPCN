/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

int main(int argc, char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
    //Au = f
    int nbpoints, la; //la = n
    int ku, kl, kv, lab; //kl sous-diagonales, 
                         //ku sur-diagonales, 
                         //kv = ligne en plus dans AB car utiliser
                         //pour stocker L et U
                         //lab = nb de lignes de AB
    int *ipiv; // tableau des pivot de la factorisation LU
    int info;  // Retour de dgbsv
    int NRHS; // nb de colonne de f
    double T0, T1; // Condition initiales T(0) = T0 et T(1) = T1
    double *RHS, *EX_SOL, *X; //RHS = f, EX_SOL = solution analytique, grille 1D des valeurs entre T0 et T1 espacé d'un pas constant
    double *AB; //matrice General Band

    double temp, relres; // utiliser pour le calcul de l'erreur relative

    NRHS = 1;
    nbpoints = 102;
    la = nbpoints-2;
    T0 = -5.0;
    T1 = 5.0;

    printf("--------- Poisson 1D ---------\n\n");
    RHS = (double*)malloc(sizeof(double) * la);
    EX_SOL = (double*)malloc(sizeof(double) * la);
    X = (double*)malloc(sizeof(double) * la);

    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    write_vec(RHS, &la, "RHS.dat");
    write_vec(EX_SOL, &la, "EX_SOL.dat");
    write_vec(X, &la, "X_grid.dat");

    kv = 1;
    ku = 1;
    kl = 1;
    lab = kv+kl+ku+1;

    AB = (double*)malloc(sizeof(double) * lab * la);

    info = 0;

    /* working array for pivot used by LU Factorization */
    ipiv = (int*)calloc(la, sizeof(int));

    int row = 1; //

    if(row == 1) // LAPACK_ROW_MAJOR
    {
        set_GB_operator_rowMajor_poisson1D(AB, &lab, &la, &kv);
        write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");
        //print_rowMajor(AB, &lab, &la);
        info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, la, kl, ku, NRHS, AB, la, ipiv, RHS, NRHS);
    } 
    else // LAPACK_COL_MAJOR
    { 
        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");
        //print_colMajor(AB, &lab, &la);
        info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);
    }    

    printf("\n INFO DGBSV = %d\n", info);

    write_xy(RHS, X, &la, "SOL.dat");

    /* Relative residual */
    temp = cblas_ddot(la, RHS, 1, RHS,1); // dotprod entre vecteur
    temp = sqrt(temp);
    cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1); // daxpy
    relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
    relres = sqrt(relres);
    relres = relres / temp;

    printf("\nThe relative residual error is relres = %e\n\n\n",relres);
    
    print_rowMajor(AB, &lab, &la);

    

    // kv = 1;
    // ku = 1;
    // kl = 1;
    // lab = kv+kl+ku+1;
    // set_GB_operator_rowMajor_poisson1D(AB, &lab, &la, &kv);

    // print_rowMajor(AB, &lab, &la);
    // info = LAPACKE_dgbtrf(LAPACK_ROW_MAJOR, la, la, kl, ku, AB, lab, ipiv);
    // printf("\n");
    // print_rowMajor(AB, &lab, &la);

    // printf("\n INFO DGBTRF = %d\n", info);
    




    free(RHS);
    free(EX_SOL);
    free(X);
    free(AB);
    free(ipiv);

    printf("\n\n--------- End -----------\n");

    /*
        
    */

    return 0;
}
