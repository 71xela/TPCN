#include "lib_poisson1D.h"

void lu_tridiag(double *AB, int *lab, int *la, int *kv)
{
	for (int j = 1; j < *la; j++)
    {
        AB[j-1 + (*kv+2) * *la] = AB[j-1 + (*kv+2) * *la]/AB[j-1 + (*kv+1) * *la];
        AB[j + (*kv+1) * *la] = AB[j + (*kv+1) * *la] - AB[j-1 + (*kv+2) * *la] * AB[j + *kv * *la];
    }
}

int main(int argc, char const *argv[])
{
	int la = 100; //la = n
    int ku, kl, kv, lab; //kl sous-diagonales, 
                         //ku sur-diagonales, 
                         //kv = ligne en plus dans AB car utiliser
                         //pour stocker L et U
                         //lab = nb de lignes de AB
    double *AB, *AB_copy; //matrice General Band
    int info;  // Retour de dgbtrf
    int *ipiv; // tableau des pivot de la factorisation LU

    kv = 1;
    ku = 1;
    kl = 1;
    lab = ku + kl + kv + 1;

    AB = (double*)malloc(sizeof(double) * lab * la);
    AB_copy = (double*)malloc(sizeof(double) * lab * la);

    /* working array for pivot used by LU Factorization */
    ipiv = (int*)calloc(la, sizeof(int));

    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la, &kv);
    lu_tridiag(AB, &lab, &la, &kv);

    set_GB_operator_rowMajor_poisson1D(AB_copy, &lab, &la, &kv);
    info = LAPACKE_dgbtrf(LAPACK_ROW_MAJOR, la, la, kl, ku, AB_copy, la, ipiv);
    
    printf("\n INFO DGBTRF = %d\n", info);

    for(int i = 0; i < la * lab; i++)
    	AB[i] -= AB_copy[i];

    print_rowMajor(AB, &lab, &la);

	return 0;
}