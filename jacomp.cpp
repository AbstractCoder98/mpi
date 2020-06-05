#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include <iostream>
#include <omp.h>

using namespace std;

int hilos = 4;
//#include "timing.h"

/* --
 * Realice barridos de iteración de Jacobi en un problema de Poisson 1D
 * *
 * -u '' = f
 * *
 * discretizado por n + 1 puntos de malla igualmente espaciados en [0,1].
 * u está sujeto a las condiciones de contorno de Dirichlet especificadas en
 * las entradas u [0] yu [n] del vector inicial.

 */
//jacobi(nsteps, n, u, f);
void jacobi(int nsweeps, int n, double* u, double* f)
{
    int i, sweep;
    double h  = 1.0 / n;
    double h2 = h*h;
    double* utmp = (double*) malloc( (n+1) * sizeof(double) );

    /* Rellenar condiciones de contorno en utmp */
    utmp[0] = u[0];
    utmp[n] = u[n];

    for (sweep = 0; sweep < nsweeps; sweep += 2) {
        #pragma omp parallel for
            /* Datos antiguos en u; nuevos datos en utmp */
            for (i = 1; i < n; ++i){
                utmp[i] = (u[i-1] + u[i+1] + h2*f[i])/2;
		}
            /* Datos antiguos en utmp; nuevos datos en u */
            for (i = 1; i < n; ++i){
                u[i] = (utmp[i-1] + utmp[i+1] + h2*f[i])/2;
		}

    }

    free(utmp);
}

int main(int argc, char** argv)
{
    int i;
    int n, nsteps;
    double* u;
    double* f;
    double h;
    unsigned tstart, tend;
    //char* fname;

    int prueba = atoi(argv[1]);

    /* Procesar argumentos */
    n      = 10000;
    nsteps = 10000;
    //fname  = (argc > 3) ? argv[3] : NULL;
    h      = 1.0/n;

    /* Asignar e inicializar matrices */
    u = (double*) malloc( (n+1) * sizeof(double) );
    f = (double*) malloc( (n+1) * sizeof(double) );
    //llenamos U con ceros
    memset(u, 0, (n+1) * sizeof(double));
    for (i = 0; i <= n; ++i)
        f[i] = i * h;

    /*
    cout << "\nImprimiendo matriz Antes:\n";
    for (int i = 0; i < n; ++i)
    {
        cout << "\t" << *(u+i);
        cout << endl;
    }

    /* Ejecute el solucionador */
    tstart = clock();
    jacobi(nsteps, n, u, f);
    tend = clock();

    //Calculo del tiempo de ejecución del programa
    float totalTime = (double(tend-tstart)/CLOCKS_PER_SEC);

    /*
    cout << "\nImprimiendo matriz Despues:\n";
    for (int i = 0; i < n; ++i)
    {
        cout << "\t" << *(u+i);
        cout << endl;
    }
    */
    //cout << "\n-----------------------------------------------------:\n";

    /* Ejecute el solucionador */    
    /*printf("Prueba: %d\n"
            "n: %d\n"
           "nsteps: %d\n"
           "Elapsed time: %f s\n", 
           prueba, n, nsteps, totalTime);*/
    /* Ejecute el solucionador */    
    cout << totalTime << endl;

    /* Escribe los resultados */
    //if (fname)
    //   write_solution(n, u, fname);

    free(f);
    free(u);
    return 0;
}
