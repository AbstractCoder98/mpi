#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <ctime>
#include <iostream>
#include <omp.h>

using namespace std;


//En este apartado es donde enviaremos la informacion a cada uno de los nodos,
//y a su vez seran recibidos para reunificar la informacion
void cambiarDatos(double* u, int n, int rank, int size)
{
    if (size == 1)
        return;

    //Definimos la variable de estado para dar a conocer al sistema el estado de
    //los mensajes.
    MPI_Status status;

    //Si estamos en el primer proceso => master
    if (rank == 0){
        // Sendrecv(que dato quiero enviar, cuantos datos quiero enviar, tipo de dato, 
        //a donde lo quiero enviar, tag == 0 [HASTA AQUI ES COMO SI HICIERAMOS UN SEND], 
        //que quiero recibir, cuantos datos quiero recibir, el tipo de dato que se recibe, 
        //a quien le recibo, tag, la parte del comunicador, el status)
        MPI_Sendrecv(&u[n-1], 1, MPI_DOUBLE, rank + 1, 0, &u[n], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
    }
    //Si estamos en el ultimo proceso => Wn3, entonces:
    else if (rank == size - 1){
        MPI_Sendrecv(&u[1], 1, MPI_DOUBLE, rank - 1, 0, &u[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
    }

    //Si no, se supone que estamos tratando con los espacios de la mitad
    else{
        MPI_Sendrecv(&u[1], 1, MPI_DOUBLE, rank - 1, 0, &u[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(&u[n-1], 1, MPI_DOUBLE, rank + 1, 0, &u[n], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
    }
}

void jacobi(int nsweeps, int n, double* u, double* f, double h2, 
            int rank, int size){
    //defino unas variables para los ciclos
    int i, sweep;
    //recupero el valor de h
    double h  = 1.0 / n;

    //generamos el vector utmp para manejarlo como mediador entre el vector u y f,
    //este vector se encarga de guardar los valores que iran al nuevo u de forma
    //temporal; en este punto si dedicamos todo el espacio del vector (n)
    double* utmp = (double*) malloc( (n+1) * sizeof(double) );

    //defino las variables de contexto las cuales se mantendran quietas, al inicio
    //y al final de utmp, en las demas posiciones guardaremos los valores de cambio
    utmp[0] = u[0];
    utmp[n] = u[n];

    //realizamos el ciclo para cambiar de viejo u a nuevo u
    for (sweep = 0; sweep < nsweeps; sweep += 2) {
	cambiarDatos(u, n, rank, size);
        utmp[0] = u[0];
        utmp[n] = u[n];

        //empiezo a llenar utmp con respecto a lo que tengo en el viejo u y en f
	#pragma omp parallel for
    		for (i = 1; i < n; ++i)
               		utmp[i] = (u[i-1] + u[i+1] + h2*f[i])/2;

        /* Exchange ghost cells */
        cambiarDatos(utmp, n, rank, size);
        u[0] = utmp[0];
        u[n] = utmp[n];

        //los datoss viejos de utmp los empiezo a pasar al nuevo u
        #pragma omp parallel for
		for (i = 1; i < n; ++i)
        		u[i] = (utmp[i-1] + utmp[i+1] + h2*f[i])/2;
    }

    //libero la memoria donde se encontraba utmp
    free(utmp);
}

int main(int argc, char** argv)
{
    //Definimos las variables de conteo
    int i;
    int n, nsteps;
    //declaramos los dos vectores que van a alojar informaciòn
    double* uloc;
    double* floc;
    double h;
    //definimos las variables para tomar el tiempo
    double tstart, tend;
    char* fname;
    //declaro la variable que alojara el rango del procesador donde nos encontremos,
    //ademas de la variable size que cuenta cuantos procesos tenemos
    int rank, size;
    //declaro las variables para saber con que cantidad de datos dividire mis vectores
    //para aplicar jacobi (nper), ademas defino donde empezara a leer en u y f (ioffset)
    //y donde terminara mi lectura (nloc), que por lo general define el tamanio de los
    //subvectores
    int ioffset, nper, nloc;

    /* INICIALIZAMOS MPI */
    MPI_Init(&argc, &argv);

    // Inicializo las variables de rank (que define en que proceso me encuentro 
    // actualmente) y size (para identificar la cantidad de procesos).

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //declaro las variables que llegaran desde consola de tamanio para los vectores (n), 
    //el numero de cambios que hare al vector (nsteps), un archivo para guardar info
    //(fname), y h lo cual define un dato 
    n      = (argc > 1) ? atoi(argv[1]) : 100;
    nsteps = (argc > 2) ? atoi(argv[2]) : 100;
    fname  = (argc > 3) ? argv[3] : NULL;
    h      = 1.0/n;

    /* Print a diagnostic message */
    if (rank == 0)
	cout << "Procesos: " << size << endl;
        cout << "Cantidad: " << n << endl;
	cout << "Sweeps: " << nsteps << endl;

    //defino entre cuantos pedazos separare mi vector para que sea trabajado por los
    //nodos del cluster, ademas defino desde que punto iniciara el procesamiento 
    //(ioffset) y cuando dara fin (nloc)
    nper    = (n+size-1)/size;
    ioffset = rank*nper;
    nloc    = (rank == size-1) ? (n-ioffset) : nper+1;

    //designo memoria para los pedazos a trabajar para cada proceso, por lo que le doy
    //espacio en memoria a u y f
    uloc = (double*) malloc( (nloc+1) * sizeof(double) );
    floc = (double*) malloc( (nloc+1) * sizeof(double) );

    //ademas genero datos en esos espacios para u y f, en el caso de u lo llenamos de 
    //ceros y a f le asignamos otros valores
    memset(uloc, 0, (nloc+1) * sizeof(double));
    for (i = 0; i <= nloc; ++i)
        floc[i] = (ioffset+i) * h;

    //Desde este punto correremos el algoritmo de Jacobi el cual recibira (cuantos 
    //cambios haremos a u, tope del subvector, pedazo de u, pedazo de f, el valor 
    //de h*h, en que proceso estamos parado, y la cantidad total de procesos)
    tstart = MPI_Wtime();
    jacobi(nsteps, nloc, uloc, floc, h*h, rank, size);
    tend = MPI_Wtime();

    //Calculo del tiempo de ejecución del programa
    float totalTime = (double(tend-tstart)/CLOCKS_PER_SEC);
    cout << totalTime << endl;

    //finalmente liberamos la memoria dedicada al vector f y u
    free(floc);
    free(uloc);

    MPI_Finalize();
    return 0;
}
