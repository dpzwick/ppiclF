#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <mpi.h>
#include <PPICLC.h>

int main(int argc, char *argv[])
{

    int np, nid, icomm;
    int imethod, ndim, iendian, npart;
    double y    [PPICLF_LPART][PPICLF_LRS]; // Opposite ordering
    double rprop[PPICLF_LPART][PPICLF_LRP]; // Opposite ordering

    int nstep, iostep;
    double dt, time;
    int ierr;

    // Init MPI
    ierr  = MPI_Init(&argc, &argv);
    icomm = MPI_Comm_c2f(MPI_COMM_WORLD);
    ierr  = MPI_Comm_rank(MPI_COMM_WORLD, &nid);
    ierr  = MPI_Comm_size(MPI_COMM_WORLD, &np);

    // Pass to library to Init MPI
    ppiclc_comm_InitMPI(&icomm, 
                        &nid  , 
			&np   );

    // Set initial conditions and parameters for particles
    imethod = 1;
    ndim    = 2;
    iendian = 0;
    npart   = 250;
    srand(nid+1); // init random numbers
    for(int i=0; i<npart; i++)
    {
       y[i][PPICLF_JX-1]  = rand() / double(RAND_MAX);
       y[i][PPICLF_JY-1]  = rand() / double(RAND_MAX);
       y[i][PPICLF_JVX-1] = 0.0;
       y[i][PPICLF_JVY-1] = 0.0;

       rprop[i][PPICLF_R_JTAUP-1] = 1.0/9.8;
    }
    ppiclc_solve_InitParticle(&imethod    , 
                              &ndim       , 
			      &iendian    , 
			      &npart      , 
			      &y[0][0]    , 
			      &rprop[0][0]);

    // Integrate particles in time
    nstep  = 1000;
    iostep = 100;
    dt     = 1E-4;
    for(int istep=1; istep<=nstep; istep++)
    {
       time = (istep-1)*dt;
       ppiclc_solve_IntegrateParticle(&istep ,
                                      &iostep,
				      &dt    ,
				      &time  );
    }

    // Finalize MPI
    ierr = MPI_Finalize();

    return 0;
}
