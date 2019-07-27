#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <math.h>
#include <mpi.h>
#include <PPICLC.h>

// The following links the common block /ucollision/ ksp, erest between c and fortran
extern"C" {
    extern struct{
       double ksp;
       double erest;
    } ucollision;
}

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
    ndim    = 3;
    iendian = 0;
    npart   = 200;
    srand(nid+1); // init random numbers
    double dp_min = 0.05;
    double dp_max = 0.07;
    double rhop   = 2500.0;
    for(int i=0; i<npart; i++)
    {
       y[i][PPICLF_JX-1]  = -0.2 + 0.4*rand() / double(RAND_MAX);
       y[i][PPICLF_JY-1]  = -0.2 + 0.6*rand() / double(RAND_MAX);
       y[i][PPICLF_JZ-1]  = -0.2 + 0.4*rand() / double(RAND_MAX);
       y[i][PPICLF_JVX-1] = 0.0;
       y[i][PPICLF_JVY-1] = 0.0;
       y[i][PPICLF_JVZ-1] = 0.0;

       rprop[i][PPICLF_R_JRHOP-1] = rhop;
       rprop[i][PPICLF_R_JDP-1]   = dp_min + (dp_max-dp_min)*rand() / double(RAND_MAX);
       rprop[i][PPICLF_R_JVOLP-1] = M_PI / 6.0 * pow(rprop[i][PPICLF_R_JDP-1],3.0);
    }
    ppiclc_solve_InitParticle(&imethod    , 
                              &ndim       , 
			      &iendian    , 
			      &npart      , 
			      &y[0][0]    , 
			      &rprop[0][0]);

    // Set min bin size to be largest particle diameter
    ppiclc_solve_InitNeighborBin(&dp_max);

    // Read file with boundary conditions
    char bndry[16+1] = "ppiclf_tank.vtk";
    ppiclc_io_ReadWallVTK(bndry);

    // For user implemented collision model
    ucollision.ksp   = 10000.0;
    ucollision.erest = 0.1;
    double rmij1     = M_PI/6.0*pow(dp_min,3.0)*rhop;
    int nres         = 20;
    double rmij      = 1.0/(1.0/rmij1 + 1.0/rmij1);
    double dt_c_max  = sqrt(rmij/ucollision.ksp*(pow(log(ucollision.erest),2) + pow(M_PI,2.0)))/nres;

    // Integrate particles in time
    nstep  = 20000;
    iostep = 25;
    dt     = dt_c_max;
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
