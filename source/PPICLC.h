#include "PPICLF_USER.h"
#include "PPICLF_STD.h"

extern"C" {

    // ppiclf_comm.f
    void ppiclc_comm_InitMPI            (int     *comm   , 
                                         int     *id     , 
					 int     *np     );
    void ppiclc_comm_InitOverlapMesh    (int     *ncell  , 
                                         int     *lx1    , 
					 int     *ly1    , 
					 int     *lz1    , 
					 double  *xgrid  , 
					 double  *ygrid  , 
					 double  *zgrid  );

    // ppiclf_io.f
    void ppiclc_io_ReadParticleVTU       (char   *filein );
    void ppiclc_io_ReadWallVTK           (char   *filein );

    // ppiclf_solve.f
    void ppiclc_solve_InitNeighborBin    (double *width  );
    void ppiclc_solve_InitTargetBins     (char   *dir    ,
                                          int    *n      ,
                                          int    *balance);
    void ppiclc_solve_InitPeriodicX      (double *xl     , 
                                          double *xr     );
    void ppiclc_solve_InitPeriodicY      (double *yl     , 
                                          double *yr     );
    void ppiclc_solve_InitPeriodicZ      (double *zl     , 
                                          double *zr     );
    void ppiclc_solve_InitGaussianFilter (double *filt   , 
                                          double *alpha  , 
					  int    *ngrid  );
    void ppiclc_solve_InitBoxFilter      (double *filt   , 
                                          int    *ngrid  ,
					  int    *sngl_elem);
    void ppiclc_solve_InitParticle       (int    *imethod, 
                                          int    *ndim   , 
					  int    *iendian, 
					  int    *npart  , 
					  double *y      , 
					  double *rprop  );
    void ppiclc_solve_AddParticles       (int    *npart  , 
					  double *y      , 
					  double *rprop  );
    void ppiclc_solve_IntegrateParticle  (int    *istep  , 
                                          int    *iostep , 
					  double *dt     , 
					  double *time   );
}
