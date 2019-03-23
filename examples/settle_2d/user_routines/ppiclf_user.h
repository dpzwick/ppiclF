C> Maximum number of particles per processor
#define PPICLF_LPART 50000

C> Max number of equations being solved for each particle
#define PPICLF_LRS 4

C> Pointers to PPICLF_Y(i,*) and PPICLF_YDOT(i,*)
C>    -CANNOT exceed PPICLF_LRS
C>    -Entirely user defined
C>    -X,Y,Z coordinates of particle MUST be ordered 1,2,3
#define PPICLF_JX  1
#define PPICLF_JY  2
#define PPICLF_JVX 3
#define PPICLF_JVY 4

C> Max number of additional properties for each particle
#define PPICLF_LRP 3

C> Pointers to PPICLF_RPROP(i,*)
C>    -CANNOT exceed PPICLF_LRP
C>    -Entirely user defined
#define PPICLF_R_JRHOP 1
#define PPICLF_R_JDP   2
#define PPICLF_R_JVOLP 3

C> Max size of local SubBin grid
#define PPICLF_BX1 50
#define PPICLF_BY1 50

C> Number of fields being projected
#define PPICLF_LRP_PRO 1

C> Pointers to projected field in PPICLF_PRO_FLD(*,*,*,*,i)
C>    -CANNOT exceed PPICLF_LRP_PRO
C>    -Entirely user defined
#define PPICLF_P_JPHIP 1   
