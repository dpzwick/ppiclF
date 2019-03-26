C> Maximum number of particles per processor
#define PPICLF_LPART 50000

C> Max number of equations being solved for each particle
#define PPICLF_LRS 6

C> Pointers to PPICLF_Y(i,*) and PPICLF_YDOT(i,*)
C>    -CANNOT exceed PPICLF_LRS
C>    -Entirely user defined
C>    -X,Y,Z coordinates of particle MUST be ordered 1,2,3
#define PPICLF_JX  1
#define PPICLF_JY  2
#define PPICLF_JZ  3
#define PPICLF_JVX 4
#define PPICLF_JVY 5
#define PPICLF_JVZ 6

C> Max number of additional properties for each particle
#define PPICLF_LRP 3

C> Pointers to PPICLF_RPROP(i,*)
C>    -CANNOT exceed PPICLF_LRP
C>    -Entirely user defined
#define PPICLF_R_JRHOP 1
#define PPICLF_R_JDP   2
#define PPICLF_R_JVOLP 3
