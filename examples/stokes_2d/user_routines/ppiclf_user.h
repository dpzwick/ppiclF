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
#define PPICLF_LRP 1

C> Pointers to PPICLF_RPROP(i,*)
C>    -CANNOT exceed PPICLF_LRP
C>    -Entirely user defined
#define PPICLF_R_JTAUP 1
