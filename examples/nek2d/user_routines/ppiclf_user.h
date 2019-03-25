C Auto-generated ppiclf_user.h file with usergen.py at:
C    2019-03-25 09:49:35.745897
C 
C Max number of particles per rank 
#define PPICLF_LPART 20000

C Number of equations solved for each particle 
#define PPICLF_LRS 4

C Naming of equations solved for each particle 
#define PPICLF_JX 1
#define PPICLF_JY 2
#define PPICLF_JVX 3
#define PPICLF_JVY 4

C Number of additional properties for each particle 
#define PPICLF_LRP 4

C Naming of additional properties for each particle 
#define PPICLF_R_JRHOP 1
#define PPICLF_R_JDP 2
#define PPICLF_R_JVOLP 3
#define PPICLF_R_JPHIP 4

C Size of external overlap grid 
#define PPICLF_LEX 6
#define PPICLF_LEY 6
#define PPICLF_LEE 25

C Number of interpolated fields 
#define PPICLF_LRP_INT 1

C Number of projected fields 
#define PPICLF_LRP_PRO 1

C Naming of projected fields 
#define PPICLF_P_JPHIP 1

C Max size of SubBin grid 
#define PPICLF_BX1 50
#define PPICLF_BY1 50

