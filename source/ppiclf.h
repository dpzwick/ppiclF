C> Maximum number of particles per processor
#ifndef PPICLF_LPART
#define PPICLF_LPART 100000
#endif

C> Max number of additional properties for each particle
#ifndef PPICLF_LRP
#define PPICLF_LRP 1
#endif

C> Max number of secondary additional properties for each particle
C>    -User never touches
#define PPICLF_LRP2 4

C> Max number of integer properties for each particle
C>    -User never touches
#ifndef PPICLF_LIP
#define PPICLF_LIP 11
#endif

C> Max number of ghost particles per processor
C>    -User never touches
#define PPICLF_LPART_GP 26*PPICLF_LPART

C> Number of fields being interpolated
#ifndef PPICLF_LRP_INT
#define PPICLF_LRP_INT 1
#define PPICLF_INTERP 0
#else
#define PPICLF_INTERP 1
#endif

C> Number of fields being projected
#ifndef PPICLF_LRP_PRO
#define PPICLF_LRP_PRO 1
#define PPICLF_PROJECT 0
#else
#define PPICLF_PROJECT 1
#endif

C> Max number of real ghost particle properties
C>    -3 (X,Y,Z) + Projected properties
C>    -User never touches
#ifndef PPICLF_LRP_GP
#define PPICLF_LRP_GP 3+PPICLF_LRP_PRO
#endif

C> Max number of integer ghost particle properties
C>    -User never touches
#ifndef PPICLF_LIP_GP
#define PPICLF_LIP_GP 5
#endif

C> Max size of external overlap in x
#ifndef PPICLF_LEX
#define PPICLF_LEX 1
#endif

C> Max size of external overlap in y
#ifndef PPICLF_LEY
#define PPICLF_LEY 1
#endif

C> Max size of external overlap in z
#ifndef PPICLF_LEZ
#define PPICLF_LEZ 1
#endif

C> Max size of external overlap elements
#ifndef PPICLF_LEE
#define PPICLF_LEE 1
#endif

C> Max number of bins per processor
C>    -User never touches
#ifndef PPICLF_BMAX
#define PPICLF_BMAX 1
#endif

C> SubBins per processor in x
#ifndef PPICLF_BX1
#define PPICLF_BX1 1
#endif

C> SubBins per processor in y
#ifndef PPICLF_BY1
#define PPICLF_BY1 1
#endif

C> SubBins per processor in z
#ifndef PPICLF_BZ1
#define PPICLF_BZ1 1
#endif

C> Overlap grid max integer properties
C>    -User never touches
#ifndef PPICLF_LRMAX
#define PPICLF_LRMAX 6
#endif

