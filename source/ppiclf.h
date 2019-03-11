c Maximum number of real particles on a processor
#ifndef PPICLF_LPART
#define PPICLF_LPART 100000
#endif

c Number of secondary real properties for a particle
#define PPICLF_LRP2 4

c Number of integer properties for a particle
#ifndef PPICLF_LIP
#define PPICLF_LIP 11
#endif

c Maximum number of ghost particles
#define PPICLF_LPART_GP 26*PPICLF_LPART

#ifndef PPICLF_LRP_PRO
#define PPICLF_LRP_PRO 1
#endif
#if PPICLF_LRP_PRO == 0
#undef PPICLF_LRP_PRO
#define PPICLF_LRP_PRO 1
#endif

#ifndef PPICLF_LRP_GP
#define PPICLF_LRP_GP 3+PPICLF_LRP_PRO
#endif

#ifndef PPICLF_LELT
#define PPICLF_LELT 1
#endif

#ifndef PPICLF_LX1
#define PPICLF_LX1 1
#endif

#ifndef PPICLF_LY1
#define PPICLF_LY1 1
#endif

#ifndef PPICLF_LZ1
#define PPICLF_LZ1 1
#endif

c max bins per processor
#ifndef PPICLF_BMAX
#define PPICLF_BMAX 1
#endif

c max gridpts per processor
#ifndef PPICLF_BX1
#define PPICLF_BX1 1
#endif
#ifndef PPICLF_BY1
#define PPICLF_BY1 1
#endif
#ifndef PPICLF_BZ1
#define PPICLF_BZ1 1
#endif
