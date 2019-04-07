#ifndef PPICLF_LPART
#define PPICLF_LPART 100000
#endif

#ifndef PPICLF_LRP
#define PPICLF_LRP 1
#endif

#define PPICLF_LRP2 4

#ifndef PPICLF_LIP
#define PPICLF_LIP 11
#endif

#define PPICLF_LPART_GP 26*PPICLF_LPART

#ifndef PPICLF_INTERP
#define PPICLF_INTERP 1
#endif
#ifndef PPICLF_LRP_INT
#undef PPICLF_INTERP
#define PPICLF_INTERP 0
#define PPICLF_LRP_INT 1
#endif

#ifndef PPICLF_PROJECT
#define PPICLF_PROJECT 1
#endif
#ifndef PPICLF_LRP_PRO
#undef PPICLF_PROJECT
#define PPICLF_PROJECT 0
#define PPICLF_LRP_PRO 1
#endif

#ifndef PPICLF_LRP_GP
#define PPICLF_LRP_GP PPICLF_LRS+PPICLF_LRP+PPICLF_LRP_PRO
#endif

#ifndef PPICLF_LIP_GP
#define PPICLF_LIP_GP 5
#endif

#ifndef PPICLF_LEX
#define PPICLF_LEX 1
#endif

#ifndef PPICLF_LEY
#define PPICLF_LEY 1
#endif

#ifndef PPICLF_LEZ
#define PPICLF_LEZ 1
#endif

#ifndef PPICLF_LEE
#define PPICLF_LEE 1
#endif

#ifndef PPICLF_BMAX
#define PPICLF_BMAX 1
#endif

#ifndef PPICLF_BX1
#define PPICLF_BX1 1
#endif

#ifndef PPICLF_BY1
#define PPICLF_BY1 1
#endif

#ifndef PPICLF_BZ1
#define PPICLF_BZ1 1
#endif

#ifndef PPICLF_LRMAX
#define PPICLF_LRMAX 6
#endif

#ifndef PPICLF_LWALL
#define PPICLF_LWALL 20
#endif

