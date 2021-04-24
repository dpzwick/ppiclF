#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
c Particle options
      LOGICAL PPICLF_RESTART, PPICLF_OVERLAP, PPICLF_LCOMM
     >       ,PPICLF_LINIT, PPICLF_LFILT, PPICLF_LINTP, PPICLF_LPROJ
     >       ,PPICLF_LSUBBIN, PPICLF_LSUBSUBBIN
     >       ,PPICLF_LFILTGAUSS, PPICLF_LFILTBOX
      COMMON /PPICLF_OPT_PARAM_L/ PPICLF_RESTART, PPICLF_OVERLAP
     >                           ,PPICLF_LCOMM, PPICLF_LINIT
     >                           ,PPICLF_LFILT, PPICLF_LINTP
     >                           ,PPICLF_LPROJ, PPICLF_LSUBBIN
     >                           ,PPICLF_LSUBSUBBIN
     >                           ,PPICLF_LFILTGAUSS, PPICLF_LFILTBOX
      DATA PPICLF_LCOMM /.false./
      DATA PPICLF_RESTART /.false./

      INTEGER*4 PPICLF_NDIM, PPICLF_IMETHOD, PPICLF_IPERIODIC(3)
     >         ,PPICLF_NGRIDS, PPICLF_CYCLE, PPICLF_IOSTEP
     >         ,PPICLF_IENDIAN, PPICLF_IWALLM, PPICLF_IOCALLD
      COMMON /PPICLF_OPT_PARAM_I/ PPICLF_NDIM, PPICLF_IMETHOD
     >                           ,PPICLF_IPERIODIC, PPICLF_NGRIDS
     >                           ,PPICLF_CYCLE, PPICLF_IOSTEP
     >                           ,PPICLF_IENDIAN, PPICLF_IWALLM
     >                           ,PPICLF_IOCALLD

      REAL*8 PPICLF_FILTER, PPICLF_ALPHA, PPICLF_RK3COEF(3,3), PPICLF_DT
     >      ,PPICLF_TIME, PPICLF_D2CHK(3)
     >      ,PPICLF_PREV_TIMES(PPICLF_PREV_MAX_TIMES)
      COMMON /PPICLF_OPT_PARAM_R/ PPICLF_FILTER, PPICLF_ALPHA
     >                           ,PPICLF_RK3COEF, PPICLF_DT
     >                           ,PPICLF_TIME, PPICLF_D2CHK
     >                           ,PPICLF_PREV_TIMES
