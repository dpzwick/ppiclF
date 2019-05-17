import time
import datetime
import sys

# Get User Input
lpart = 100000
lpartin = int(raw_input("Max number of particles per rank [" + str(lpart) + "]: ") or lpart)

ndim = 2
ndim = int(raw_input("Problem dimensions [" + str(ndim) + "]: ") or ndim)

lrs = ndim
lrs = int(raw_input("Number of equations solved for each particle (>= ndim) [" + str(lrs) + "]: ") or lrs)

lrs_name = []
for i in range(lrs):
   if i < ndim:
      if i == 0:
         dum_str = "X"
      elif i == 1:
         dum_str = "Y"
      elif i == 2:
         dum_str = "Z"
   else:
      dum_str = "P" + str(i+1)

   dum_str = raw_input("Equation " + str(i+1) + " name [" + dum_str + "]: ") or dum_str
   lrs_name.append(dum_str.upper())

lrp = 1
lrp = int(raw_input("Number of properties for each particle [" + str(lrp) + "]: ") or lrp)

lrp_name = []
for i in range(lrp):
   dum_str = "P" + str(i+1)
   dum_str = raw_input("Property " + str(i+1) + " name [" + dum_str + "]: ") or dum_str
   lrp_name.append(dum_str.upper())

yes = {'yes','y', 'ye', ''}
no = {'no','n'}

lrp_int = 0
lrp_int_name = []
lrp_pro = 0
lex = 2
ley = 2
lez = 2
lee = 1
choiceo = raw_input("Overlap external mesh? (y/n) [n]: ").lower() or "n"
if choiceo in yes:

   lex = int(raw_input("Number of overlap points in x [" + str(lex) + "]: ") or lex)
   if ndim > 1:
      ley = int(raw_input("Number of overlap points in y [" + str(ley) + "]: ") or ley)
   if ndim > 2:
      lez = int(raw_input("Number of overlap points in z [" + str(lez) + "]: ") or lez)
   lee = int(raw_input("Max number of overlap cells per rank [" + str(lee) + "]: ") or lee)

   lrp_int = int(raw_input("Number of interpolated fields [" + str(lrp_int) + "]: ") or lrp_int)
   lrp_int = max(lrp_int,0)
   if lrp_int > 0:
      for i in range(lrp_int):
         dum_str = raw_input("Interpolated property " + str(i+1) + " name (from previous properties): ") 
         lrp_int_name.append(dum_str.upper())
         if lrp_int_name[i] not in lrp_name:
            sys.exit("ERROR: Interpolated property name must be from list of previous properties");


   lrp_pro = int(raw_input("Number of projected fields [" + str(lrp_pro) + "]: ") or lrp_pro)
   lrp_pro = max(lrp_pro,0)


choices = 'no'
lrp_pro_name = []
if lrp_pro > 0:
   for i in range(lrp_pro):
      dum_str = "P" + str(i+1)
      dum_str = raw_input("Projected property " + str(i+1) + " name [" + dum_str + "]: ") or dum_str
      lrp_pro_name.append(dum_str.upper())

# Create Header File
header_file = "PPICLF_USER.h"
with open (header_file, "w") as f:
      if lpartin != lpart:
         f.write ("#define PPICLF_LPART " + str(lpart) + "\n");
      f.write ("#define PPICLF_LRS " + str(lrs) + "\n");
      f.write ("#define PPICLF_LRP " + str(lrp) + "\n");
      if choiceo in yes:
         f.write ("#define PPICLF_LEE " + str(lee) + "\n");
         f.write ("#define PPICLF_LEX " + str(lex) + "\n");
	 if ndim > 1:
            f.write ("#define PPICLF_LEY " + str(ley) + "\n");
	 if ndim > 2:
            f.write ("#define PPICLF_LEZ " + str(lez) + "\n");
      if lrp_int != 0:
         f.write ("#define PPICLF_LRP_INT " + str(lrp_int) + "\n");
      if lrp_pro != 0:
         f.write ("#define PPICLF_LRP_PRO " + str(lrp_pro) + "\n");

      f.write("\n")

      for i in range(lrs):
         f.write("#define PPICLF_J" + lrs_name[i] + " " + str(i+1) + "\n")
      for i in range(lrp):
         f.write("#define PPICLF_R_J" + lrp_name[i] + " " + str(i+1) + "\n")
      if lrp_pro != 0:
         for i in range(lrp_pro):
            f.write("#define PPICLF_P_J" + lrp_pro_name[i] + " " + str(i+1) + "\n")

# Create Fortran File
fortran_file = "ppiclf_user.f"
with open (fortran_file, "w") as f:
      f.write ("! ----------------------------------------------------------------------\n");
      f.write ("! \n");
      f.write ("! Auto-generated ppiclf_user.f file with usergen.py at:\n");
      f.write ("!    " + str(datetime.datetime.now()) + "\n");
      f.write ("! \n");
      f.write ("! Note: REPLACE ANY ? VARIABLES IN THIS FILE WITH DESIRED PROPERTIES\n");
      f.write ("! \n");
      f.write ("! ----------------------------------------------------------------------\n");
      f.write ("      subroutine ppiclf_user_SetYdot\n");
      f.write ("!\n");
      f.write ("      implicit none\n");
      f.write ("!\n");
      f.write ("#include \"PPICLF.h\"\n");
      f.write ("!\n");
      f.write ("! Internal:\n");
      f.write ("!\n");
      f.write ("      integer*4 i\n");
      f.write ("!\n");

      f.write ("! evaluate ydot\n")
      f.write ("      do i=1,ppiclf_npart\n")
      f.write ("\n")
      f.write ("         ! Operations done here to eventually evaluate ppiclf_ydot\n")
      f.write ("\n")
      for i in range(lrs):
         f.write ("         ppiclf_ydot(PPICLF_J" + lrs_name[i] + ",i) = ?\n")
      f.write ("      enddo\n")
      f.write ("\n")
      f.write ("      return\n")
      f.write ("      end\n")

      f.write ("! ----------------------------------------------------------------------\n");
      f.write ("      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)\n");
      f.write ("!\n");
      f.write ("      implicit none\n");
      f.write ("!\n");
      f.write ("! Input:\n");
      f.write ("!\n");
      f.write ("      real*8 y    (PPICLF_LRS)\n")
      f.write ("      real*8 ydot (PPICLF_LRS)\n")
      f.write ("      real*8 ydotc(PPICLF_LRS)\n")
      f.write ("      real*8 rprop(PPICLF_LRP)\n")
      f.write ("!\n");
      f.write ("! Output:\n");
      f.write ("!\n");
      f.write ("      real*8 map  (PPICLF_LRP_PRO)\n")
      f.write ("!\n");
      for i in range(lrp_pro):
         f.write ("      map(PPICLF_P_J" + lrp_pro_name[i] + ") = ?\n")
      f.write ("\n")
      f.write ("      return\n")
      f.write ("      end\n")

      f.write ("! ----------------------------------------------------------------------\n");
      f.write ("      subroutine ppiclf_user_EvalNearestNeighbor\n");
      f.write ("     >                                         (i,j,yi,rpropi,yj,rpropj)\n");
      f.write ("!\n");
      f.write ("      implicit none\n");
      f.write ("!\n");
      f.write ("#include \"PPICLF.h\"\n");
      f.write ("!\n");
      f.write ("! Input:\n");
      f.write ("!\n");
      f.write ("      integer*4 i\n")
      f.write ("      integer*4 j\n")
      f.write ("      real*8 yi    (PPICLF_LRS)\n")
      f.write ("      real*8 rpropi(PPICLF_LRP)\n")
      f.write ("      real*8 yj    (PPICLF_LRS)\n")
      f.write ("      real*8 rpropj(PPICLF_LRP)\n")
      f.write ("!\n");
      f.write ("\n")
      f.write ("      return\n")
      f.write ("      end\n")
      f.write ("! ----------------------------------------------------------------------\n");

# Create External Fortran File
external_file = "external.f"
with open (external_file, "w") as f:
      f.write ("! ----------------------------------------------------------------------\n");
      f.write ("! \n");
      f.write ("! Auto-generated ppiclf_user.h file with usergen.py at:\n");
      f.write ("!    " + str(datetime.datetime.now()) + "\n");
      f.write ("! \n");
      f.write ("! Note: REPLACE ANY ? VARIABLES IN THIS FILE WITH DESIRED PROPERTIES\n");
      f.write ("! \n");
      f.write ("! ----------------------------------------------------------------------\n");
      f.write ("! \n");

      f.write ("! ----------------------------------------------------------------------\n");
      f.write ("! Initialization routine order from external code:\n");
      f.write ("! ----------------------------------------------------------------------\n");
      f.write ("      ! Set-up ppiclF MPI:\n");
      f.write ("      call ppiclf_comm_InitMPI(comm?\n");
      f.write ("     >                        ,id?\n");
      f.write ("     >                        ,np?\n");
      f.write ("! \n");
      f.write ("      ! User sets initial conditions for system, then calls:\n");
      f.write ("      call ppiclf_solve_InitParticle(imethod?\n");
      f.write ("     >                              ,ndim?\n");
      f.write ("     >                              ,iendian?\n");
      f.write ("     >                              ,npart?\n");
      f.write ("     >                              ,y?\n");
      f.write ("     >                              ,rprop?)\n");
      f.write ("! \n");

      if lrp_pro > 0:
         f.write ("      ! User sets up any filter parameters (i.e., Box or Gaussian):\n");
         f.write ("      call ppiclf_solve_InitGaussianFilter(f?,a?,iflg?)\n");
         f.write ("! \n");
      if choiceo in yes:
         f.write ("      ! User sets up overlap mesh:\n");
         f.write ("      call ppiclf_comm_InitOverlapMesh(ncell?\n");
         f.write ("     >                                ,lx?\n");
         f.write ("     >                                ,ly?\n");
         f.write ("     >                                ,lz?\n");
         f.write ("     >                                ,xgrid?\n");
         f.write ("     >                                ,ygrid?\n");
         f.write ("     >                                ,zgrid?)\n");
         f.write ("! \n");
      f.write ("      ! User may specify distance to do a nearest neighbor search:\n");
      f.write ("      ! call ppiclf_solve_InitNeighborBin(W?)\n");
      f.write ("! \n");
      f.write ("      ! User may add boundaries or periodicity below:\n");
      f.write ("! \n");

      f.write ("! ----------------------------------------------------------------------\n");
      f.write ("! Solve routine order from external code:\n");
      f.write ("! ----------------------------------------------------------------------\n");
      f.write ("#include \"PPICLF.h\"\n");

      if lrp_int != 0:
         f.write ("! Interpolate Fields\n")
	 for i in range(lrp_int):
            f.write ("      call ppiclf_solve_InterpFieldUser(PPICLF_R_J" + lrp_int_name[i] + "," + lrp_int_name[i].lower() + "_fld?)\n")
         f.write ("! \n");

      f.write ("! Integrate system\n")
      f.write ("      call ppiclf_solve_IntegrateParticle(istep?\n");
      f.write ("     >                                   ,iostep?\n");
      f.write ("     >                                   ,dt?\n");
      f.write ("     >                                   ,time?)\n");
      f.write ("! ----------------------------------------------------------------------\n");

