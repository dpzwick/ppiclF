import time
import datetime

# Get User Input
lpart = 100000
lpart = int(raw_input("Max number of particles per rank [" + str(lpart) + "]: ") or lpart)

ndim = 1
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
lrp_pro = 0
lex = 2
ley = 2
lez = 2
lee = 1
choiceo = raw_input("Overlap external mesh? (y/n) [n]: ").lower()
if choiceo in yes:

   lex = int(raw_input("Number of overlap points in x [" + str(lex) + "]: ") or lex)
   if ndim > 1:
      ley = int(raw_input("Number of overlap points in y [" + str(ley) + "]: ") or ley)
   if ndim > 2:
      lez = int(raw_input("Number of overlap points in z [" + str(lez) + "]: ") or lez)
   lee = int(raw_input("Max number of overlap cells per rank [" + str(lee) + "]: ") or lee)

   lrp_int = int(raw_input("Number of interpolated fields [" + str(lrp_int) + "]: ") or lrp_int)
   lrp_int = max(lrp_int,0)

   lrp_pro = int(raw_input("Number of projected fields [" + str(lrp_pro) + "]: ") or lrp_pro)
   lrp_pro = max(lrp_pro,0)


choices = 'no'
bx1 = 1
by1 = 1
bz1 = 1
lrp_pro_name = []
if lrp_pro > 0:

   if lrp_pro > 0:
      for i in range(lrp_pro):
         dum_str = "P" + str(i+1)
         dum_str = raw_input("Projected property " + str(i+1) + " name [" + dum_str + "]: ") or dum_str
         lrp_pro_name.append(dum_str.upper())

   choices = raw_input("Project to SubBins? (y/n) [n]: ").lower()

   if choices in yes:
      bx1 = int(raw_input("Max size of SubBin grid in x [" + str(bx1) + "]: ") or bx1)
      if ndim > 1:
         by1 = int(raw_input("Max size of SubBin grid in y [" + str(by1) + "]: ") or by1)
      if ndim > 2:
         bz1 = int(raw_input("Max size of SubBin grid in z [" + str(bz1) + "]: ") or bz1)

# Create Header File
header_file = "ppiclf_user.h"
with open (header_file, "w") as f:
      f.write ("C Auto-generated ppiclf_user.h file with usergen.py at:\n");
      f.write ("C    " + str(datetime.datetime.now()) + "\n");
      f.write ("C \n");

      f.write ("C Max number of particles per rank \n");
      f.write ("#define PPICLF_LPART " + str(lpart) + "\n");
      f.write("\n")

      f.write ("C Number of equations solved for each particle \n");
      f.write ("#define PPICLF_LRS " + str(lrs) + "\n");
      f.write("\n")

      f.write ("C Naming of equations solved for each particle \n");
      for i in range(lrs):
         f.write("#define PPICLF_J" + lrs_name[i] + " " + str(i+1) + "\n")
      f.write("\n")

      f.write ("C Number of additional properties for each particle \n");
      f.write ("#define PPICLF_LRP " + str(lrp) + "\n");
      f.write("\n")

      f.write ("C Naming of additional properties for each particle \n");
      for i in range(lrp):
         f.write("#define PPICLF_R_J" + lrp_name[i] + " " + str(i+1) + "\n")
      f.write("\n")

      if choiceo in yes:
         f.write ("C Size of external overlap grid \n");
         f.write ("#define PPICLF_LEX " + str(lex) + "\n");
	 if ndim > 1:
            f.write ("#define PPICLF_LEY " + str(ley) + "\n");
	 if ndim > 2:
            f.write ("#define PPICLF_LEZ " + str(lez) + "\n");
         f.write ("#define PPICLF_LEE " + str(lee) + "\n");
         f.write("\n")

      if lrp_int != 0:
         f.write ("C Number of interpolated fields \n");
         f.write ("#define PPICLF_LRP_INT " + str(lrp_int) + "\n");
         f.write("\n")

      if lrp_pro != 0:
         f.write ("C Number of projected fields \n");
         f.write ("#define PPICLF_LRP_PRO " + str(lrp_pro) + "\n");
         f.write("\n")

         f.write ("C Naming of projected fields \n");
         for i in range(lrp_pro):
            f.write("#define PPICLF_P_J" + lrp_pro_name[i] + " " + str(i+1) + "\n")
         f.write("\n")

      if choices in yes:
         f.write ("C Max size of SubBin grid \n");
         f.write ("#define PPICLF_BX1 " + str(bx1) + "\n");
	 if ndim > 1:
            f.write ("#define PPICLF_BY1 " + str(by1) + "\n");
	 if ndim > 2:
            f.write ("#define PPICLF_BZ1 " + str(bz1) + "\n");
         f.write("\n")

# Create Fortran File
fortran_file = "ppiclf_user.f"
with open (fortran_file, "w") as f:
      f.write ("C ----------------------------------------------------------------------\n");
      f.write ("C \n");
      f.write ("C Auto-generated ppiclf_user.f file with usergen.py at:\n");
      f.write ("C    " + str(datetime.datetime.now()) + "\n");
      f.write ("C \n");
      f.write ("C Note: REPLACE ANY ? VARIABLES IN THIS FILE WITH DESIRED PROPERTIES\n");
      f.write ("C \n");
      f.write ("C ----------------------------------------------------------------------\n");
      f.write ("      subroutine ppiclf_user_SetYdot(time,y,ydot)\n");
      f.write ("#include \"PPICLF\"\n");
      f.write ("\n")
      f.write ("      real time\n")
      f.write ("      real y(*)\n")
      f.write ("      real ydot(*)\n")
      f.write ("\n")
      if lrp_int != 0:
         f.write ("C Interpolate Fields\n")
         f.write ("      call ppiclf_solve_InitInterp\n")
	 for i in range(lrp_int):
            f.write ("         call ppiclf_solve_InterpField(PPICLF_R_J?,field?" + str(i+1) +")\n")
         f.write ("      call ppiclf_solve_FinalizeInterp\n")
         f.write ("\n")

      f.write ("C Evaluate Ydot\n")
      f.write ("      do i=1,ppiclf_npart\n")
      f.write ("         ! striding in y and ydot\n")
      f.write ("         j = PPICLF_LRS*(i-1)\n")
      f.write ("\n")
      f.write ("         ! User may perform a nearest neighbor search to \n")
      f.write ("         ! activate EvalNearestNeighbor():\n")
      f.write ("         ! call ppiclf_solve_NearestNeighbor(i)\n")
      f.write ("\n")
      f.write ("         ! Operations done here to eventually evaluate ydot\n")
      for i in range(lrs):
         f.write ("         ydot(PPICLF_J" + lrs_name[i] + " +j) = ?\n")
      f.write ("      enddo\n")
      f.write ("\n")
      f.write ("      return\n")
      f.write ("      end\n")

      f.write ("C ----------------------------------------------------------------------\n");
      f.write ("      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)\n");
      f.write ("\n")
      f.write ("      real map(*)\n")
      f.write ("      real y(*)\n")
      f.write ("      real ydot(*)\n")
      f.write ("      real ydotc(*)\n")
      f.write ("      real rprop(*)\n")
      f.write ("\n")
      f.write ("C Set which propreties get projected\n")
      for i in range(lrp_pro):
         f.write ("      map(PPICLF_P_J" + lrp_pro_name[i] + ") = ?\n")
      f.write ("\n")
      f.write ("      return\n")
      f.write ("      end\n")

      f.write ("C ----------------------------------------------------------------------\n");
      f.write ("      subroutine ppiclf_user_EvalNearestNeighbor\n");
      f.write ("     >                                         (i,j,yi,rpropi,yj,rpropj)\n");
      f.write ("#include \"PPICLF\"\n");
      f.write ("\n")
      f.write ("      integer i\n")
      f.write ("      integer j\n")
      f.write ("      real yi(*)     ! PPICLF_LRS\n")
      f.write ("      real rpropi(*) ! PPICLF_LRP\n")
      f.write ("      real yj(*)     ! PPICLF_LRS\n")
      f.write ("      real rpropj(*) ! PPICLF_LRP\n")
      f.write ("\n")
      f.write ("C - Nearest neighbor search is performed for particle i with solution\n")
      f.write ("C   vector yi and property vector rpropi\n")
      f.write ("C - Neighbor j particle has solution vector yj and proprty vector rpropj\n")
      f.write ("C - When j == 0, the j particle is not a particle but a wall with the\n")
      f.write ("C   first ndim values of yj being the nearest point on the wall\n")
      f.write ("C - If information needs to be computed for particle i, it can safely\n")
      f.write ("C   be stored/summed in ppiclf_ydotc(*,i) and evaluated in SetYdot()\n")
      f.write ("C   since ppiclf_ydotc(*,*) is cleared before each call to SetYdot()\n")
      f.write ("\n")
      f.write ("      return\n")
      f.write ("      end\n")
      f.write ("C ----------------------------------------------------------------------\n");

# Create External Fortran File
external_file = "external.f"
with open (external_file, "w") as f:
      f.write ("C ----------------------------------------------------------------------\n");
      f.write ("C \n");
      f.write ("C Auto-generated ppiclf_user.h file with usergen.py at:\n");
      f.write ("C    " + str(datetime.datetime.now()) + "\n");
      f.write ("C \n");
      f.write ("C Note: REPLACE ANY ? VARIABLES IN THIS FILE WITH DESIRED PROPERTIES\n");
      f.write ("C \n");
      f.write ("C ----------------------------------------------------------------------\n");

      f.write ("C Initialization routine order from external code:\n");
      f.write ("C     ! Set-up ppiclF MPI:\n");
      f.write ("C     call ppiclf_comm_InitMPI(comm?,nid?,np?)\n");
      f.write ("C \n");
      f.write ("C     ! User sets initial conditions for particle, then calls:\n");
      f.write ("C     call ppiclf_solve_InitParticle(method?,ndim?,endian?,npart?,y?)\n");
      f.write ("C \n");

      if lrp_pro > 0:
         f.write ("C     ! User sets up any filter parameters (i.e., Box or Gaussian):\n");
         f.write ("C     call ppiclf_solve_InitGaussianFilter(w?,a?,n?)\n");
         f.write ("C     call ppiclf_solve_InitBoxFilter(w?,n?)\n");
         f.write ("C \n");
      if choiceo in yes:
         f.write ("C     ! User sets up overlap mesh:\n");
         f.write ("C     call ppiclf_commm_InitOverlapMesh(ne?,lx?,ly?,lz?,x?,y?,z?)\n");
         f.write ("C \n");
      f.write ("C     ! User may specify distance to do a nearest neighbor search:\n");
      f.write ("C     call ppiclf_solve_InitNeighborBin(w?)\n");
      f.write ("C \n");
      f.write ("C     ! User may add walls or periodicity below:\n");
      f.write ("C \n");

      f.write ("C ----------------------------------------------------------------------\n");
      f.write ("C Call each time step from external code:\n");
      f.write ("C #include \"PPICLF\"\n");
      f.write ("C     call ppiclf_solve_IntegrateParticle(istep?,iostep?,dt?,time?\n");
      f.write ("C    >                                   ,ppiclf_y,ppiclf_ydot)\n");
      f.write ("C ----------------------------------------------------------------------\n");

