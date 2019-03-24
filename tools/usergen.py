# Get User Input
lpart = 100000
lpart = int(raw_input("Max number of particles per rank [" + str(lpart) + "]: ") or lpart)

ndim = 1
ndim = int(raw_input("Problem dimensions [" + str(ndim) + "]: ") or ndim)

lrs = 1
lrs = int(raw_input("Number of equations solved for each particle [" + str(lrs) + "]: ") or lrs)

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
      dum_str = str(i+1)

   dum_str = raw_input("Equation " + str(i+1) + " name [" + dum_str + "]: ") or dum_str
   lrs_name.append(dum_str.upper())

lrp = 1
lrp = int(raw_input("Number of properties for each particle [" + str(lrp) + "]: ") or lrp)

lrp_name = []
for i in range(lrp):
   dum_str = "P" + str(i+1)
   dum_str = raw_input("Property " + str(i+1) + " name [" + dum_str + "]: ") or dum_str
   lrp_name.append(dum_str.upper())

# Create Header File
header_file = "ppiclf_user.h"
with open (header_file, "w") as f:
      f.write ("C Max number of particles per rank \n");
      f.write ("#define PPICLF_LPART " + str(lpart) + "\n\n");

      f.write ("C Number of equations solved for each particle \n");
      f.write ("#define PPICLF_LRS " + str(lrs) + "\n\n");

      f.write ("C Naming of equations solved for each particle \n");
      for i in range(lrs):
         f.write("#define PPICLF_J" + lrs_name[i] + "\n")
