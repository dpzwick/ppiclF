ENAME            = test.out
CC               = mpiCC
FC               = mpif77
PPICLF_LOCATION  = $(HOME)/ppiclF
CFLAGS           = 
FFLAGS           = 

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF=$(PPICLF_LOCATION)/source
PPICLF_IFLAGS+=-I$(SOURCE_ROOT_PPICLF)
USR_LFLAGS+=-L$(SOURCE_ROOT_PPICLF) -lppiclF

# STD C++ LIBRARY SINCE WE LINK WITH $(FC)
USR_LFLAGS+=-lstdc++

SRC = test.cpp
OBJ = test.o

# Make commands
default: getObjs linkObjs

linkObjs:  $(OBJ)
	$(FC) $(FFLAGS) -o $(ENAME) $(OBJ) $(USR_LFLAGS) 
	@echo "                         "
	@echo "*************************"
	@echo "*** USER LINK SUCCESS ***"
	@echo "*************************"
	@echo "                         "

getObjs: $(SRC)
	$(CC) $(CFLAGS) -c $(SRC) $(PPICLF_IFLAGS)
	@echo "                           "
	@echo "***************************"
	@echo "*** USER OBJECT SUCCESS ***"
	@echo "***************************"
	@echo "                           "

# Cleaning stuff
clean:
	@rm $(ENAME) *.o
