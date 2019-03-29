ENAME            = test.out
FC               = mpif77
FFLAGS          += -fdefault-real-8                  \
                   -fdefault-double-8                \
	           -cpp                              \
          	   -fbacktrace                       \
 	           -Wall                             \
		   -z muldefs
#,--allow-multiple-definition 
PPICLF_LOCATION  = $(HOME)/libraries/ppiclF

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF=$(PPICLF_LOCATION)/source
PPICLF_IFLAGS+=-I$(SOURCE_ROOT_PPICLF)
USR_LFLAGS+=-L$(SOURCE_ROOT_PPICLF) -lppiclF

SRC = test.f
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
	$(FC) $(FFLAGS) -c $(SRC) $(PPICLF_IFLAGS) 
	@echo "                           "
	@echo "***************************"
	@echo "*** USER OBJECT SUCCESS ***"
	@echo "***************************"
	@echo "                           "

# Cleaning stuff
clean:
	@rm $(ENAME) *.o
