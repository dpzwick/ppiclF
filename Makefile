FC               = mpif77
FFLAGS          += -fdefault-real-8   \
                   -fdefault-double-8 \
	           -cpp               \
          	   -fbacktrace        \
	           -Wall
INSTALL_LOCATION = .

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF=$(INSTALL_LOCATION)/source
PPICLF_IFLAGS+=-I$(SOURCE_ROOT_PPICLF)

# GSLIB LIBRARY
GSLIB_IFLAGS = -DPREFIX=gslib_   \
               -DFPREFIX=fgslib_ \
	       -DGLOBAL_LONG_LONG
SOURCE_ROOT_GSLIB=$(INSTALL_LOCATION)/3rd_party/gslib
GSLIB_IFLAGS+=-I$(SOURCE_ROOT_GSLIB)/include

SRC = $(SOURCE_ROOT_PPICLF)/ppiclf_user.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_op.f \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.f       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.f
OBJ = $(SOURCE_ROOT_PPICLF)/ppiclf_user.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_op.o \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.o

# Make commands
default: makeThird getObjs libObjs 

libObjs: $(OBJ)
	@ar crv $(SOURCE_ROOT_PPICLF)/libppiclF.a $(OBJ) 
	@echo "                       "
	@echo "***********************"
	@echo "*** LIBRARY SUCCESS ***"
	@echo "***********************"
	@echo "                       "

getObjs: $(SRC)
	$(FC) $(FFLAGS) -c $(SRC) $(GSLIB_IFLAGS) $(PPICLF_IFLAGS)
	mv *.o $(SOURCE_ROOT_PPICLF)
	@echo "                              "
	@echo "******************************"
	@echo "*** LIBRARY OBJECT SUCCESS ***"
	@echo "******************************"
	@echo "                              "

makeThird: $(SOURCE_ROOT_GSLIB)/install
	cd $(SOURCE_ROOT_GSLIB); \
	./install
	@echo "                       "
	@echo "***********************"
	@echo "*** INSTALLED GSLIB ***"
	@echo "***********************"
	@echo "                       "

clean:
	@rm -r $(SOURCE_ROOT_PPICLF)/*.o         \
	       $(SOURCE_ROOT_PPICLF)/libppiclF.a \
	       $(SOURCE_ROOT_GSLIB)/gslib        \
	       $(SOURCE_ROOT_GSLIB)/lib          \
	       $(SOURCE_ROOT_GSLIB)/include      \
	       $(SOURCE_ROOT_GSLIB)/*.tar.gz
