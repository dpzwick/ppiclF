FC               = mpif77 # Valid MPI Fortran compiler
CC               = mpicc  # Valid MPI C compiler
FFLAGS          += -cpp               \
             	   -fbacktrace        \
#	           -Wall              \
#                  -fdefault-real-8   
#FFLAGS          += -DPPICLC
 	           
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
SOURCE_ROOT_GSLIB_OBJ=$(SOURCE_ROOT_GSLIB)/gslib/src
GSLIB_IFLAGS+=-I$(SOURCE_ROOT_GSLIB)/include

SRC = $(SOURCE_ROOT_PPICLF)/ppiclf_user.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_op.f       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.f       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.f    \
      $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.f      
OBJ = $(SOURCE_ROOT_PPICLF)/ppiclf_user.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_op.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.o    \
      $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.o

# Make commands
default: makeThird getObjs libObjs 

libObjs: $(OBJ)
	@ar crv $(SOURCE_ROOT_PPICLF)/libppiclF.a $(OBJ) $(SOURCE_ROOT_GSLIB_OBJ)/*.o
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
	./install $(CC) $(FC)
	@echo "                       "
	@echo "***********************"
	@echo "*** INSTALLED GSLIB ***"
	@echo "***********************"
	@echo "                       "

clean:
	 rm -r $(SOURCE_ROOT_PPICLF)/*.o                        
	 rm -r $(SOURCE_ROOT_PPICLF)/libppiclF.a                
	 rm -r $(SOURCE_ROOT_GSLIB)/gslib                       
	 rm -r $(SOURCE_ROOT_GSLIB)/lib                         
	 rm -r $(SOURCE_ROOT_GSLIB)/include                     
	 rm -r $(SOURCE_ROOT_GSLIB)/*.tar.gz                    
	 rm -r $(INSTALL_LOCATION)/short_tests/*.pyc            
	 rm -r $(INSTALL_LOCATION)/short_tests/lib/*.pyc        
	 rm -r $(INSTALL_LOCATION)/short_tests/test1/test1.log*
	 rm -r $(INSTALL_LOCATION)/short_tests/test1/*.vtu
	 cd $(INSTALL_LOCATION)/short_tests/test1/; \
	make clean

