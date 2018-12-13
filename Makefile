FC = mpif77
FFLAGS = 
ENAME = test.out

GSLIB_PREFIX=gslib_
GSLIB_FPREFIX=fgslib_
GSLIB_IFLAGS=-DPREFIX=$(GSLIB_PREFIX) -DFPREFIX=$(GSLIB_FPREFIX) -DGLOBAL_LONG_LONG

SOURCE_ROOT_GSLIB=3rd_party
GSLIB_IFLAGS+=-I$(SOURCE_ROOT_GSLIB)/include
USR_LFLAGS+=-L$(SOURCE_ROOT_GSLIB)/lib -lgs

FFLAGS+=-fdefault-real-8 -fdefault-double-8 -cpp

#src = test.f lib/lpm_comm.f lib/lpm_io.f lib/lpm_solve.f
#obj = test.o lib/lpm_comm.o lib/lpm_io.o lib/lpm_solve.o

src = test.f 
obj = test.o 

default: getObjs linkObjs 

linkObjs:  $(obj)
	$(FC) $(FFLAGS) -o $(ENAME) $(obj) $(USR_LFLAGS)
	@echo LINK SUCCESS 

getObjs: $(src)
	$(FC) $(FFLAGS) -c $(src) $(GSLIB_IFLAGS)
	@echo OBJECT SUCCESS

clean:
	rm $(ENAME) *.o
