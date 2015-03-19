#modify the following options:
#set FC to be the fortran compiler (mpi wrapper generally)
FC=mpif90

#set the compilation options for the fortran compiler
FCOPTS = -O3 -fopenmp
FCOPTS = -fimplicit-none -Wall -fbounds-check -fbacktrace -fopenmp

#set the location of the libraries
LIB_LIST = -L$(METISDIR) -l$(METISLIB) -lblas -llapack

#set the option for placing module files
MODFLAG=-J

#JULESDIR=$(PWD)
#This line for doing the build from NetBeans
JULESDIR=/home/xw904346/SatWin/jules/JULES_trunk/
#set the path for all further actions
DA_DIR = $(JULESDIR)/src/control/standalone/da_filter
RUNDIR = /home/xw904346/SatWin/jules/harvard/


all: EMPIRE DO_COPY #ALLTESTS TEST_R TEST_Q TEST_HQHTR #TEST_H

COPY_FILES=$(DA_DIR)/controllers/pf_parameters.dat
DO_COPY:
	cp -f $(COPY_FILES) $(RUNDIR)

OBJ_DIR=$(DA_DIR)/obj/
BIN_DIR=$(DA_DIR)/bin/
MODLOC:=$(OBJ_DIR)
SR_FILTS=$(DA_DIR)/filters/
SR_UTILS=$(DA_DIR)/utils/
SR_CONTS=$(DA_DIR)/controllers/
SR_DATAS=$(DA_DIR)/data/
SR_TESTS=$(DA_DIR)/tests/
SR_OPERS=$(DA_DIR)/operations/
#OBJ_LIST= sizes.o pf_couple.o Qdata.o Rdata.o equivalent_weights_filter.o comms.o gen_rand.o random_d.o proposal_filter.o histogram.o pf_control.o data_io.o model_specific.o operator_wrappers.o quicksort.o resample.o diagnostics.o perturb_particle.o genQ.o sir_filter.o stochastic_model.o tests.o letkf_analysis.o deterministic_model.o
OBJ_LIST= 			    \
       Qdata.o				    \
       Rdata.o				    \
       empire_mod.o				    \
       equivalent_weights_step.o			    \
       gen_rand.o				    \
       random_d.o				    \
       proposal_filter.o				    \
       histogram.o				    \
       pf_control.o				    \
       data_io.o				    \
       model_specific.o				    \
       operator_wrappers.o			    \
       quicksort.o resample.o			    \
       diagnostics.o				    \
       perturb_particle.o			    \
       genQ.o				    \
       sir_filter.o				    \
       stochastic_model.o			    \
       tests.o				    \
       letkf_analysis.o				    \
       deterministic_model.o			    \
       empire.o

#this is the list the exe is dependent upon: obj/*.o
OBJECTS=$(addprefix $(OBJ_DIR),$(OBJ_LIST))

#set the module file location to obj/
FCOPTS+=$(MODFLAG) $(MODLOC)

#$(OBJ_DIR)sizes.o: $(SR_CONTS)sizes.f90
#	$(FC) $(FCOPTS) -c $(SR_CONTS)sizes.f90 -o $@
	
$(OBJ_DIR)histogram.o: $(SR_UTILS)histogram.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)histogram.f90 -o $@
	@echo histogram

$(OBJ_DIR)stochastic_model.o: $(SR_FILTS)stochastic_model.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)stochastic_model.f90 -o $@

$(OBJ_DIR)deterministic_model.o: $(SR_FILTS)deterministic_model.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)deterministic_model.f90 -o $@

$(OBJ_DIR)sir_filter.o: $(SR_FILTS)sir_filter.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)sir_filter.f90 -o $@

$(OBJ_DIR)genQ.o: $(SR_UTILS)genQ.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)genQ.f90 -o $@

$(OBJ_DIR)Qdata.o: $(SR_DATAS)Qdata.f90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_DATAS)Qdata.f90 -o $@

$(OBJ_DIR)Rdata.o: $(SR_DATAS)Rdata.f90
	$(FC) $(FCOPTS) -c $(SR_DATAS)Rdata.f90 -o $@

$(OBJ_DIR)model_specific.o: model_specific.f90 $(OBJ_DIR)pf_control.o 
	$(FC) $(FCOPTS) -c model_specific.f90 -o $@
	@echo model_specific

$(OBJ_DIR)operator_wrappers.o: $(SR_OPERS)operator_wrappers.f90	\
			   $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)operator_wrappers.f90 -o $@

$(OBJ_DIR)pf_control.o: $(SR_CONTS)pf_control.f90  $(OBJ_DIR)histogram.o
	$(FC) $(FCOPTS) -c $(SR_CONTS)pf_control.f90 -o $@

$(OBJ_DIR)perturb_particle.o: $(SR_OPERS)perturb_particle.f90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)perturb_particle.f90 -o $@

$(OBJ_DIR)data_io.o: $(SR_UTILS)data_io.f90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_UTILS)data_io.f90 -o $@

$(OBJ_DIR)proposal_filter.o: $(SR_FILTS)proposal_filter.f90 $(OBJ_DIR)random_d.o
	$(FC) $(FCOPTS) -c $(SR_FILTS)proposal_filter.f90 -o $@
	@echo proposal_filter

#$(OBJ_DIR)comms.o: $(SR_UTILS)comms.f90 $(OBJ_DIR)pf_control.o
#	$(FC) $(FCOPTS) -c $(SR_UTILS)comms.f90 -o $@

$(OBJ_DIR)equivalent_weights_step.o: $(SR_FILTS)equivalent_weights_step.f90 \
				   $(OBJ_DIR)random_d.o
	$(FC) $(FCOPTS) -c $(SR_FILTS)equivalent_weights_step.f90 -o $@

$(OBJ_DIR)random_d.o: $(SR_UTILS)random_d.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)random_d.f90 -o $@

$(OBJ_DIR)resample.o: $(SR_OPERS)resample.f90 $(OBJ_DIR)pf_control.o $(OBJ_DIR)random_d.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)resample.f90 -o $@

$(OBJ_DIR)diagnostics.o: $(SR_UTILS)diagnostics.f90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_UTILS)diagnostics.f90 -o $@

$(OBJ_DIR)quicksort.o: $(SR_UTILS)quicksort.f90
	$(FC) $(FCOPTS) -c $(SR_UTILS)quicksort.f90 -o $@

#$(OBJ_DIR)pf_couple.o: $(SR_CONTS)pf_couple.f90 $(OBJ_DIR)comms.o $(OBJ_DIR)pf_control.o
#	$(FC) $(FCOPTS) -c $(SR_CONTS)pf_couple.f90 -o $@

#$(OBJ_DIR)empi_wrappers.o: empi_wrappers.F90 
#	$(FC) $(FCOPTS) -c empi_wrappers.F90 -o $@
	
$(OBJ_DIR)empire_mod.o: empire_mod.F90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c empire_mod.F90 -o $@
	
$(OBJ_DIR)empire.o: empire.F90 $(OBJ_DIR)empire_mod.o $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c empire.F90 -o $@

$(OBJ_DIR)alltests.o: $(SR_TESTS)alltests.f90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)alltests.f90 -o $@

$(OBJ_DIR)test_h.o: $(SR_TESTS)test_h.f90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)test_h.f90 -o $@

$(OBJ_DIR)test_r.o: $(SR_TESTS)test_r.f90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)test_r.f90 -o $@

$(OBJ_DIR)tests.o: $(SR_TESTS)tests.f90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)tests.f90 -o $@

$(OBJ_DIR)test_q.o: $(SR_TESTS)test_q.f90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)test_q.f90 -o $@

$(OBJ_DIR)test_hqhtr.o: $(SR_TESTS)test_hqhtr.f90 $(OBJ_DIR)pf_control.o
	$(FC) $(FCOPTS) -c $(SR_TESTS)test_hqhtr.f90 -o $@

$(OBJ_DIR)gen_rand.o: $(SR_OPERS)gen_rand.f90 $(OBJ_DIR)random_d.o
	$(FC) $(FCOPTS) -c $(SR_OPERS)gen_rand.f90 -o $@


#ENKF SECTION:
$(OBJ_DIR)letkf_analysis.o: $(SR_FILTS)letkf_analysis.f90
	$(FC) $(FCOPTS) -c $(SR_FILTS)letkf_analysis.f90 -o $@




EMPIRE: $(OBJECTS) 
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(JULESDIR)empire.exe $(OBJECTS) $(LIB_LIST)


OBJS_ALLTEST = $(shell echo $(OBJECTS) | sed 's/empire/alltests/g') 
ALLTESTS: $(OBJS_ALLTEST)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN_DIR)alltests $(OBJS_ALLTEST) $(LIB_LIST)

OBJS_TEST_H = $(shell echo $(OBJECTS) | sed 's/empire/test_h/g') 
TEST_H: $(OBJS_TEST_H)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN_DIR)test_h $(OBJS_TEST_H) $(LIB_LIST)

OBJS_TEST_R = $(shell echo $(OBJECTS) | sed 's/empire/test_r/g') 
TEST_R: $(OBJS_TEST_R)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN_DIR)test_r $(OBJS_TEST_R) $(LIB_LIST)

OBJS_TEST_Q = $(shell echo $(OBJECTS) | sed 's/empire/test_q/g') 
TEST_Q: $(OBJS_TEST_Q)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN_DIR)test_q $(OBJS_TEST_Q) $(LIB_LIST)

OBJS_TEST_HQHTR = $(shell echo $(OBJECTS) | sed 's/empire/test_hqhtr/g') 
TEST_HQHTR: $(OBJS_TEST_HQHTR)
	$(FC) $(FCOPTS) $(LOADOPTS) -o $(BIN_DIR)test_hqhtr $(OBJS_TEST_HQHTR) $(LIB_LIST)


docs: doc_latex

doc_html: FORCE
	doxygen .Doxyfile

doc_latex: doc_html
	cd doc/latex && make

FORCE:


clean:
	rm -f obj/* bin/*