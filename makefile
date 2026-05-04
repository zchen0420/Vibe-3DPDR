# Makefile for 3DPDR
# Written by T. G. Bisbas
# University College London

F90               = gfortran
CC                = gcc
CPPFLAGS          = -cpp
SUNDIALS_PREFIX  ?= .bin/sundials
INCLUDES          = -I$(SUNDIALS_PREFIX)/include
LIBRARIES         = -L$(SUNDIALS_PREFIX)/lib
LIBS              = -lsundials_cvode -lsundials_nvecserial -lm
OPTIMISE          = 0

BUILD_DIR         = build
OBJ_DIR           = $(BUILD_DIR)/obj
MOD_DIR           = $(BUILD_DIR)/mod
TEST_DIR          = $(BUILD_DIR)/tests
LOG_DIR           = $(BUILD_DIR)/logs
CHECK_LOG         = $(LOG_DIR)/run_3dpdr_check.log

DIMENSIONS        = 1
NETWORK           = REDUCED
DUST              = 2
GUESS_TEMP        = 1
THERMALBALANCE    = 1
TEMP_FIX          = 1
CO_FIX            = 1
H2FORM            = 1

ifeq ($(THERMALBALANCE),1)
  CFLAGS += -DTHERMALBALANCE
endif
ifeq ($(GUESS_TEMP),1)
  CFLAGS += -DGUESS_TEMP
endif
ifeq ($(DIMENSIONS),1)
  CFLAGS += -DPSEUDO_1D
endif
ifeq ($(DIMENSIONS),2)
  CFLAGS += -DPSEUDO_2D
endif
ifeq ($(OPTIMISE),0)
  OPT += -O0
endif
ifeq ($(OPTIMISE),1)
  OPT += -O1
endif
ifeq ($(OPTIMISE),2)
  OPT += -O2
endif
ifeq ($(OPTIMISE),3)
  OPT += -O3
endif
ifeq ($(OPTIMISE),4)
  OPT += -fast
endif
ifeq ($(DUST),1)
  CFLAGS += -DDUST
endif
ifeq ($(DUST),2)
  CFLAGS += -DDUST2
endif
ifeq ($(CO_FIX),1)
  CFLAGS += -DCO_FIX
endif
ifeq ($(H2FORM),1)
  CFLAGS += -DH2FORM
endif
ifeq ($(TEMP_FIX),1)
  CFLAGS += -DTEMP_FIX
endif
ifeq ($(NETWORK),REDUCED)
  CFLAGS += -DREDUCED
endif
ifeq ($(NETWORK),FULL)
  CFLAGS += -DFULL
endif
ifeq ($(NETWORK),MYNETWORK)
  CFLAGS += -DMYNETWORK
endif

FFLAGS            = $(CPPFLAGS) $(OPT) $(CFLAGS) -J$(MOD_DIR) -I$(MOD_DIR)
C_COMPILE_FLAGS   = -O0 $(INCLUDES) $(CFLAGS)

VPATH             = src/app src/io src/init src/evolution src/physics src/physics/chemistry .

MODULE_SRC += definitions.F90 healpix_types.F90 healpix_state.F90 shielding_tables.F90
MODULE_SRC += physical_parameters.F90 photorate_interfaces.F90 chemistry_controls.F90
MODULE_SRC += coolants.F90 modules.F90
MODULE_SRC += runtime_config.F90 grid_io.F90 spatial_index.F90 memory.F90
MODULE_SRC += excitation.F90 convergence.F90 radiation.F90 columns.F90
MODULE_SRC += dark_region.F90 iteration_chemistry.F90 level_population_solver.F90
MODULE_SRC += evolution_setup.F90 thermal_balance.F90 iteration_convergence.F90
MODULE_SRC += output.F90 coolant_io.F90 chemistry_io.F90 initial_conditions.F90
MODULE_SRC += geometry_setup.F90 particle_storage.F90

CODE_F90_SRC += healpix.F90 input_parameters.F90 solvlevpop.F90 read_species.F90
CODE_F90_SRC += read_rates.F90 heapsort.F90 calc_reac_rates.F90 shield.F90
CODE_F90_SRC += spline.F90 escape_probability.F90 eval_points.F90
CODE_F90_SRC += sub_calculate_heating.F90 find_Ccoeff.F90 read_input.F90
ifeq ($(DUST),2)
CODE_F90_SRC += dust_t.F90
endif
ifeq ($(H2FORM),1)
CODE_F90_SRC += h2_form.F90
endif

CHEM_C_SRC += calculate_abundances.c
ifeq ($(NETWORK),REDUCED)
CHEM_C_SRC += odes_reduced.c jacobian_reduced.c
endif
ifeq ($(NETWORK),FULL)
CHEM_C_SRC += odes_full.c jacobian_full.c
endif
ifeq ($(NETWORK),MYNETWORK)
CHEM_C_SRC += odes_mynetwork.c jacobian_mynetwork.c
endif

MODULE_OBJ        = $(patsubst %.F90,$(OBJ_DIR)/%.o,$(MODULE_SRC))
CODE_F90_OBJ      = $(patsubst %.F90,$(OBJ_DIR)/%.o,$(CODE_F90_SRC))
CHEM_C_OBJ        = $(patsubst %.c,$(OBJ_DIR)/%.o,$(CHEM_C_SRC))
OBJ               = $(MODULE_OBJ) $(CODE_F90_OBJ) $(CHEM_C_OBJ)
MAIN_OBJ          = $(OBJ_DIR)/3DPDR.o

TEST_EXECS        = $(TEST_DIR)/test_excitation $(TEST_DIR)/test_coolants
TEST_EXECS       += $(TEST_DIR)/test_convergence $(TEST_DIR)/test_radiation $(TEST_DIR)/test_columns

.PHONY: all run unit-test check clean compress

all: 3DPDR

$(OBJ_DIR) $(MOD_DIR) $(TEST_DIR) $(LOG_DIR):
	mkdir -p $@

$(OBJ_DIR)/definitions.o: definitions.F90 | $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/healpix_types.o: healpix_types.F90 $(OBJ_DIR)/definitions.o | $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/modules.o: modules.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o | $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: %.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/modules.o | $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: %.c | $(OBJ_DIR)
	$(CC) $(C_COMPILE_FLAGS) -c $< -o $@

$(OBJ_DIR)/3DPDR.o: 3DPDR.F90 \
	$(OBJ_DIR)/modules.o \
	$(OBJ_DIR)/memory.o \
	$(OBJ_DIR)/runtime_config.o \
	$(OBJ_DIR)/grid_io.o \
	$(OBJ_DIR)/spatial_index.o \
	$(OBJ_DIR)/coolant_io.o \
	$(OBJ_DIR)/chemistry_io.o \
	$(OBJ_DIR)/initial_conditions.o \
	$(OBJ_DIR)/geometry_setup.o \
	$(OBJ_DIR)/particle_storage.o \
	$(OBJ_DIR)/evolution_setup.o \
	$(OBJ_DIR)/iteration_chemistry.o \
	$(OBJ_DIR)/level_population_solver.o \
	$(OBJ_DIR)/thermal_balance.o \
	$(OBJ_DIR)/iteration_convergence.o \
	$(OBJ_DIR)/output.o \
	| $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -c $< -o $@

3DPDR: $(OBJ) $(MAIN_OBJ)
	$(F90) $(OPT) $(CFLAGS) $(LIBRARIES) -o 3DPDR $(OBJ) $(MAIN_OBJ) $(LIBS)

run: 3DPDR
	./3DPDR

unit-test: $(TEST_EXECS)
	$(TEST_DIR)/test_excitation
	$(TEST_DIR)/test_coolants
	$(TEST_DIR)/test_convergence
	$(TEST_DIR)/test_radiation
	$(TEST_DIR)/test_columns

check: 3DPDR unit-test | $(LOG_DIR)
	./3DPDR > $(CHECK_LOG) 2>&1
	grep -q "RESULT status=converged iterations=227" $(CHECK_LOG)

$(TEST_DIR)/test_excitation: tests/test_excitation.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/excitation.o | $(TEST_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -o $@ tests/test_excitation.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/excitation.o

$(TEST_DIR)/test_coolants: tests/test_coolants.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/coolants.o | $(TEST_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -o $@ tests/test_coolants.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/coolants.o

$(TEST_DIR)/test_convergence: tests/test_convergence.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o $(OBJ_DIR)/excitation.o $(OBJ_DIR)/convergence.o | $(TEST_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -o $@ tests/test_convergence.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o $(OBJ_DIR)/excitation.o $(OBJ_DIR)/convergence.o

$(TEST_DIR)/test_radiation: tests/test_radiation.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/modules.o $(OBJ_DIR)/radiation.o | $(TEST_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -o $@ tests/test_radiation.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/modules.o $(OBJ_DIR)/radiation.o

$(TEST_DIR)/test_columns: tests/test_columns.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/modules.o $(OBJ_DIR)/columns.o | $(TEST_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -o $@ tests/test_columns.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/modules.o $(OBJ_DIR)/columns.o

$(OBJ_DIR)/healpix_state.o: healpix_state.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o
$(OBJ_DIR)/shielding_tables.o: shielding_tables.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o
$(OBJ_DIR)/physical_parameters.o: physical_parameters.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o
$(OBJ_DIR)/photorate_interfaces.o: photorate_interfaces.F90 $(OBJ_DIR)/definitions.o $(OBJ_DIR)/healpix_types.o $(OBJ_DIR)/shielding_tables.o
$(OBJ_DIR)/chemistry_controls.o: chemistry_controls.F90 $(OBJ_DIR)/healpix_types.o
$(OBJ_DIR)/coolants.o: coolants.F90 $(OBJ_DIR)/healpix_types.o
$(OBJ_DIR)/runtime_config.o: runtime_config.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/shielding_tables.o $(OBJ_DIR)/physical_parameters.o
$(OBJ_DIR)/grid_io.o: grid_io.F90 $(OBJ_DIR)/modules.o
$(OBJ_DIR)/spatial_index.o: spatial_index.F90 $(OBJ_DIR)/modules.o
$(OBJ_DIR)/memory.o: memory.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/shielding_tables.o $(OBJ_DIR)/coolants.o
$(OBJ_DIR)/excitation.o: excitation.F90 $(OBJ_DIR)/healpix_types.o
$(OBJ_DIR)/convergence.o: convergence.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o $(OBJ_DIR)/excitation.o
$(OBJ_DIR)/radiation.o: radiation.F90 $(OBJ_DIR)/modules.o
$(OBJ_DIR)/columns.o: columns.F90 $(OBJ_DIR)/modules.o
$(OBJ_DIR)/dark_region.o: dark_region.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o $(OBJ_DIR)/columns.o $(OBJ_DIR)/convergence.o
$(OBJ_DIR)/iteration_chemistry.o: iteration_chemistry.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o $(OBJ_DIR)/columns.o
$(OBJ_DIR)/level_population_solver.o: level_population_solver.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o $(OBJ_DIR)/convergence.o
$(OBJ_DIR)/evolution_setup.o: evolution_setup.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o $(OBJ_DIR)/chemistry_controls.o $(OBJ_DIR)/columns.o $(OBJ_DIR)/dark_region.o $(OBJ_DIR)/radiation.o $(OBJ_DIR)/iteration_chemistry.o $(OBJ_DIR)/level_population_solver.o
$(OBJ_DIR)/thermal_balance.o: thermal_balance.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o
$(OBJ_DIR)/iteration_convergence.o: iteration_convergence.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o $(OBJ_DIR)/convergence.o
$(OBJ_DIR)/output.o: output.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o
$(OBJ_DIR)/coolant_io.o: coolant_io.F90 $(OBJ_DIR)/modules.o
$(OBJ_DIR)/chemistry_io.o: chemistry_io.F90 $(OBJ_DIR)/modules.o
$(OBJ_DIR)/initial_conditions.o: initial_conditions.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/physical_parameters.o
$(OBJ_DIR)/geometry_setup.o: geometry_setup.F90 $(OBJ_DIR)/modules.o $(OBJ_DIR)/healpix_state.o
$(OBJ_DIR)/particle_storage.o: particle_storage.F90 $(OBJ_DIR)/modules.o
clean:
	rm -rf $(BUILD_DIR)
	rm -f 3DPDR 3DPDR.o test_columns test_convergence test_coolants test_excitation test_radiation
	rm -f *.mod *.o src/*.o src/*/*.o src/*/*/*.o cvode.log main.log fort.* HEALPix_vectors.dat run_3dpdr_check.log V1*.fin

compress:
	tar cvzf 3DPDR.tgz configs/*.params data/*.dat data/*.d makefile src/app/*.F90 src/io/*.F90 src/init/*.F90 src/evolution/*.F90 src/physics/*.F90 src/physics/chemistry/*.c
