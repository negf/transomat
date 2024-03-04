# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include Makefile.arch


.DEFAULT:
	-echo $@ does not exist.
all: $(BIN_DIR)/transomat
$(OBJ_DIR)/analysis_mod.o: src/analysis_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/k_on_demand_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/analysis_mod.f90 -o $(OBJ_DIR)/analysis_mod.o
$(OBJ_DIR)/check_sys.o: src/check_sys.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/check_sys.f90 -o $(OBJ_DIR)/check_sys.o
$(OBJ_DIR)/conquest_mod.o: src/conquest_mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/conquest_mod.f90 -o $(OBJ_DIR)/conquest_mod.o
$(OBJ_DIR)/control.o: src/control.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/control.f90 -o $(OBJ_DIR)/control.o
$(OBJ_DIR)/dft_sigma_mod.o: src/dft_sigma_mod.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/dft_sigma_mod.f90 -o $(OBJ_DIR)/dft_sigma_mod.o
$(OBJ_DIR)/dftu_mod.o: src/dftu_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/k_on_demand_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/petsc_wrapper_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/dftu_mod.f90 -o $(OBJ_DIR)/dftu_mod.o
$(OBJ_DIR)/eigenchannel_mod.o: src/eigenchannel_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/slepc_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/eigenchannel_mod.f90 -o $(OBJ_DIR)/eigenchannel_mod.o
$(OBJ_DIR)/error_handler_mod.o: src/error_handler_mod.f90 $(OBJ_DIR)/kinds.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/error_handler_mod.f90 -o $(OBJ_DIR)/error_handler_mod.o
$(OBJ_DIR)/fix_mat_mod.o: src/fix_mat_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/fix_mat_mod.f90 -o $(OBJ_DIR)/fix_mat_mod.o
$(OBJ_DIR)/ft_mod.o: src/ft_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/ft_mod.f90 -o $(OBJ_DIR)/ft_mod.o
$(OBJ_DIR)/globals.o: src/globals.f90 $(OBJ_DIR)/kinds.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/globals.f90 -o $(OBJ_DIR)/globals.o
$(OBJ_DIR)/init_control.o: src/init_control.f90 $(OBJ_DIR)/globals.o $(OBJ_DIR)/control.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/init_control.f90 -o $(OBJ_DIR)/init_control.o
$(OBJ_DIR)/init_ep_coupling_mod.o: src/init_ep_coupling_mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/init_ep_coupling_mod.f90 -o $(OBJ_DIR)/init_ep_coupling_mod.o
$(OBJ_DIR)/init_petsc.o: src/init_petsc.f90 $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/init_petsc.f90 -o $(OBJ_DIR)/init_petsc.o
$(OBJ_DIR)/init_sys.o: src/init_sys.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/read_sysinfo.o $(OBJ_DIR)/lattice_mod.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/conquest_mod.o $(OBJ_DIR)/kmat_mod.o $(OBJ_DIR)/slepc_mod.o $(OBJ_DIR)/init_ep_coupling_mod.o $(OBJ_DIR)/phonon_mod.o $(OBJ_DIR)/loe_mod.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/k_on_demand_mod.o $(OBJ_DIR)/dftu_mod.o $(OBJ_DIR)/fix_mat_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/init_sys.f90 -o $(OBJ_DIR)/init_sys.o
$(OBJ_DIR)/integration_weights.o: src/integration_weights.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/integration_weights.f90 -o $(OBJ_DIR)/integration_weights.o
$(OBJ_DIR)/integrator.o: src/integrator.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/math_stuff.o $(OBJ_DIR)/surf_gf.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/integration_weights.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/pexsi_wrapper.o $(OBJ_DIR)/ft_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/integrator.f90 -o $(OBJ_DIR)/integrator.o
$(OBJ_DIR)/integrator_scalar.o: src/integrator_scalar.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/math_stuff.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/integrator.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/integrator_scalar.f90 -o $(OBJ_DIR)/integrator_scalar.o
$(OBJ_DIR)/j_dump_mod.o: src/j_dump_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/write_conquest_dump.mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/j_dump_mod.f90 -o $(OBJ_DIR)/j_dump_mod.o
$(OBJ_DIR)/k_on_demand_mod.o: src/k_on_demand_mod.f90 $(OBJ_DIR)/globals.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/petsc_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/k_on_demand_mod.f90 -o $(OBJ_DIR)/k_on_demand_mod.o
$(OBJ_DIR)/kinds.o: src/kinds.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/kinds.f90 -o $(OBJ_DIR)/kinds.o
$(OBJ_DIR)/kmat_from_diag.o: src/kmat_from_diag.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/k_on_demand_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/slepc_mod.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/dftu_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/kmat_from_diag.f90 -o $(OBJ_DIR)/kmat_from_diag.o
$(OBJ_DIR)/kmat_mod.o: src/kmat_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/write_conquest_dump.mod.o $(OBJ_DIR)/misc.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/kmat_mod.f90 -o $(OBJ_DIR)/kmat_mod.o
$(OBJ_DIR)/lattice_mod.o: src/lattice_mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/math_stuff.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/lattice_mod.f90 -o $(OBJ_DIR)/lattice_mod.o
$(OBJ_DIR)/loe_mod.o: src/loe_mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/phonon_mod.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/integrator_scalar.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/loe_mod.f90 -o $(OBJ_DIR)/loe_mod.o
$(OBJ_DIR)/math_stuff.o: src/math_stuff.f90 $(OBJ_DIR)/kinds.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/math_stuff.f90 -o $(OBJ_DIR)/math_stuff.o
$(OBJ_DIR)/misc.o: src/misc.f90 $(OBJ_DIR)/kinds.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/misc.f90 -o $(OBJ_DIR)/misc.o
$(OBJ_DIR)/petsc_mod.o: src/petsc_mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/math_stuff.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/petsc_mod.f90 -o $(OBJ_DIR)/petsc_mod.o
$(OBJ_DIR)/petsc_test.o: src/petsc_test.f90 $(OBJ_DIR)/slepc_mod.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/petsc_test.f90 -o $(OBJ_DIR)/petsc_test.o
$(OBJ_DIR)/petsc_wrapper_mod.o: src/petsc_wrapper_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/timing_mod.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/slepc_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/petsc_wrapper_mod.f90 -o $(OBJ_DIR)/petsc_wrapper_mod.o
$(OBJ_DIR)/pexsi_wrapper.o: src/pexsi_wrapper.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/globals.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/pexsi_wrapper.f90 -o $(OBJ_DIR)/pexsi_wrapper.o
$(OBJ_DIR)/phonon_mod.o: src/phonon_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/slepc_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/phonon_mod.f90 -o $(OBJ_DIR)/phonon_mod.o
$(OBJ_DIR)/read_sysinfo.o: src/read_sysinfo.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/lattice_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/read_sysinfo.f90 -o $(OBJ_DIR)/read_sysinfo.o
$(OBJ_DIR)/scf_mod.o: src/scf_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/kmat_mod.o $(OBJ_DIR)/timing_mod.o $(OBJ_DIR)/k_on_demand_mod.o $(OBJ_DIR)/j_dump_mod.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/scf_mod.f90 -o $(OBJ_DIR)/scf_mod.o
$(OBJ_DIR)/slepc_mod.o: src/slepc_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/slepc_mod.f90 -o $(OBJ_DIR)/slepc_mod.o
$(OBJ_DIR)/surf_gf.o: src/surf_gf.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/math_stuff.o $(OBJ_DIR)/timing_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/surf_gf.f90 -o $(OBJ_DIR)/surf_gf.o
$(OBJ_DIR)/timing_mod.o: src/timing_mod.f90 $(OBJ_DIR)/kinds.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/timing_mod.f90 -o $(OBJ_DIR)/timing_mod.o
$(OBJ_DIR)/transmission_mod.o: src/transmission_mod.f90 $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/slepc_mod.o $(OBJ_DIR)/eigenchannel_mod.o $(OBJ_DIR)/k_on_demand_mod.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/transmission_mod.f90 -o $(OBJ_DIR)/transmission_mod.o
$(OBJ_DIR)/transomat.o: src/transomat.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/init_sys.o $(OBJ_DIR)/init_petsc.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/scf_mod.o $(OBJ_DIR)/loe_mod.o $(OBJ_DIR)/phonon_mod.o $(OBJ_DIR)/transmission_mod.o $(OBJ_DIR)/kmat_from_diag.o $(OBJ_DIR)/kmat_mod.o $(OBJ_DIR)/analysis_mod.o $(OBJ_DIR)/dftu_mod.o $(OBJ_DIR)/timing_mod.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/transomat.f90 -o $(OBJ_DIR)/transomat.o
$(OBJ_DIR)/write_conquest_dump.mod.o: src/write_conquest_dump.mod.f90 $(OBJ_DIR)/kinds.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/error_handler_mod.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	src/write_conquest_dump.mod.f90 -o $(OBJ_DIR)/write_conquest_dump.mod.o
SRC = src/petsc_wrapper_mod.f90 src/j_dump_mod.f90 src/globals.f90 src/loe_mod.f90 src/init_petsc.f90 src/ft_mod.f90 src/read_sysinfo.f90 src/fix_mat_mod.f90 src/analysis_mod.f90 src/integrator_scalar.f90 src/math_stuff.f90 src/eigenchannel_mod.f90 src/petsc_test.f90 src/scf_mod.f90 src/timing_mod.f90 src/conquest_mod.f90 src/kinds.f90 src/transomat.f90 src/lattice_mod.f90 src/dftu_mod.f90 src/phonon_mod.f90 src/check_sys.f90 src/write_conquest_dump.mod.f90 src/init_sys.f90 src/surf_gf.f90 src/dft_sigma_mod.f90 src/control.f90 src/petsc_mod.f90 src/error_handler_mod.f90 src/pexsi_wrapper.f90 src/slepc_mod.f90 src/k_on_demand_mod.f90 src/transmission_mod.f90 src/kmat_from_diag.f90 src/init_ep_coupling_mod.f90 src/misc.f90 src/integration_weights.f90 src/init_control.f90 src/integrator.f90 src/kmat_mod.f90
OBJ = $(OBJ_DIR)/petsc_wrapper_mod.o $(OBJ_DIR)/j_dump_mod.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/loe_mod.o $(OBJ_DIR)/init_petsc.o $(OBJ_DIR)/ft_mod.o $(OBJ_DIR)/read_sysinfo.o $(OBJ_DIR)/fix_mat_mod.o $(OBJ_DIR)/analysis_mod.o $(OBJ_DIR)/integrator_scalar.o $(OBJ_DIR)/math_stuff.o $(OBJ_DIR)/eigenchannel_mod.o $(OBJ_DIR)/petsc_test.o $(OBJ_DIR)/scf_mod.o $(OBJ_DIR)/timing_mod.o $(OBJ_DIR)/conquest_mod.o $(OBJ_DIR)/kinds.o $(OBJ_DIR)/transomat.o $(OBJ_DIR)/lattice_mod.o $(OBJ_DIR)/dftu_mod.o $(OBJ_DIR)/phonon_mod.o $(OBJ_DIR)/check_sys.o $(OBJ_DIR)/write_conquest_dump.mod.o $(OBJ_DIR)/init_sys.o $(OBJ_DIR)/surf_gf.o $(OBJ_DIR)/dft_sigma_mod.o $(OBJ_DIR)/control.o $(OBJ_DIR)/petsc_mod.o $(OBJ_DIR)/error_handler_mod.o $(OBJ_DIR)/pexsi_wrapper.o $(OBJ_DIR)/slepc_mod.o $(OBJ_DIR)/k_on_demand_mod.o $(OBJ_DIR)/transmission_mod.o $(OBJ_DIR)/kmat_from_diag.o $(OBJ_DIR)/init_ep_coupling_mod.o $(OBJ_DIR)/misc.o $(OBJ_DIR)/integration_weights.o $(OBJ_DIR)/init_control.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/kmat_mod.o
OFF = src/k_on_demand_mod.f90 src/misc.f90 src/init_sys.f90 src/loe_mod.f90 src/petsc_test.f90 src/petsc_mod.f90 src/init_petsc.f90 src/lattice_mod.f90 src/init_control.f90 src/ft_mod.f90 src/check_sys.f90 src/fix_mat_mod.f90 src/integrator.f90 src/dftu_mod.f90 src/control.f90 src/eigenchannel_mod.f90 src/pexsi_wrapper.f90 src/conquest_mod.f90 src/read_sysinfo.f90 src/kinds.f90 src/scf_mod.f90 src/analysis_mod.f90 src/phonon_mod.f90 src/timing_mod.f90 src/error_handler_mod.f90 src/globals.f90 src/dft_sigma_mod.f90 src/slepc_mod.f90 src/j_dump_mod.f90 src/write_conquest_dump.mod.f90 src/integration_weights.f90 src/transomat.f90 src/init_ep_coupling_mod.f90 src/math_stuff.f90 src/integrator_scalar.f90 src/kmat_from_diag.f90 src/transmission_mod.f90 src/kmat_mod.f90 src/surf_gf.f90 src/petsc_wrapper_mod.f90
clean: neat
	-rm -f .transomat.cppdefs $(OBJ) $(MOD_DIR)/*.mod $(BIN_DIR)/transomat
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
$(OBJ_DIR):
	mkdir $(OBJ_DIR)
$(MOD_DIR):
	mkdir $(MOD_DIR)
$(BIN_DIR):
	mkdir $(BIN_DIR)
$(BIN_DIR)/transomat: $(OBJ) | $(BIN_DIR) $(OBJ_DIR) 
	$(LD) $(OBJ) -o $(BIN_DIR)/transomat  $(LDFLAGS)
