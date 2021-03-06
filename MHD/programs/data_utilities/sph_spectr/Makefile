#
#
#

SPH_UTILS_DIR = $$(DATA_UTILS_DIR)/sph_spectr

TARGET_SPH_ENE_SPEC =      ene_sph_spec
TARGET_SPH_ENE_LAYER =     ene_sph_layer
TARGET_PICK_SPH_VOL_ENE =  pick_sph_vol_ene
TARGET_PICK_SPH_ENE =      pick_sph_ene

TARGET_PICK_SPH_SPEC =     pick_sph_spec
TARGET_PICK_GAUSS_COEF =   pick_gauss_coefs

TARGET_RENAME_SPH_FLD =    rename_spectr_field
TARGET_DIFF_SPH_SPEC =     diff_spectr_data
TARGET_ZM_SPH_SPEC =       zonal_mean_sph
TARGET_T_AVE_SPEC =        t_ave_sph_spec
TARGET_SPH_NEUTRAL =       sph_neutral_points
TARGET_TAVE_SPH_ENE_SPEC = t_ave_sph_mean_square
TARGET_PSUM_SPH_ENE_SPEC = part_sum_sph_mean_square

TARGET_T_AVE_GAUSS_COEFS =   t_ave_gauss_coefs
TARGET_T_AVE_PICK_SPH_SPEC = t_ave_picked_sph_coefs
TARGET_T_AVE_PICK_SPH_ENE =  t_ave_picked_ene_sph
TARGET_T_AVE_NUSSELT =       t_ave_nusselt

TARGET_MAX_SPH_ENE =         sph_pwr_max_mode
TARGET_SPH_ULI_LSCALE =      sph_uli_lengscale

TARGET_SPH_ENE_CHECK =      sph_ene_check


LIBS_SPH_UTIL = \
'$$(LIB_VIZ)'    \
'$$(LIB_PARA)'   \
'$$(LIB_SOLVER)' \
'$$(LIB_MESH)'    \
'$$(LIB_FFTPACK)' \
'$$(LIB_ISPACK)'

LIB_FILES_SPH_UTIL = \
'$$(LIB_VIZ_FILE)'    \
'$$(LIB_PARA_FILE)'   \
'$$(LIB_SOLVER_FILE)' \
'$$(LIB_MESH_FILE)'   \
'$$(LIB_FFTPACK_FILE)' \
'$$(LIB_ISPACK_FILE)'

ifdef FORTRAN_C_CONNECT
  LIBS_SPH_UTIL+=        '$$(LIB_KEMO_C)'
  LIB_FILES_SPH_UTIL+=   '$$(LIB_KEMO_C_FILE)'
endif

SOURCES = $(shell ls *.f90)

MOD_TAVE_SPH_SPEC = \
m_spheric_data_sph_spetr.o \
m_ctl_data_4_sph_utils.o \
m_ctl_params_sph_utils.o \
\
cal_t_ave_sph_spectr_data.o \
analyzer_time_ave_sph.o \
main_time_ave_sph.o

MOD_SPH_ENE_SPEC = \
m_spheric_data_sph_spetr.o \
m_ctl_data_4_sph_utils.o \
m_ctl_params_sph_utils.o \
\
analyzer_ene_sph_shell.o \
main_ene_sph_shell.o

MOD_SPH_ENE_LAYER = \
m_spheric_data_sph_spetr.o \
m_ctl_data_4_sph_utils.o \
m_ctl_params_sph_utils.o \
\
analyzer_ene_sph_layer.o \
main_ene_sph_layer.o

MOD_PICK_SPH_VOL_ENE = \
m_spheric_data_sph_spetr.o \
m_ctl_data_4_sph_utils.o \
m_ctl_params_sph_utils.o \
\
pickup_sph_rms_spectr.o \
analyzer_pick_rms_vol.o \
main_pick_rms_vol.o

MOD_PICK_SPH_ENE = \
m_spheric_data_sph_spetr.o \
m_ctl_data_4_sph_utils.o \
m_ctl_params_sph_utils.o \
\
pickup_sph_rms_spectr.o \
analyzer_pick_rms_sph.o \
main_pick_rms_sph.o

MOD_PICK_SPH_SPEC = \
m_spheric_data_sph_spetr.o \
m_ctl_data_4_sph_utils.o \
m_ctl_params_sph_utils.o \
\
analyzer_pickup_mode_sph.o \
main_pickup_mode_sph.o

MOD_PICK_GAUSS_COEF = \
m_spheric_data_sph_spetr.o \
m_ctl_data_4_sph_utils.o \
m_ctl_params_sph_utils.o \
\
analyzer_pickup_gauss_coefs.o \
main_pickup_gauss_coefs.o

MOD_RENAME_SPH_FLD =  \
t_ctl_data_sph_data_utils.o \
t_ctl_param_sph_data_utils.o \
\
analyzer_rename_sph_field.o \
main_rename_sph_field.o

MOD_DIFF_SPH_SPEC =  \
t_ctl_data_sph_data_utils.o \
t_ctl_param_sph_data_utils.o \
\
analyzer_diff_sph_spectr.o \
main_diff_sph_spectr.o

MOD_ZM_SPH_SPEC = \
m_spheric_data_sph_spetr.o \
m_ctl_data_4_sph_utils.o \
m_ctl_params_sph_utils.o \
\
analyzer_zonal_mean_sph.o \
main_zonal_mean_sph.o

MOD_SPH_NEUTRAL = \
neutral_point_by_pick_sph.o \
m_spheric_data_sph_spetr.o \
m_neutral_pt_by_pick_sph.o


MOD_TAVE_SPH_ENE_SPEC = \
t_average_sph_ene_spec.o \
t_read_sph_spectra.o \
m_tave_sph_ene_spectr.o \
sph_mean_square_IO_select.o \
simple_sph_spectr_head_IO.o \
simple_sph_spectr_data_IO.o \
cal_tave_sph_ene_spectr.o

MOD_MAX_SPH_ENE = \
maxmode_sph_ene_spec.o \
t_read_sph_spectra.o \
m_maxmode_sph_ene_spectr.o \
sph_mean_square_IO_select.o \
simple_sph_spectr_head_IO.o \
simple_sph_spectr_data_IO.o

MOD_SPH_ULI_LSCALE = \
sph_uli_lengh_scale.o \
t_read_sph_spectra.o \
m_sph_uli_lengh_scale.o \
sph_mean_square_IO_select.o \
simple_sph_spectr_head_IO.o \
simple_sph_spectr_data_IO.o

MOD_PSUM_SPH_ENE_SPEC = \
t_read_sph_spectra.o \
m_part_sum_sph_ene_spectr.o \
sph_mean_square_IO_select.o \
simple_sph_spectr_head_IO.o \
simple_sph_spectr_data_IO.o \
part_sum_sph_ene_spec.o

MOD_T_AVE_GAUSS_COEFS = \
tave_picked_gauss_coefs.o

MOD_T_AVE_PICK_SPH = \
tave_picked_sph_spec_data.o

MOD_T_AVE_PICK_SPH_ENE = \
tave_picked_sph_rms_data.o

MOD_T_AVE_NUSSELT = \
t_average_nusselt.o

MOD_SPH_ENE_CHECK = \
compare_sph_mean_square.o

#
#  ------------------------------------------------------------------
#

dir_list:
	@echo 'SPH_UTILS_DIR = $(SPH_UTILS_DIR)'     >> $(MAKENAME)

target_list:
	@echo 'TARGET_SPH_ENE_SPEC =     $$(BUILDDIR)/$(TARGET_SPH_ENE_SPEC)' \
	>> $(MAKENAME)
	@echo 'TARGET_SPH_ENE_LAYER =    $$(BUILDDIR)/$(TARGET_SPH_ENE_LAYER)' \
	>> $(MAKENAME)
	@echo 'TARGET_PICK_SPH_VOL_ENE = $$(BUILDDIR)/$(TARGET_PICK_SPH_VOL_ENE)' \
	>> $(MAKENAME)
	@echo 'TARGET_PICK_SPH_ENE =     $$(BUILDDIR)/$(TARGET_PICK_SPH_ENE)' \
	>> $(MAKENAME)
	@echo 'TARGET_PICK_SPH_SPEC =    $$(BUILDDIR)/$(TARGET_PICK_SPH_SPEC)' \
	>> $(MAKENAME)
	@echo 'TARGET_PICK_GAUSS_COEF =  $$(BUILDDIR)/$(TARGET_PICK_GAUSS_COEF)' \
	>> $(MAKENAME)
	@echo 'TARGET_RENAME_SPH_FLD =   $$(BUILDDIR)/$(TARGET_RENAME_SPH_FLD)' \
	>> $(MAKENAME)
	@echo 'TARGET_DIFF_SPH_SPEC =    $$(BUILDDIR)/$(TARGET_DIFF_SPH_SPEC)' \
	>> $(MAKENAME)
	@echo 'TARGET_ZM_SPH_SPEC =      $$(BUILDDIR)/$(TARGET_ZM_SPH_SPEC)' \
	>> $(MAKENAME)
	@echo 'TARGET_T_AVE_SPEC =       $$(BUILDDIR)/$(TARGET_T_AVE_SPEC)' \
	>> $(MAKENAME)
	@echo 'TARGET_SPH_NEUTRAL =      $$(BUILDDIR)/$(TARGET_SPH_NEUTRAL)' \
	>> $(MAKENAME)
	@echo 'TARGET_TAVE_SPH_ENE_SPEC = $$(BUILDDIR)/$(TARGET_TAVE_SPH_ENE_SPEC)' \
	>> $(MAKENAME)
	@echo 'TARGET_T_AVE_GAUSS_COEFS = $$(BUILDDIR)/$(TARGET_T_AVE_GAUSS_COEFS)' \
	>> $(MAKENAME)
	@echo 'TARGET_T_AVE_PICK_SPH_SPEC = $$(BUILDDIR)/$(TARGET_T_AVE_PICK_SPH_SPEC)' \
	>> $(MAKENAME)
	@echo 'TARGET_T_AVE_PICK_SPH_ENE = $$(BUILDDIR)/$(TARGET_T_AVE_PICK_SPH_ENE)' \
	>> $(MAKENAME)
	@echo 'TARGET_T_AVE_NUSSELT = $$(BUILDDIR)/$(TARGET_T_AVE_NUSSELT)' \
	>> $(MAKENAME)
	@echo 'TARGET_MAX_SPH_ENE = $$(BUILDDIR)/$(TARGET_MAX_SPH_ENE)' \
	>> $(MAKENAME)
	@echo 'TARGET_SPH_ULI_LSCALE = $$(BUILDDIR)/$(TARGET_SPH_ULI_LSCALE)' \
	>> $(MAKENAME)
	@echo 'TARGET_SPH_ENE_CHECK = $$(BUILDDIR)/$(TARGET_SPH_ENE_CHECK)' \
	>> $(MAKENAME)
	@echo 'TARGET_PSUM_SPH_ENE_SPEC = $$(BUILDDIR)/$(TARGET_PSUM_SPH_ENE_SPEC)' \
	>> $(MAKENAME)
	@echo >> $(MAKENAME)

target_task:
	@echo sph_data_util: \
	'$$(TARGET_SPH_ENE_SPEC)       $$(TARGET_SPH_ENE_LAYER)' \
	'$$(TARGET_PICK_SPH_VOL_ENE)   $$(TARGET_PICK_SPH_ENE)'  \
	'$$(TARGET_PICK_SPH_SPEC)'                               \
	'$$(TARGET_RENAME_SPH_FLD)     $$(TARGET_DIFF_SPH_SPEC)'  \
	'$$(TARGET_ZM_SPH_SPEC)        $$(TARGET_T_AVE_SPEC)'      \
	'$$(TARGET_SPH_NEUTRAL)        $$(TARGET_TAVE_SPH_ENE_SPEC)' \
	'$$(TARGET_PICK_GAUSS_COEF)    $$(TARGET_T_AVE_PICK_SPH_SPEC)' \
	'$$(TARGET_T_AVE_GAUSS_COEFS)  $$(TARGET_T_AVE_PICK_SPH_ENE)' \
	'$$(TARGET_MAX_SPH_ENE)        $$(TARGET_SPH_ULI_LSCALE)'     \
	'$$(TARGET_T_AVE_NUSSELT)      $$(TARGET_SPH_ENE_CHECK)'      \
	'$$(TARGET_PSUM_SPH_ENE_SPEC) ' >> $(MAKENAME)
	@echo >> $(MAKENAME)
	@echo '$$(TARGET_SPH_ENE_SPEC)': '$$(MOD_SPH_ENE_SPEC)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_SPH_ENE_SPEC)' \
	'$$(MOD_SPH_ENE_SPEC)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_SPH_ENE_LAYER)': '$$(MOD_SPH_ENE_LAYER)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_SPH_ENE_LAYER)' \
	'$$(MOD_SPH_ENE_LAYER)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_PICK_SPH_VOL_ENE)': '$$(MOD_PICK_SPH_VOL_ENE)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_PICK_SPH_VOL_ENE)' \
	'$$(MOD_PICK_SPH_VOL_ENE)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_PICK_SPH_ENE)': '$$(MOD_PICK_SPH_ENE)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_PICK_SPH_ENE)' \
	'$$(MOD_PICK_SPH_ENE)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_PICK_SPH_SPEC)': '$$(MOD_PICK_SPH_SPEC)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_PICK_SPH_SPEC)' \
	'$$(MOD_PICK_SPH_SPEC)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_PICK_GAUSS_COEF)': '$$(MOD_PICK_GAUSS_COEF)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_PICK_GAUSS_COEF)' \
	'$$(MOD_PICK_GAUSS_COEF)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_RENAME_SPH_FLD)': '$$(MOD_RENAME_SPH_FLD)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_RENAME_SPH_FLD)' \
	'$$(MOD_RENAME_SPH_FLD)' \
	'-L. $$(LIBS_SPH_UTIL)  $$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_DIFF_SPH_SPEC)': '$$(MOD_DIFF_SPH_SPEC)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_DIFF_SPH_SPEC)' \
	'$$(MOD_DIFF_SPH_SPEC)' \
	'-L. $$(LIBS_SPH_UTIL)  $$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_ZM_SPH_SPEC)': '$$(MOD_ZM_SPH_SPEC)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_ZM_SPH_SPEC)' \
	'$$(MOD_ZM_SPH_SPEC)' \
	'-L. $$(LIBS_SPH_UTIL)  $$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_T_AVE_SPEC)': '$$(MOD_TAVE_SPH_SPEC)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_T_AVE_SPEC)' \
	'$$(MOD_TAVE_SPH_SPEC)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_SPH_NEUTRAL)': '$$(MOD_SPH_NEUTRAL)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_SPH_NEUTRAL)' \
	'$$(MOD_SPH_NEUTRAL)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_TAVE_SPH_ENE_SPEC)': '$$(MOD_TAVE_SPH_ENE_SPEC)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_TAVE_SPH_ENE_SPEC)' \
	'$$(MOD_TAVE_SPH_ENE_SPEC)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_T_AVE_GAUSS_COEFS)': '$$(MOD_T_AVE_GAUSS_COEFS)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_T_AVE_GAUSS_COEFS)' \
	'$$(MOD_T_AVE_GAUSS_COEFS)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_T_AVE_PICK_SPH_SPEC)': '$$(MOD_T_AVE_PICK_SPH)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_T_AVE_PICK_SPH_SPEC)' \
	'$$(MOD_T_AVE_PICK_SPH)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_PSUM_SPH_ENE_SPEC)': '$$(MOD_PSUM_SPH_ENE_SPEC)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_PSUM_SPH_ENE_SPEC)' \
	'$$(MOD_PSUM_SPH_ENE_SPEC)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_T_AVE_PICK_SPH_ENE)': '$$(MOD_T_AVE_PICK_SPH_ENE)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_T_AVE_PICK_SPH_ENE)' \
	'$$(MOD_T_AVE_PICK_SPH_ENE)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_T_AVE_NUSSELT)': '$$(MOD_T_AVE_NUSSELT)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_T_AVE_NUSSELT)' \
	'$$(MOD_T_AVE_NUSSELT)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_MAX_SPH_ENE)': '$$(MOD_MAX_SPH_ENE)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_MAX_SPH_ENE)' \
	'$$(MOD_MAX_SPH_ENE)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_SPH_ULI_LSCALE)': '$$(MOD_SPH_ULI_LSCALE)' \
	'$$(LIB_FILES_SPH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_SPH_ULI_LSCALE)' \
	'$$(MOD_SPH_ULI_LSCALE)' \
	'-L. $$(LIBS_SPH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_SPH_ENE_CHECK)': '$$(MOD_SPH_ENE_CHECK)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_SPH_ENE_CHECK)' \
	'$$(MOD_SPH_ENE_CHECK)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)


lib_name:
	@echo  LIBS_SPH_UTIL=  \\          >> $(MAKENAME)
	@echo  $(LIBS_SPH_UTIL)            >> $(MAKENAME)
	@echo  LIB_FILES_SPH_UTIL=  \\     >> $(MAKENAME)
	@echo  $(LIB_FILES_SPH_UTIL)       >> $(MAKENAME)
	@echo  '#' >> $(MAKENAME)

mod_list:
	@echo  MOD_SPH_ENE_SPEC=  \\     >> $(MAKENAME)
	@echo  $(MOD_SPH_ENE_SPEC)       >> $(MAKENAME)
	@echo  MOD_SPH_ENE_LAYER=  \\    >> $(MAKENAME)
	@echo  $(MOD_SPH_ENE_LAYER)      >> $(MAKENAME)
	@echo  MOD_PICK_SPH_VOL_ENE=  \\ >> $(MAKENAME)
	@echo  $(MOD_PICK_SPH_VOL_ENE)   >> $(MAKENAME)
	@echo  MOD_PICK_SPH_ENE=  \\     >> $(MAKENAME)
	@echo  $(MOD_PICK_SPH_ENE)       >> $(MAKENAME)
	@echo  MOD_PICK_SPH_SPEC=  \\    >> $(MAKENAME)
	@echo  $(MOD_PICK_SPH_SPEC)      >> $(MAKENAME)
	@echo  MOD_PICK_GAUSS_COEF=  \\  >> $(MAKENAME)
	@echo  $(MOD_PICK_GAUSS_COEF)    >> $(MAKENAME)
	@echo  MOD_RENAME_SPH_FLD=  \\   >> $(MAKENAME)
	@echo  $(MOD_RENAME_SPH_FLD)     >> $(MAKENAME)
	@echo  MOD_DIFF_SPH_SPEC=  \\    >> $(MAKENAME)
	@echo  $(MOD_DIFF_SPH_SPEC)      >> $(MAKENAME)
	@echo  MOD_ZM_SPH_SPEC=  \\      >> $(MAKENAME)
	@echo  $(MOD_ZM_SPH_SPEC)        >> $(MAKENAME)
	@echo  MOD_TAVE_SPH_SPEC=  \\    >> $(MAKENAME)
	@echo  $(MOD_TAVE_SPH_SPEC)      >> $(MAKENAME)
	@echo  MOD_SPH_NEUTRAL=  \\      >> $(MAKENAME)
	@echo  $(MOD_SPH_NEUTRAL)        >> $(MAKENAME)
	@echo  MOD_TAVE_SPH_ENE_SPEC=  \\  >> $(MAKENAME)
	@echo  $(MOD_TAVE_SPH_ENE_SPEC)    >> $(MAKENAME)
	@echo  MOD_T_AVE_GAUSS_COEFS=  \\  >> $(MAKENAME)
	@echo  $(MOD_T_AVE_GAUSS_COEFS)    >> $(MAKENAME)
	@echo  MOD_T_AVE_PICK_SPH=  \\     >> $(MAKENAME)
	@echo  $(MOD_T_AVE_PICK_SPH)       >> $(MAKENAME)
	@echo  MOD_T_AVE_PICK_SPH_ENE=  \\ >> $(MAKENAME)
	@echo  $(MOD_T_AVE_PICK_SPH_ENE)   >> $(MAKENAME)
	@echo  MOD_MAX_SPH_ENE=  \\        >> $(MAKENAME)
	@echo  $(MOD_MAX_SPH_ENE)          >> $(MAKENAME)
	@echo  MOD_SPH_ULI_LSCALE=  \\     >> $(MAKENAME)
	@echo  $(MOD_SPH_ULI_LSCALE)       >> $(MAKENAME)
	@echo  MOD_PSUM_SPH_ENE_SPEC=  \\  >> $(MAKENAME)
	@echo  $(MOD_PSUM_SPH_ENE_SPEC)    >> $(MAKENAME)
	@echo  MOD_T_AVE_NUSSELT=  \\      >> $(MAKENAME)
	@echo  $(MOD_T_AVE_NUSSELT)        >> $(MAKENAME)
	@echo  MOD_SPH_ENE_CHECK=  \\      >> $(MAKENAME)
	@echo  $(MOD_SPH_ENE_CHECK)        >> $(MAKENAME)
	@echo   >> $(MAKENAME)

module:
	@$(MAKE_MOD_DEP) '$(MAKENAME)' '$$(SPH_UTILS_DIR)' $(SOURCES)
