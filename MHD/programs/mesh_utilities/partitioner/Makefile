#
#
#

PART_MAINDIR =     $$(MESHER_DIR)/partitioner

TARGET_PART_KEMO =  part_kemo
TARGET_ITP_PART =   part_interpolate_table

SOURCES = $(shell ls *.f90)

MOD_PARTITIONER = \
partitioner.o \
set_partition_by_fine_mesh.o \
grouping_for_partition.o

MOD_ITP_PART = \
distribute_para_itp_table.o


#
#  -------------------------------------------------------------------------
#

dir_list:
	@echo 'PART_MAINDIR =     $$(MESHER_DIR)/partitioner'     >> $(MAKENAME)

target_list:
	@echo 'TARGET_PART_KEMO = $$(BUILDDIR)/$(TARGET_PART_KEMO)' >> $(MAKENAME)
	@echo 'TARGET_ITP_PART =  $$(BUILDDIR)/$(TARGET_ITP_PART)' >> $(MAKENAME)
	@echo >> $(MAKENAME)

target_task:
	@echo partitioners: \
	'$$(TARGET_PART_KEMO) $$(TARGET_ITP_PART) ' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_PART_KEMO)': '$$(MOD_PARTITIONER)' \
	'$$(LIB_FILES_MESH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_PART_KEMO)' \
	'$$(MOD_PARTITIONER)' '-L. $$(LIBS_MESH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)
	@echo '$$(TARGET_ITP_PART)': '$$(MOD_ITP_PART)' \
	'$$(LIB_FILES_MESH_UTIL)' \
	>> $(MAKENAME)
	@echo '	''$$(F90)' '$$(F90FLAGS)' -o '$$(TARGET_ITP_PART)' \
	'$$(MOD_ITP_PART)' '-L. $$(LIBS_MESH_UTIL)' \
	'$$(F90LIBS)' >> $(MAKENAME)
	@echo '' >> $(MAKENAME)


lib_name:

mod_list:
	@echo  MOD_PARTITIONER =  \\ >> $(MAKENAME)
	@echo  $(MOD_PARTITIONER)    >> $(MAKENAME)
	@echo  MOD_ITP_PART =  \\ >> $(MAKENAME)
	@echo  $(MOD_ITP_PART)    >> $(MAKENAME)
	@echo   >> $(MAKENAME)


module:
	@$(MAKE_MOD_DEP) '$(MAKENAME)' '$$(PART_MAINDIR)' $(SOURCES)

clean:
	rm -f *.o *~ *.mod *.diag *.L
