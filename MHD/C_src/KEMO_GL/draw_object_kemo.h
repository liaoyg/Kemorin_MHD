
/* draw_object_kemo.h */

#ifndef DEAW_OBJECT_KEMO_
#define DEAW_OBJECT_KEMO_

#include "m_kemoviewer_menu.h"
#include "m_gl_transfer_matrix.h"
#include "draw_patches_4_PSF.h"
#include "draw_isolines_4_PSF.h"
#include "draw_isolines_4_map.h"
#include "draw_node_by_ico_c.h"
#include "draw_grids_4_mesh.h"
#include "draw_patch_4_mesh_c.h"
#include "draw_colorbar_gl.h"
#include "sort_by_patch_distance.h"
#include "init_gl_lighting_c.h"
#include "draw_mapgrid.h"
#include "draw_coastline.h"


/* prototypes */

int draw_objects_4_map(struct psf_data **psf_s, struct mesh_menu_val *mesh_m,
                        struct psf_menu_val **psf_m, struct kemo_array_control *psf_a,
                        struct view_element *view_s, struct buffer_for_gl *gl_buf);
	
void draw_nodes_4_domain(struct viewer_mesh *mesh_s, struct mesh_menu_val *mesh_m,
						 struct buffer_for_gl *gl_buf);
void draw_grids_4_domain(struct viewer_mesh *mesh_s, struct mesh_menu_val *mesh_m,
						 struct buffer_for_gl *gl_buf);
void draw_patches_4_domain(struct viewer_mesh *mesh_s, struct mesh_menu_val *mesh_m, 
						   struct buffer_for_gl *gl_buf);
void draw_transparent_4_domain(struct viewer_mesh *mesh_s, struct mesh_menu_val *mesh_m,
			struct view_element *view_s, struct buffer_for_gl *gl_buf);

#endif
