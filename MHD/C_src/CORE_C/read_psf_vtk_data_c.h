//
//  read_psf_vtk_data_c.h
//  Kemoview_Cocoa
//
//  Created by Hiroaki Matsui on 13/01/12.
//
//

#ifndef READ_PSF_VTK_DATA_C__
#define READ_PSF_VTK_DATA_C__

#include "kemosrc_param_c.h"
#include "m_psf_data_4_viewer_c.h"
#include "skip_comment_c.h"

/* prototypes */

int read_psf_vtg(const char *file_head, struct psf_data *viz_s);
int read_psf_vtd(const char *file_head, int istep, struct psf_data *viz_s);
int read_kemoview_vtk(const char *file_head, struct psf_data *viz_s);

#endif
