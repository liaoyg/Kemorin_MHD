
/* take_normal_psf_c.h */

#ifndef TAKE_NORMAL_PSF_C_
#define TAKE_NORMAL_PSF_C_

#include <math.h>
#include <stdio.h>

#include "m_psf_data_4_viewer_c.h"
#include "cal_surface_center_normal_c.h"
#include "coordinate_converter_c.h"

/* prototype */
void take_normal_psf(struct psf_data *viz_s);
void take_minmax_psf(struct psf_data *viz_s);

void take_length_fline(struct psf_data *viz_s);
void take_minmax_fline(struct psf_data *viz_s);

#endif
