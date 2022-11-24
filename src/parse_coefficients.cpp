#include "egravity.hpp"
#include "icgemio.hpp"
#include <cstdio>
#include <cassert>
#include <algorithm>

int dso::parse_gravity_model(const char *model_fn, int degree, int order,
                             dso::GravityField &grav,
                             bool denormalize) noexcept {
  dso::Icgem gfc(model_fn);

  // parse the header ...
  if (gfc.parse_header()) {
    fprintf(stderr,
            "[ERROR] Failed to parse icgem header for %s (traceback: %s)!\n",
            model_fn, __func__);
    return 1;
  }

  if (!(degree <= gfc.degree() && order <= degree)) {
    fprintf(stderr,
            "[ERROR]  invalid degree/order %d/%d for input gravity model %s "
            "(traceback: %s)\n",
            degree, order, model_fn, __func__);
    return 1;
  }

  // resize HarmonicCoeffs to fit input ...
  // If we have a TVG, count that first
  if (gfc.degree_tv_stop>0) {
    // we have a tvg
    assert(gfc.min_degree_tv == 0);
    if (degree<=gfc.degree_tv_stop) {
      // only TVG, no static part
      grav._tvg.resize(degree);
    } else {
      // TVG + static gravity
      grav._tvg.resize(gfc.degree_tv_stop);
      grav._static.resize(degree - gfc.degree_tv_stop + 1);
    }
  } else {
    // no TVG, only parse static field
    assert(gfc.min_degree_static == 0);
    grav._static.resize(degree);
    grav._tvg.resize(0);
  }
  
  // do we have harmonic coefficients ?
  grav._per.resize(
      gfc.harmonics.size(),
      (gfc.degree_tv_stop > 0) ? (std::min(gfc.degree_tv_stop, degree)) : (0));
  grav._per.copy_harmonics(gfc.harmonics);

  // parse data; store coefficients to harmonics
  if (gfc.parse_data(degree, order, &grav)) {
    fprintf(stderr,
            "[ERROR] Failed to parse harmonic coefficients from file %s "
            "(traceback: %s)\n",
            model_fn, __func__);
    return 1;
  }

  // if needed denormalize coefficients
  //if (denormalize)
  //  harmonics.denormalize();

  // all done
  return 0;
}
