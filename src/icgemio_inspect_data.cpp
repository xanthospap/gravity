#include "icgemio.hpp"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

/// @brief Check if a given line (string) starts with a given pattern
/// @param[in] pattern A null-terminated string, the pattern we search for
/// @param[in] line    A null-terminated string, the string we search (against
///                    the pattern)
/// @return Returns true if the first N characters of line are exactly equal 
/// to the pattern string; N is the number of characters in pattern, excluding 
/// the NULL character.
namespace {
bool starts_with(const char *pattern, const char *line) noexcept {
  if (line)
    return !std::strncmp(pattern, line, std::strlen(pattern));
  return false;
}
}

int dso::Icgem::inspect_data() noexcept {
  int error = 0;

  // we will go directly to the data-section block, so this must be already 
  // set! (aka, we should already have read the header)
  if (!data_section_pos) {
    fprintf(stderr,
            "[ERROR] Failed inspecting data for icgem file %s; did you forget "
            "to read the header? (traceback: %s)\n",
            filename.c_str(), __func__);
    return 1;
  }

  // open file ...
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening icgem file %s (traceback: %s)\n",
            filename.c_str(), __func__);
    return 1;
  }
  
  // go to the data section
  fin.seekg(data_section_pos);

  // reset all member variables characterizing the model
  max_degree_static_start=0;
  max_degree_static_stop=0;
  max_order_static_start=0;
  max_order_static_stop=0;
  max_degree_tv_start=0;
  max_degree_tv_stop=0;
  max_order_tv_start=0;
  max_order_tv_stop=0;
  harmonics.clear();
  
  char line[dso::Icgem::max_data_line];
  char *start, *end;

  // start reading ....
  fin.getline(line, max_data_line);
  while (fin.good()) {
    
    // gfc lines are for static-gravity field (if L=0=M, no effect ...)
    if (starts_with("gfc ", line)) {
      // expecting columns: degree, order, Clm, Slm, [...]; note that it
      // (seldom) happens that the doubles are written in fortran format ...
      start = line + 4;
      int ll = std::strtol(start, &end, 10);
      if (end == start) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      start = end;
      int mm = std::strtol(start, &end, 10);
      if (end == start) {
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      if (!max_degree_static_start && ll) max_degree_static_start = ll;
      if (ll > max_degree_static_stop) max_degree_static_stop = ll;
      if (!max_order_static_start && mm) max_order_static_start = mm;
      if (mm > max_order_static_stop) max_order_static_stop = mm;
    
    } else if (starts_with("gfct", line)) {
      // expecting columns: degree, order, Clm, Slm, [...]; note that it
      // (seldom) happens that the doubles are written in fortran format ...
      start = line + 4;
      int ll = std::strtol(start, &end, 10);
      if (end == start) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      start = end;
      int mm = std::strtol(start, &end, 10);
      if (end == start) {
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      if (!max_degree_tv_start && ll) max_degree_tv_start = ll;
      if (ll > max_degree_tv_stop) max_degree_tv_stop = ll;
      if (!max_order_tv_start && mm) max_order_tv_start = mm;
      if (mm > max_order_tv_stop) max_order_tv_stop = mm;
    
    } else if (starts_with("trnd", line)) {
      // expecting that this 'trnd' field should match the current degree and 
      // order of the --already read-- TVG coefficients
      start = line + 4;
      int ll = std::strtol(start, &end, 10);
      if (end == start) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      start = end;
      int mm = std::strtol(start, &end, 10);
      if (end == start) {
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      if ((ll != max_degree_tv_stop) || (mm != max_order_tv_stop)) {
        fprintf(stderr,
                "[ERROR] Reading line of type \'trnd\' but order/degree do not "
                "match with previous TVG coefficients read!\n");
        fprintf(stderr,
                "[ERROR] Current TVG degree and order: %d/%d, icgem file: %s "
                "(traceback: %s)\n",
                max_degree_tv_stop, max_order_tv_stop, filename.c_str(),
                __func__);
        return 1;
      }

    } else if (starts_with("acos", line) || starts_with("asin", line)) {
      // expecting that this 'acos'/'asin' field should match the current 
      // degree and order of the --already read-- TVG coefficients
      // example line:
      // acos   1    0  1.98940208316E-10  0.00000000000E+00 2.4920E-11 0.0000E+00 19500101.0000 19930115.0546 1.0
      start = line + 4;
      int ll = std::strtol(start, &end, 10);
      if (end == start) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      start = end;
      int mm = std::strtol(start, &end, 10);
      if (end == start) {
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      // parse values untill period
      double yperiod;
      for (int i=0; i<7; i++) {
        start = end;
        yperiod = std::strtod(start, &end);
        if (start == end) ++error;
      }
      if (error) {
        fprintf(stderr,
                "[ERROR] Failed parsing components in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      if ((ll != max_degree_tv_stop) || (mm != max_order_tv_stop)) {
        fprintf(stderr,
                "[ERROR] Reading line of type \'acos/asin\' but order/degree do not "
                "match with previous TVG coefficients read!\n");
        fprintf(stderr,
                "[ERROR] Current TVG degree and order: %d/%d, icgem file: %s "
                "(traceback: %s)\n",
                max_degree_tv_stop, max_order_tv_stop, filename.c_str(),
                __func__);
        return 1;
      }

      // the period is only allowd to not-exist, if this line describes the 
      // degree/order 1/0 coefficients
      bool period_exists = (std::find(harmonics.begin(), harmonics.end(),
                                      yperiod) == harmonics.end())
                               ? false
                               : true;
      if (ll == 1 && mm == 0) {
        if (!period_exists) harmonics.push_back(yperiod);
      } else {
        if (!period_exists) {
          fprintf(stderr, "[ERROR] Reading line of type \'acos/asin\' with unknown harmonic period!\n");
          fprintf(stderr, "[ERROR] Period %.3f/year not listed, line is \"%s\", in file %s (traceback: %s)\n", yperiod, line, filename.c_str(), __func__);
          return 1;
        }
      }
    } else {
      fprintf(stderr, "[WRNNG] ICGEM line skipped: \'%s\' (file %s)\n", line, filename.c_str());
    }
  }

    return 0;
}
