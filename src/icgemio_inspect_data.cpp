#include "icgemio.hpp"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <vector>
#if defined(__GNUC__) && (__GNUC__ >= 11)
#  include <charconv>
#else
#  include <cstdlib>
#  include <cerrno>
#endif 

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

const char *skip_ws(const char *str) noexcept {
  while (*str && *str==' ') ++str;
  return str;
}
} // namespace

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
  degree_static_start = 0;
  degree_static_stop = 0;
  order_static_start = 0;
  order_static_stop = 0;
  degree_tv_start = 0;
  degree_tv_stop = 0;
  order_tv_start = 0;
  order_tv_stop = 0;
  // current degree(ll) and order(mm) of TVG (while parsing)
  int tvg_ll(-1), tvg_mm(-1);
  // clear out any harmonics
  harmonics.clear();

  char line[dso::Icgem::max_data_line];
  char *start;
  [[maybe_unused]]char *end;

  #if defined(__GNUC__) && (__GNUC__ < 11)
  // clear errno before we start reading ...
  errno = 0;
  #endif

  // start reading ....
  fin.getline(line, max_data_line);
  while (fin.good()) {

    // size of current line
    const int sz = std::strlen(line);

    // gfc lines are for static-gravity field (if L=0=M, no effect ...)
    if (starts_with("gfc ", line)) {
      // expecting columns: degree, order, Clm, Slm, [...]; note that it
      // (seldom) happens that the doubles are written in fortran format ...
      start = line + 4;
      #if defined(__GNUC__) && (__GNUC__ < 11)
      int ll = std::strtol(start, &end, 10);
      if ((end == start) || errno) {
      #else
      int ll;
      auto ccres = std::from_chars(skip_ws(start), line+sz, ll);
      if (ccres.ec != std::errc{}) {
      #endif
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      #if defined(__GNUC__) && (__GNUC__ < 11)
      start = end;
      int mm = std::strtol(start, &end, 10);
      if (end == start) {
      #else
      int mm;
      ccres = std::from_chars(skip_ws(ccres.ptr), line+sz, mm);
      if (ccres.ec != std::errc{}) {
      #endif
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      // check if we have a max static degree (start/stop)
      if (!degree_static_start && ll)
        degree_static_start = ll;
      if (ll > degree_static_stop)
        degree_static_stop = ll;
      // check if we have a max static order (start/stop)
      if (!order_static_start && mm)
        order_static_start = mm;
      if (mm > order_static_stop)
        order_static_stop = mm;

    // gfct lines are for tvg field (if L=0=M, no effect ...)
    } else if (starts_with("gfct", line)) {
      // expecting columns: degree, order, Clm, Slm, [...]; note that time is
      // of no interest here
      start = line + 4;
      #if defined(__GNUC__) && (__GNUC__ < 11)
      int ll = std::strtol(start, &end, 10);
      if ((end == start) || errno) {
      #else
      int ll;
      auto ccres = std::from_chars(skip_ws(start), line+sz, ll);
      if (ccres.ec != std::errc{}) {
      #endif
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      #if defined(__GNUC__) && (__GNUC__ < 11)
      start = end;
      int mm = std::strtol(start, &end, 10);
      if (end == start) {
      #else
      int mm;
      ccres = std::from_chars(skip_ws(ccres.ptr), line+sz, mm);
      if (ccres.ec != std::errc{}) {
      #endif
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      // assign current ll/mm of tvg
      tvg_ll = ll;
      tvg_mm = mm;

      // check if we have a max tvg degree (start/stop)
      if (!degree_tv_start && ll)
        degree_tv_start = ll;
      if (ll > degree_tv_stop)
        degree_tv_stop = ll;
      // check if we have a max tvg order (start/stop)
      if (!order_tv_start && mm)
        order_tv_start = mm;
      if (mm > order_tv_stop)
        order_tv_stop = mm;

    // trnd lines are for trend/drift (if L=0=M, no effect ...)
    } else if (starts_with("trnd", line)) {
      // expecting that this 'trnd' field should match the current degree and
      // order of the --already read-- TVG coefficients
      start = line + 4;
      #if defined(__GNUC__) && (__GNUC__ < 11)
      int ll = std::strtol(start, &end, 10);
      if ((end == start) || errno) {
      #else
      int ll;
      auto ccres = std::from_chars(skip_ws(start), line+sz, ll);
      if (ccres.ec != std::errc{}) {
      #endif
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      #if defined(__GNUC__) && (__GNUC__ < 11)
      start = end;
      int mm = std::strtol(start, &end, 10);
      if (end == start) {
      #else
      int mm;
      ccres = std::from_chars(skip_ws(ccres.ptr), line+sz, mm);
      if (ccres.ec != std::errc{}) {
      #endif
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      if ((ll != tvg_ll) || (mm != tvg_mm)) {
        fprintf(stderr,
                "[ERROR] Reading line of type \'trnd\' but order/degree do not "
                "match with previous TVG coefficients read (%d,%d)!\n",
                ll, mm);
        fprintf(stderr,
                "[ERROR] Current TVG degree and order: %d/%d, icgem file: %s "
                "(traceback: %s)\n",
                tvg_ll, tvg_mm, filename.c_str(), __func__);
        return 1;
      }

    // acos/asin are for periodic terms
    } else if (starts_with("acos", line) || starts_with("asin", line)) {
      // expecting that this 'acos'/'asin' field should match the current
      // degree and order of the --already read-- TVG coefficients
      // example line:
      // acos   1    0  1.98940208316E-10  0.00000000000E+00 2.4920E-11
      // 0.0000E+00 19500101.0000 19930115.0546 1.0
      #if defined(__GNUC__) && (__GNUC__ < 11)
      int ll = std::strtol(start, &end, 10);
      if ((end == start) || errno) {
      #else
      int ll;
      auto ccres = std::from_chars(skip_ws(start), line+sz, ll);
      if (ccres.ec != std::errc{}) {
      #endif
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      #if defined(__GNUC__) && (__GNUC__ < 11)
      start = end;
      int mm = std::strtol(start, &end, 10);
      if (end == start) {
      #else
      int mm;
      ccres = std::from_chars(skip_ws(ccres.ptr), line+sz, mm);
      if (ccres.ec != std::errc{}) {
      #endif
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      // parse values untill period
      double yperiod;
      for (int i = 0; i < 7; i++) {
        #if defined(__GNUC__) && (__GNUC__ < 11)
        start = end;
        yperiod = std::strtod(start, &end);
        if ((start == end) || errno)
          ++error;
        #else
        ccres = std::from_chars(skip_ws(ccres.ptr), line+sz, yperiod);
        error += (ccres.ec != std::errc{});
        #endif
      }
      if (error) {
        fprintf(stderr,
                "[ERROR] Failed parsing components in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      if ((ll != tvg_ll) || (mm != tvg_mm)) {
        fprintf(stderr, "[ERROR] Reading line of type \'acos/asin\' but "
                        "order/degree do not "
                        "match with previous TVG coefficients read!\n");
        fprintf(stderr,
                "[ERROR] Current TVG degree and order: %d/%d, icgem file: %s "
                "(traceback: %s)\n",
                tvg_ll, tvg_mm, filename.c_str(),
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
        if (!period_exists)
          harmonics.push_back(yperiod);
      } else {
        if (!period_exists) {
          fprintf(stderr, "[ERROR] Reading line of type \'acos/asin\' with "
                          "unknown harmonic period!\n");
          fprintf(stderr,
                  "[ERROR] Period %.3f/year not listed, line is \"%s\", in "
                  "file %s (traceback: %s)\n",
                  yperiod, line, filename.c_str(), __func__);
          return 1;
        }
      }
    } else {
      fprintf(stderr, "[WRNNG] ICGEM line skipped: \'%s\' (file %s)\n", line,
              filename.c_str());
    }

    // read next line ...
    fin.getline(line, max_data_line);
  }

  if (!fin.eof()) {
    fprintf(
        stderr,
        "[ERROR] Failed to reach EOF, parsing error! file %s (traceback: %s)\n",
        filename.c_str(), __func__);
    return 5;
  }

  return 0;
}
