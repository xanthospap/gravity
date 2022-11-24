#ifndef __ICGEM_POTENTIAL_IO_HPP__
#define __ICGEM_POTENTIAL_IO_HPP__

/// Basic handling of icgem files (holding gravity potential models).
/// Information can be found on the International Centre for Global Earth 
/// Models (ICGEM) website, http://icgem.gfz-potsdam.de/home
/// see Ince, E. S., Barthelmes, F., Reißland, S., Elger, K., Förste, C.,
/// Flechtner, F., Schuh, H. (2019): ICGEM – 15 years of successful collection 
/// and distribution of global gravitational models, associated services and 
/// future plans.- Earth System Science Data, 11, pp.647 - 674, DOI : 
/// http:// doi.org/10.5194/essd-11-647-2019.

#include <fstream>
#include <cstring>
#include <vector>
#include "harmonic_coeffs.hpp"

namespace dso {

/// @brief A class to hold the reading/parsing of ICGEM gravity models.
/// Download a gfc file from http://icgem.gfz-potsdam.de/tom_longtime and parse
/// it via this class. Note that the implementation is still incomplete and 
/// not all models/parameters are read (aka only parameters of type 'gfc' can be
/// parsed but some models also have more coefficients).
class Icgem {
public:
  typedef std::ifstream::pos_type pos_type;

  ///< approximate max data line length, see
  ///< http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf
  static constexpr const std::size_t max_data_line = 512;

private:
    std::string filename;

    ///< Position in stream where data starts
    pos_type data_section_pos{0};
    
    ///< Header information
    char product_type[64] = {'\0'};
    char modelname[128] = {'\0'};
    char tide_system[64] = "unknown";
    char norm[64] = "fully_normalized";
    char errors[64] = {'\0'};
    double earth_gravity_constant{0e0};
    double radius{0e0};
    int max_degree{0};

public:  
    ///< model/data information
    int degree_static_start{0};
    int degree_static_stop{0};
    int order_static_start{0};
    int order_static_stop{0};
    int degree_tv_start{0};
    int degree_tv_stop{0};
    int order_tv_start{0};
    int order_tv_stop{0};
    std::vector<double> harmonics;

    Icgem(const char *fn) noexcept : filename(fn) {};
    Icgem(const Icgem&) noexcept = delete;
    Icgem& operator=(const Icgem&) noexcept = delete;
    ~Icgem() noexcept {};

    /// @brief Return max degree of the file
    int degree() const noexcept { return std::max(degree_static_stop,degree_tv_stop); }
    /// @brief Return max order of the file
    int order() const noexcept {return std::max(order_static_stop,order_tv_stop);}
    double earth_radius() const noexcept { return radius; }
    double gm() const noexcept { return earth_gravity_constant; }
  
    bool is_normalized() const noexcept {
      return !std::strcmp(norm, "fully_normalized");
    }

    #ifdef DEBUG
    void print_details();
    #endif

    /// @brief Read the file header and assign basic information.
    int parse_header() noexcept;

    /// @brief Inspect an ICGEM-format file; parse the file and extract info
    /// @note The header of the file should have already been read through
    ///       (see parse_header())
    /// @return Anything other than 0 denotes an error
    int inspect_data() noexcept;

    /// @brief Parse harmonic coefficients up to degree l and order m.
    /// Note that only data/values with a key value of 'gfc' are read. Some
    /// models however, include parameters with more keys.
    /// The function will also assign the model's constants (aka GM and Re) as
    /// well as the normalization status of the coefficients.
    ///
    /// @param[in] l Max degree of S/C harmonic coeffients values to read and
    ///            store.
    /// @param[in] k Max order of S/C harmonic coeffients values to read and
    ///            store (k <= l)
    /// @param[out] coeffs Pointer an instance of type HarmonicCoeffs where the
    ///            S/C harmonic coefficients are to be stored. Note that this
    ///            instance should have been allocated with enough space.
    /// @see http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf
    int parse_data(int l, int k, HarmonicCoeffs *coeffs) noexcept;
}; //Icgem

} // dso

#endif
