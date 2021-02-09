//
// Created by pascal on 2/9/21.
//

#include "CDFGeneration.h"
#include <cmath>

constexpr size_t cdf_size = 100; // points in a cumulative distribution function

namespace CDFGeneration {
/**
 * \brief Get ID (if it exists) of the Commulative Distribution Function (CFD)
 * for a particular temperature (bin) \param temperature temperature for the CFD
 * \return ID of the CFD
 */
    int GetCDFId(const std::set<double> &temperatureList, double temperature) {
        if (!temperatureList.empty()) {
            auto lowerBound = std::lower_bound(temperatureList.begin(),
                                               temperatureList.end(), temperature);
            if (lowerBound == temperatureList.begin())
                return -1;
            --lowerBound; // even temperatureList.end() can be a bound

            if (std::abs(temperature - *lowerBound) > 1E-5) {
                return std::distance(temperatureList.begin(), lowerBound);
            }
        }
        return -1;
    }

/**
 * \brief Generate a new Commulative Distribution Function (CFD) for a
 * particular temperature (bin) \param temperature for the CFD \return Previous
 * size of temperatures vector, which determines new ID
 */
    std::pair<int, std::vector<CDF_p>>
    GenerateNewCDF(std::set<double> &temperatureList, const double temperature,
                   const double gasMass) {
        size_t i = temperatureList.size();
        temperatureList.emplace(temperature);
        std::vector<CDF_p> cdf_v = Generate_CDF(temperature, gasMass, cdf_size);
        return std::make_pair((int) i, cdf_v);
    }

/**
 * \brief Generate cumulative distribution function (CFD) for the velocity
 * \param gasTempKelvins gas temperature in Kelvin
 * \param gasMassGramsPerMol molar gas mass in grams per mol
 * \param size amount of points/bins of the CFD
 * \return CFD as a Vector containing a pair of double values (x value =
 * speed_bin, y value = cumulated value)
 */
    std::vector<std::pair<double, double>>
    Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol, size_t size) {
        std::vector<std::pair<double, double>> cdf;
        cdf.reserve(size);
        constexpr double Kb = 1.38E-23;
        constexpr double R = 8.3144621;
        const double a = std::sqrt(
                Kb * gasTempKelvins /
                (gasMassGramsPerMol * 1.67E-27)); // distribution a parameter. Converting
        // molar mass to atomic mass

        // Generate cumulative distribution function
        double mostProbableSpeed =
                std::sqrt(2.0 * R * gasTempKelvins / (gasMassGramsPerMol / 1000.0));
        double binSize =
                4.0 * mostProbableSpeed /
                (double) size; // distribution generated between 0 and 4*V_prob

        for (size_t i = 0; i < size; i++) {
            double x = (double) i * binSize;
            double x_square_per_2_a_square =
                    std::pow(x, 2.0) / (2.0 * std::pow(a, 2.0));
            cdf.emplace_back(
                    std::make_pair(x, 1.0 - std::exp(-x_square_per_2_a_square) *
                                            (x_square_per_2_a_square + 1.0)));
        }

        return cdf;
    }
} // namespace CDFGeneration