#include <cmath>
#include "CDFGeneration.h"
#include <cmath>

constexpr size_t velocity_cdf_size = 100; // points in a cumulative distribution function

namespace CDFGeneration {
/**
 * \brief Get ID (if it exists) of the Commulative Distribution Function (CFD)
 * for a particular temperature (bin) \param temperature temperature for the CFD
 * \return ID of the CFD
 */

    int GetCDFId(const std::vector<double> &temperatureList, double temperature) {
        if(!temperatureList.empty()) {
            // find temp within error range of 1e-5
            //auto is_almost_equal = [temperature](int i){ return (std::abs(i - temperature) < 1E-5); };
            auto pos = std::find_if(temperatureList.begin(), temperatureList.end(), [temperature](double i){ return (std::abs(i - temperature) < 1E-5);});

            /*auto pos2 = std::find_if(temperatureList.begin(), temperatureList.end(), temperature,
                                                 [](double i,double j) {return (std::abs(i - j) < 1E-5);}
                                                 );*/
            if(pos == temperatureList.end())
                return -1;
            else
                return std::distance(temperatureList.begin(), pos);
        }
        return -1;
    }

/**
 * \brief Generate a new Commulative Distribution Function (CFD) for a
 * particular temperature (bin)
 * \param temperature for the CFD
 * \return Previous
 * size of temperatures vector, which determines new ID
 */
    std::pair<int, std::vector<IntegratedVelocityEntry>>
    GenerateNewCDF(std::vector<double> &temperatureList, double temperature,
                   double gasMass) {
        size_t i = temperatureList.size();
        temperatureList.push_back(temperature);
        auto cdf_vect = Generate_CDF(temperature, gasMass);
        return std::make_pair((int) i, cdf_vect);
    }

/**
 * \brief Generate cumulative distribution function (CFD) for the velocity
 * \param gasTempKelvins gas temperature in Kelvin
 * \param gasMassGramsPerMol molar gas mass in grams per mol
 * \param size amount of points/bins of the CFD
 * \return CFD as a Vector containing a pair of double values (x value =
 * speed_bin, y value = cumulated value)
 */
    std::vector<IntegratedVelocityEntry>
    Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol) {
        std::vector<IntegratedVelocityEntry> cdf;
        cdf.reserve(velocity_cdf_size);
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
                (double) velocity_cdf_size; // distribution generated between 0 and 4*V_prob

        for (size_t i = 0; i < velocity_cdf_size; i++) {
            double x = (double) i * binSize;
            double x_square_per_2_a_square =
                    std::pow(x, 2.0) / (2.0 * std::pow(a, 2.0));
            IntegratedVelocityEntry entry;
            entry.velocity=x;
            entry.cumulativeProbability = 1.0 - std::exp(-x_square_per_2_a_square) * (x_square_per_2_a_square + 1.0);
            cdf.emplace_back(entry);
        }

        return cdf;
    }
} // namespace CDFGeneration