/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY / Pascal BAEHR
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/

#include <Helper/MathTools.h> //IDX
#include <cmath>
#include "AnglemapGeneration.h"

namespace AnglemapGeneration {
	/*
	* \brief Generates a theta angle when generating following recorded angle map
	* \param anglemapParams: angle map parameters
	* \param anglemap: recorded angle map and distributions precalculated at initialization. Not const so we can access memory address.
	* \param lookupValue: random number between 0 and 1
	* \return double theta value, integer lower index (referencing bin midpoints: 0->first bin midpoint, -1 -> first bin before midpoint, double overshoot (oveshoot: how many bins above previous bin midpoint (0..1)
	*/
	std::tuple<double, int, double>
		GenerateThetaFromAngleMap(const AnglemapParams& anglemapParams, Anglemap& anglemap, const double lookupValue) {
		int thetaLowerIndex; //can be -1 if lookupValue lower than first CDF value (theta below first bin midpoint)
		double theta, thetaOvershoot;

		if (lookupValue < anglemap.thetaLowerRatio) { //theta in lower theta region

			thetaLowerIndex = my_lower_bound(lookupValue,
				anglemap.theta_CDF_lower); //returns midpoint index after where lookup value resides in ( -1 .. thetaLowerLimit-2 )

			if (thetaLowerIndex == -1) { //theta in the first bin, below its midpoint

				thetaOvershoot = 0.5 * lookupValue / anglemap.theta_CDF_lower[0]; //between 0 and 0.5, 0 is theta=0, 0.5 is theta=first bin midpoint
				theta = GetTheta(thetaOvershoot, anglemapParams); //between 0 and the first section midpoint
			}
			else if (thetaLowerIndex == (anglemapParams.thetaLowerRes - 1)) { //theta in last half of lower part, above last recorded midpoitn CDF
				thetaOvershoot = 0.5 * (lookupValue - anglemap.theta_CDF_lower[thetaLowerIndex])
					/ (anglemap.thetaLowerRatio - anglemap.theta_CDF_lower[thetaLowerIndex]); //between 0 and 0.5, 0 is last bin midpoint, 0.5 is theta=thetaLimit
				theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot,
					anglemapParams); //between last midpoint and thetaLimit
			}
			else { //regular section
				if (anglemap.phi_pdfsums_lowerTheta[thetaLowerIndex] == anglemap.phi_pdfsums_lowerTheta[thetaLowerIndex + 1]) {
					//Theta probability is the same for the two adjacent theta lines, hence pdf slope is 0 between the two line midpoints. Simple linear interpolation between the two lines
					thetaOvershoot = (lookupValue - anglemap.theta_CDF_lower[thetaLowerIndex]) /
						(anglemap.theta_CDF_lower[thetaLowerIndex + 1] - anglemap.theta_CDF_lower[thetaLowerIndex]); //What fraction of a theta bin [0..1[ over the lower line's midpoint
					theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams);
				}
				else {
					//2nd degree interpolation.
					//To interpolate linearly in pdf, we have to interpolate on a parabolic curve connecting CDF sampled points
					//(Integral of a linear curve y=mx+b is a parabolic curve y=m/2*x^2 + bx + c )
					// y(x) = ax^2 + bx + c            where a=m/2
					// c = y(x0) = y0 = CDF value at lower line midpoint          //Motion equivalent: initial position (s0)
					// b: pdf value at lower line midpoint                        //Motion equivalent: initial speed (v0)
					// a: half of pdf slope at lower line midpoint, m/2           //Motion equivalent: acceleration (per two)  (a/2)
					// dy := y - c                                                //Motion equivalent: distance from start    (ds = s - s0)
					// dx := x - [x at lower line midpoint]                       //Motion equivalent: elapsed time           (dt = t - t0)
					// dy = a*dx^2 + bx                                           //Motion equivalent: ds = a/2 * t^2 + v0 * t
					// Solved for dx:  dx = ( -b + sqrt(b^2 +4*a*dy) ) / (2a)     //Solution with + is the right one in all cases
					
					double thetaStep = GetTheta((double)thetaLowerIndex + 1.5, anglemapParams) -
						GetTheta((double)thetaLowerIndex + 0.5, anglemapParams); //works in both lower and higher theta part
					double c = anglemap.theta_CDF_lower[thetaLowerIndex]; //CDF value at lower line midpoint
					double b = (double)anglemap.phi_pdfsums_lowerTheta[thetaLowerIndex] / (double)anglemap.theta_CDFsum_higher / thetaStep; //theta probability (pdf value) at lower line
					double a = 0.5 * ((double)(anglemap.phi_pdfsums_lowerTheta[thetaLowerIndex + 1]) -
						(double)anglemap.phi_pdfsums_lowerTheta[thetaLowerIndex]) /
						(double)anglemap.theta_CDFsum_higher / Sqr(thetaStep); //pdf slope between lower and next line midpoints
					double dy = lookupValue - c;

					double dx = (-b + sqrt(Sqr(b) + 4 * a * dy)) /
						(2 * a); //Since b>=0 it's the + branch of the +- that is valid for us

					thetaOvershoot = dx / thetaStep; //Fraction of bin over lower line midpoint
					theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams);
				}
			}

		}
		else { //theta in higher res region. For detailed comments, see lower part

			thetaLowerIndex = anglemapParams.thetaLowerRes + my_lower_bound(lookupValue,
				anglemap.theta_CDF_higher);

			if (thetaLowerIndex == anglemapParams.thetaLowerRes - 1) { //theta in the first half of the higher res part (below recorded CDF at midpoint)
				thetaOvershoot = 0.5 * (lookupValue - anglemap.thetaLowerRatio) / (anglemap.theta_CDF_higher[0] - anglemap.thetaLowerRatio); //between 0 and 0.5, 0 is theta=thetaLimit, 0.5 is theta=first bin midpoint in higher part
				theta = GetTheta((double)anglemapParams.thetaLowerRes + thetaOvershoot, anglemapParams); //between thetaLowerRes and thetaLowerRes+0.5
			}
			else if (thetaLowerIndex == (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1)) { //theta in last half of higher res part
				thetaOvershoot = 0.5 * (lookupValue - anglemap.theta_CDF_higher[thetaLowerIndex - anglemapParams.thetaLowerRes])
					/ (1.0 - anglemap.theta_CDF_higher[thetaLowerIndex - anglemapParams.thetaLowerRes]); //between 0 and 0.5
				theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot,
					anglemapParams); //between 0 and the first section end
			}
			else { //regular section
				if (anglemap.phi_pdfsums_higherTheta[thetaLowerIndex - anglemapParams.thetaLowerRes] == anglemap.phi_pdfsums_higherTheta[thetaLowerIndex - anglemapParams.thetaLowerRes + 1]) {
					//The pdf's slope is 0, linear interpolation
					thetaOvershoot = (lookupValue - anglemap.theta_CDF_higher[thetaLowerIndex - anglemapParams.thetaLowerRes]) /
						(anglemap.theta_CDF_higher[thetaLowerIndex - anglemapParams.thetaLowerRes + 1] - anglemap.theta_CDF_higher[thetaLowerIndex - anglemapParams.thetaLowerRes]);
					theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams);
				}
				else {
					double thetaStep = GetTheta((double)thetaLowerIndex + 1.5, anglemapParams) -
						GetTheta((double)thetaLowerIndex + 0.5, anglemapParams);
					double c = anglemap.theta_CDF_higher[thetaLowerIndex - anglemapParams.thetaLowerRes]; //CDF value at lower index
					double b = (double)anglemap.phi_pdfsums_higherTheta[thetaLowerIndex - anglemapParams.thetaLowerRes] / (double)anglemap.theta_CDFsum_higher /
						thetaStep; //pdf value at lower index
					double a = 0.5 * ((double)(anglemap.phi_pdfsums_higherTheta[thetaLowerIndex - anglemapParams.thetaLowerRes + 1]) -
						(double)anglemap.phi_pdfsums_higherTheta[thetaLowerIndex - anglemapParams.thetaLowerRes]) /
						(double)anglemap.theta_CDFsum_higher / Sqr(thetaStep); //pdf slope at lower index
					double dy = lookupValue - c;

					double dx = (-b + sqrt(Sqr(b) + 4 * a * dy)) /
						(2 * a); //Since b>=0 it's the + branch of the +- that is valid for us

					thetaOvershoot = dx / thetaStep;
					theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams);
				}
			}
		}
		return { theta, thetaLowerIndex, thetaOvershoot };
	}

	/**
	* \brief Generates phi angle (azimuth) from angle map
	* \param thetaLowerIndex lower index of theta bin (line)
	* \param thetaOvershoot what fraction of a bin are we over the midpoint of the line with thetaLowerIndex [0..1[ -> used to a weighing of the lower and higher theta lines
	* \param anglemapParams parameters of the angle map
	* \param angleMap the angle map itself, with precalculated quantities by InitializeAngleMap(). Cannot be passed as const since pass its memory addresses to my_lower_bound
	* \param randomGenerator reference to the random number generator (Mersenne Twister)
	* \return phi angle
	*/
	double GeneratePhiFromAngleMap(const int thetaLowerIndex, const double thetaOvershoot, const AnglemapParams& anglemapParams, Anglemap& anglemap, double lookupValue) {
		if (anglemapParams.phiWidth == 1) return -PI + 2.0 * PI * lookupValue; //simplest case speedup, uniform phi distribution
		int phiLowerIndex;
		//The lookupValue is looked up from first midpoint to "width+1"th midpoint. It will never be in the first half bin, but can go over width by half bin due to periodic BC
		//In other words, the lookupValue of 0..1 is mapped to 0.5...width+0.5 of phi_CDFs
		//Because of this, before calling my_lower_bound, we add the CDF value of the first bin's midpoint (first CDF element of actual line)
		double weigh; //0: take previous theta line (with index thetaLowerIndex), 1: take next theta line (index: thetaLowerIndex + 1), 0..1: interpolate in-between
		if (thetaLowerIndex < (int)anglemapParams.thetaLowerRes) { //In lower part
			if (thetaLowerIndex == -1) { //theta lower than first theta line midpoint, no weighing with next line
				lookupValue += anglemap.phi_CDFs_lowerTheta[0]; //periodic BCs over -PI...PI, can be larger than 1. So add value of first midpoint
				phiLowerIndex = my_lower_bound(lookupValue, &anglemap.phi_CDFs_lowerTheta[0], anglemapParams.phiWidth); //take the phi distro belonging to first theta with 100% weight
				weigh = thetaOvershoot; // [0.5,1[, will add thetaLowerIndex=-1 then 0.5 when evaluating thetaIndex -> thetaLowerIndex=[0, .5[
			}
			else if (thetaLowerIndex == ((int)anglemapParams.thetaLowerRes - 1)) { //theta above lower part's last bin midpoint, no weighing with previous line, neither with higher part first line (thetaStep is different between lower and higher parts)
				lookupValue += anglemap.phi_CDFs_lowerTheta[thetaLowerIndex * anglemapParams.phiWidth]; //(first element of last line) periodic BCs over -PI...PI, can be larger than 1, so add first half-bin midpoint
				phiLowerIndex = my_lower_bound(lookupValue, &anglemap.phi_CDFs_lowerTheta[thetaLowerIndex * anglemapParams.phiWidth], anglemapParams.phiWidth); //take the phi CDF belonging to last line with 100% weight
				weigh = thetaOvershoot; // [0 - 0.5[ above last theta bin midpoint, will add 0.5 when evaluating thetaIndex
			}
			else {
				//Here we do a weighing both by the hit sum of the previous and next lines (w1 and w2) and also the weighs of the two lines based on thetaOvershoot (w3 and w4)
				//This is required because a line with lot of hits should dominate the interpolation of the phi distributions
				// w1: sum of hits in previous line
				// w2: sum of hits in next line
				// w3: weigh of previous line (1 - thetaOvershoot)
				// w4: weigh of next line     (thetaOvershoot)
				// result: previous value weight: w1*w3 / (w1*w3 + w2*w4)
				//         next     value weight: w2*w4 / (w1*w3 + w2*w4) <- this will be the input for weighed_lower_bound

				double div;
				div = ((double)anglemap.phi_pdfsums_lowerTheta[thetaLowerIndex] * (1.0 - thetaOvershoot) + (double)anglemap.phi_pdfsums_lowerTheta[thetaLowerIndex + 1] * thetaOvershoot); // (w1*w3 + w2*w4)
				if (div > 0.0) {
					weigh = (thetaOvershoot * (double)anglemap.phi_pdfsums_lowerTheta[thetaLowerIndex + 1]) / div;    // w2*w4 / (w1*w3 + w2*w4)
					//weigh = thetaOvershoot; //debug
					//weigh=0.99;
				}
				else {
					weigh = thetaOvershoot;
				}
				lookupValue += Weigh((double)anglemap.phi_CDFs_lowerTheta[thetaLowerIndex * anglemapParams.phiWidth],
					(double)anglemap.phi_CDFs_lowerTheta[(thetaLowerIndex + 1) * anglemapParams.phiWidth], weigh); //Shift lookup value by weighed average of first elements of CDF lines
				phiLowerIndex = weighed_lower_bound_X(lookupValue, weigh,
					&anglemap.phi_CDFs_lowerTheta[thetaLowerIndex * anglemapParams.phiWidth],
					&anglemap.phi_CDFs_lowerTheta[(thetaLowerIndex + 1) * anglemapParams.phiWidth],
					anglemapParams.phiWidth);
			}
		}
		else { //In higher part, for detailed comments see lower part
			if (thetaLowerIndex == (int)anglemapParams.thetaLowerRes - 1) { //theta above thetaLimit but below first higher part theta bin midpoint
				lookupValue += anglemap.phi_CDFs_higherTheta[0]; //periodic BCs over -PI...PI, can be larger than 1, shift lookupValue by first element
				phiLowerIndex = my_lower_bound(lookupValue, &anglemap.phi_CDFs_higherTheta[0], anglemapParams.phiWidth); //take the phi distro belonging to first theta with full weigh
				weigh = thetaOvershoot; 
			}
			else if (thetaLowerIndex == ((anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 1)) { //last theta half section
				lookupValue += anglemap.phi_CDFs_higherTheta[(thetaLowerIndex - anglemapParams.thetaLowerRes) * anglemapParams.phiWidth]; //(first element of last line) periodic BCs over -PI...PI, can be larger than 1
				phiLowerIndex = my_lower_bound(lookupValue, &anglemap.phi_CDFs_higherTheta[(thetaLowerIndex - anglemapParams.thetaLowerRes) * anglemapParams.phiWidth], anglemapParams.phiWidth); //take entirely the phi ditro belonging to last line
				weigh = thetaOvershoot; 
			}
			else {
				//Here we do a weighing both by the hit sum of the previous and next lines (w1 and w2) and also the weighs of the two lines based on thetaOvershoot (w3 and w4)
				// w1: sum of hits in previous line
				// w2: sum of hits in next line
				// w3: weigh of previous line (1 - thetaOvershoot)
				// w4: weigh of next line     (thetaOvershoot)
				// result: previous value weight: w1*w3 / (w1*w3 + w2*w4)
				//         next     value weight: w2*w4 / (w1*w3 + w2*w4) <- this will be the input for weighed_lower_bound

				double div;
				div = ((double)anglemap.phi_pdfsums_higherTheta[thetaLowerIndex - anglemapParams.thetaLowerRes] * (1.0 - thetaOvershoot) + (double)anglemap.phi_pdfsums_higherTheta[thetaLowerIndex - anglemapParams.thetaLowerRes + 1] * thetaOvershoot); // (w1*w3 + w2*w4)
				if (div > 0.0) {
					weigh = (thetaOvershoot * (double)anglemap.phi_pdfsums_higherTheta[thetaLowerIndex - anglemapParams.thetaLowerRes + 1]) / div;    // w2*w4 / (w1*w3 + w2*w4)
				}
				else {
					weigh = thetaOvershoot;
				}
				lookupValue += Weigh((double)anglemap.phi_CDFs_higherTheta[(thetaLowerIndex - anglemapParams.thetaLowerRes) * anglemapParams.phiWidth],
					(double)anglemap.phi_CDFs_higherTheta[(thetaLowerIndex - anglemapParams.thetaLowerRes + 1) * anglemapParams.phiWidth], weigh);
				phiLowerIndex = weighed_lower_bound_X(lookupValue, weigh,
					&anglemap.phi_CDFs_higherTheta[(thetaLowerIndex - anglemapParams.thetaLowerRes) * anglemapParams.phiWidth],
					&anglemap.phi_CDFs_higherTheta[(thetaLowerIndex - anglemapParams.thetaLowerRes + 1) * anglemapParams.phiWidth],
					anglemapParams.phiWidth);
			}
		}

		//Weigh and phiLowerIndex calculated, get actual phi angle
		double phi, phiOvershoot;
		double thetaIndex = (double)thetaLowerIndex + 0.5 + weigh;
		if (phiLowerIndex == -1) { //first half section
			DEBUG_BREAK; //should never happen since we shifted the lookup value with first phi bin CDF midpoint
			phiOvershoot = 0.5 * lookupValue / GetPhiCDFValue(thetaIndex, 0, anglemapParams, anglemap); //between 0 and .5
			phi = GetPhi(phiOvershoot, anglemapParams); //between 0 and the first bin midpoint
		}
		else { //regular or last section
			if (GetPhiNormalizedPdfValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap) ==
				GetPhiNormalizedPdfValue(thetaIndex, phiLowerIndex + 1, anglemapParams, anglemap)) {
				//The pdf's slope is 0 between the two phi bin midpoints, linear interpolation
				phiOvershoot = (lookupValue - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap))
					/ (GetPhiCDFValue(thetaIndex, phiLowerIndex + 1, anglemapParams, anglemap) - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap));
				phi = GetPhi((double)phiLowerIndex + 0.5 + phiOvershoot, anglemapParams);
			}
			else {

				//2nd degree interpolation, see GenerateThetaFromAngleMap comments for details
				double phiStep = 2.0 * PI / (double)anglemapParams.phiWidth;
				double c = GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap); //CDF value at lower index midpoint
				double b = GetPhiNormalizedPdfValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap) / phiStep; //pdf value at lower index midpoint
				double a = 0.5 * (GetPhiNormalizedPdfValue(thetaIndex, phiLowerIndex + 1, anglemapParams, anglemap) -
					GetPhiNormalizedPdfValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap)) / Sqr(phiStep); //pdf slope between lower and higher lines (and thus their midpoints)
				double dy = lookupValue - c;

				double D = Sqr(b) + 4 * a * dy; //Discriminant. In rare cases it might be slightly negative, then fall back to linear interpolation:
				if (D < 0) {
					phiOvershoot = (lookupValue - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap))
						/ (GetPhiCDFValue(thetaIndex, (int)IDX(phiLowerIndex + 1, anglemapParams.phiWidth), anglemapParams, anglemap)
							- GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap));
				}
				else {
					double dx = (-b + sqrt(Sqr(b) + 4 * a * dy)) /
						(2 * a); //Since b>=0 it's the + branch of the +- that is valid for us
					phiOvershoot = dx / phiStep;
				}
				phi = GetPhi((double)phiLowerIndex + 0.5 + phiOvershoot, anglemapParams);
			}
		}
		assert(phi > -PI && phi < PI);
		return phi;
	}

	/**
	* \brief Maps anglemap theta index to theta value
	* \param thetaIndex theta index from 0 to map theta size (lowerRes+higherRes), can be non-integer
	* \param anglemapParams parameters of the angle map
	* \return theta angle in rad
	*/
	double GetTheta(const double thetaIndex, const AnglemapParams& anglemapParams) {
		if ((size_t)(thetaIndex) < anglemapParams.thetaLowerRes) { // 0 < theta < limit
			return anglemapParams.thetaLimit * (thetaIndex) / (double)anglemapParams.thetaLowerRes;
		}
		else { // limit < theta < PI/2
			return anglemapParams.thetaLimit +
				(PI / 2.0 - anglemapParams.thetaLimit) * (thetaIndex - (double)anglemapParams.thetaLowerRes) /
				(double)anglemapParams.thetaHigherRes;
		}
	}

	/**
	* \brief makes phiIndex circular and converts from index to -pi...pi
	* \param phiIndex phi index
	* \param anglemapParams parameters of the angle map
	* \return phi angle
	*/

	double GetPhi(const double phiIndex, const AnglemapParams& anglemapParams)
		//makes phiIndex circular and converts from index to -pi...pi
	{
		double width = (double)anglemapParams.phiWidth;
		double correctedIndex = (phiIndex < width) ? phiIndex : phiIndex - width;
		return -PI + 2.0 * PI * correctedIndex / width;
	}

	double GetPhiNormalizedPdfValue(const double thetaIndex, const int phiLowerIndex,
		const AnglemapParams& anglemapParams, const Anglemap& anglemap)
	{
		size_t phiIndex = IDX(phiLowerIndex,anglemapParams.phiWidth);//phiLowerIndex is circularized
		if (thetaIndex < 0.5) { //First line
			return (anglemapParams.thetaLowerRes > 0) ? anglemap.phi_pdfs_lowerTheta[phiIndex] : anglemap.phi_pdfs_higherTheta[phiIndex];
		}
		else if (thetaIndex > (double)(anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) { //last line
			return (anglemapParams.thetaHigherRes > 0) ?
				anglemap.phi_pdfs_higherTheta[anglemapParams.phiWidth * (anglemapParams.thetaHigherRes - 1) + phiIndex]
				: anglemap.phi_pdfs_lowerTheta[anglemapParams.phiWidth * (anglemapParams.thetaLowerRes - 1) + phiIndex];
		}
		else { //intermediate line
			size_t thetaLowerIndex = (size_t)(thetaIndex - 0.5);
			double thetaOvershoot = thetaIndex - 0.5 - (double)thetaLowerIndex;

			double valueFromLowerpdf = (thetaLowerIndex < anglemapParams.thetaLowerRes)
				? anglemap.phi_pdfs_lowerTheta[thetaLowerIndex*anglemapParams.phiWidth + phiIndex]
				: anglemap.phi_pdfs_higherTheta[(thetaLowerIndex-anglemapParams.thetaLowerRes)* anglemapParams.phiWidth + phiIndex];
			size_t nextLineIndex = thetaLowerIndex + 1;
			double valueFromHigherpdf = (nextLineIndex < anglemapParams.thetaLowerRes)
				? anglemap.phi_pdfs_lowerTheta[nextLineIndex *anglemapParams.phiWidth + phiIndex]
				: anglemap.phi_pdfs_higherTheta[(nextLineIndex -anglemapParams.thetaLowerRes)* anglemapParams.phiWidth + phiIndex];

			return Weigh(valueFromLowerpdf, valueFromHigherpdf, thetaOvershoot);
		}
	}

	/**
	* \brief Get phi value from cummulative density function
	* \param thetaIndex theta index
	* \param phiLowerIndex lower index of the bin
	* \param anglemapParams parameters of the angle map
	* \return phi cdf value
	*/

	double GetPhiCDFValue(const double thetaIndex, const int phiLowerIndex, const AnglemapParams& anglemapParams, const Anglemap& anglemap) {
		//todo: check if circular mapping ever happens (phiLowerIndex >= anglemapParams.phiWidth)
		if (thetaIndex < (double)anglemapParams.thetaLowerRes) { //In lower part
			if (thetaIndex < 0.5) {
				return (phiLowerIndex < anglemapParams.phiWidth)
					? anglemap.phi_CDFs_lowerTheta[phiLowerIndex] : 1.0 + anglemap.phi_CDFs_lowerTheta[0]; //phi is circular -> last slice mapped to first
			}
			else if (thetaIndex > (double)anglemapParams.thetaLowerRes - 0.5) {
				return (phiLowerIndex < anglemapParams.phiWidth)
					? anglemap.phi_CDFs_lowerTheta[anglemapParams.phiWidth * (anglemapParams.thetaLowerRes - 1) + phiLowerIndex]
					: 1.0 + anglemap.phi_CDFs_lowerTheta[anglemapParams.phiWidth * (anglemapParams.thetaLowerRes - 1)]; //Clamp CDF to 1.0
			}
			else {
				size_t thetaLowerIndex = (size_t)(thetaIndex - 0.5);
				double thetaOvershoot = thetaIndex - 0.5 - (double)thetaLowerIndex;
				double valueFromLowerCDF = (phiLowerIndex < anglemapParams.phiWidth)
					? anglemap.phi_CDFs_lowerTheta[anglemapParams.phiWidth * thetaLowerIndex + phiLowerIndex]
					: 1.0 + anglemap.phi_CDFs_lowerTheta[anglemapParams.phiWidth * (thetaLowerIndex)]; //Clamp CDF to 1.0
				double valueFromHigherCDF = (phiLowerIndex < anglemapParams.phiWidth)
					? anglemap.phi_CDFs_lowerTheta[anglemapParams.phiWidth * (thetaLowerIndex + 1) + phiLowerIndex]
					: 1.0 + anglemap.phi_CDFs_lowerTheta[anglemapParams.phiWidth * (thetaLowerIndex + 1)]; //Clamp CDF to 1.0
				return Weigh(valueFromLowerCDF, valueFromHigherCDF, thetaOvershoot);
			}
		}
		else { //In higher part
			if (thetaIndex < (double)anglemapParams.thetaLowerRes + 0.5) {
				return (phiLowerIndex < anglemapParams.phiWidth)
					? anglemap.phi_CDFs_higherTheta[phiLowerIndex] : 1.0 + anglemap.phi_CDFs_higherTheta[0];
			}
			else if (thetaIndex > (double) (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) {
				return (phiLowerIndex < anglemapParams.phiWidth)
					? anglemap.phi_CDFs_higherTheta[anglemapParams.phiWidth * (anglemapParams.thetaHigherRes - 1) + phiLowerIndex]
					: 1.0 + anglemap.phi_CDFs_higherTheta[anglemapParams.phiWidth * (anglemapParams.thetaHigherRes - 1)];
			}
			else {
				size_t thetaLowerIndex = (size_t)(thetaIndex - 0.5);
				double thetaOvershoot = thetaIndex - 0.5 - (double)thetaLowerIndex;
				double valueFromLowerCDF = (phiLowerIndex < anglemapParams.phiWidth)
					? anglemap.phi_CDFs_higherTheta[anglemapParams.phiWidth * (thetaLowerIndex - anglemapParams.thetaLowerRes) + phiLowerIndex]
					: 1.0 + anglemap.phi_CDFs_higherTheta[anglemapParams.phiWidth * (thetaLowerIndex - anglemapParams.thetaLowerRes)];
				double valueFromHigherCDF = (phiLowerIndex < anglemapParams.phiWidth)
					? anglemap.phi_CDFs_higherTheta[anglemapParams.phiWidth * (thetaLowerIndex - anglemapParams.thetaLowerRes + 1) + phiLowerIndex]
					: 1.0 + anglemap.phi_CDFs_higherTheta[anglemapParams.phiWidth * (thetaLowerIndex - anglemapParams.thetaLowerRes + 1)];
				return Weigh(valueFromLowerCDF, valueFromHigherCDF, thetaOvershoot);
			}
		}
	}

	/**
	* \brief Get phi value from cumulative density function
	* \param thetaIndex theta index
	* \param anglemapParams parameters of the angle map
	* \return phi cdf summed value
	*/

	double GetPhiCDFSum(const double thetaIndex, const AnglemapParams& anglemapParams,
		const Anglemap& anglemap) {
		if (thetaIndex < (double)anglemapParams.thetaLowerRes) { //In lower part

			if (thetaIndex < 0.5) { //first half-bin lower part
				return (double)anglemap.phi_pdfsums_lowerTheta[0];
			}
			else if (thetaIndex > (double)anglemapParams.thetaLowerRes - 0.5) { //last half-bin lower part
				return (double)anglemap.phi_pdfsums_lowerTheta[anglemapParams.thetaLowerRes - 1];
			}
			else { //regular bin in lower part
				size_t thetaLowerIndex = (size_t)(thetaIndex - 0.5);
				double thetaOvershoot = thetaIndex - 0.5 - (double)thetaLowerIndex;
				double valueFromLowerSum = (double)anglemap.phi_pdfsums_lowerTheta[thetaLowerIndex];
				double valueFromHigherSum = (double)anglemap.phi_pdfsums_lowerTheta[thetaLowerIndex + 1];
				return Weigh(valueFromLowerSum, valueFromHigherSum, thetaOvershoot);
			}
		}
		else { //In higher part
			if (thetaIndex < (double)anglemapParams.thetaLowerRes + 0.5) { //first half-bin higher part
				return (double)anglemap.phi_pdfsums_higherTheta[0];
			}
			else if (thetaIndex > (double) (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) { //last half-bin higher part
				return (double)anglemap.phi_pdfsums_higherTheta[anglemapParams.thetaHigherRes - 1];
			}
			else { //regular bin in higher part
				size_t thetaLowerIndex = (size_t)(thetaIndex - 0.5);
				double thetaOvershoot = thetaIndex - 0.5 - (double)thetaLowerIndex;
				double valueFromLowerSum = (double)anglemap.phi_pdfsums_higherTheta[thetaLowerIndex - anglemapParams.thetaLowerRes];
				double valueFromHigherSum = (double)anglemap.phi_pdfsums_higherTheta[thetaLowerIndex - anglemapParams.thetaLowerRes + 1];
				return Weigh(valueFromLowerSum, valueFromHigherSum, thetaOvershoot);
			}
		}
	}
} // namespace
