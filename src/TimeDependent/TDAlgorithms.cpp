//
// Created by pascal on 7/9/22.
//

#include "TDAlgorithms.h"
#include "fmt/color.h"
#include "fmt/core.h"
#include <cmath>

namespace MFTD {

/*!
 * @brief Lookup the index of the interval related to a given key and a start
 * position for accelerated lookup
 * @param key specific moment
 * @param moments vector of time intervals
 * @param startIndex offset to only look in a subset of moments
 * @return -1 if moment doesnt relate to an interval, else index of moment (+1
 * to account for [0]== steady state)
 * TODO: Removed +1 for test runs
 */
    int LookupMomentIndex(double key,
                          const std::vector<std::pair<double, double>> &moments,
                          size_t startIndex) {

        if (!moments.empty()) {
            auto lowerBound = std::lower_bound(moments.begin() + startIndex,
                                               moments.end(), std::make_pair(key, key));
            if (lowerBound == moments.begin())
                return -1;
            --lowerBound; // even moments.end() can be a bound

            if (lowerBound->first <= key && key <= lowerBound->second) {
                return static_cast<int>(
                        std::distance(moments.begin() + startIndex, lowerBound) + startIndex);
            }
        }
        return -1;
    }

    int LookupMomentIndex(double key, const std::vector<Moment> &moments) {

        if (!moments.empty()) {
            auto lowerBound = std::lower_bound(moments.begin(), moments.end(),
                                               std::make_pair(key, key));
            if (lowerBound == moments.begin())
                return -1;
            --lowerBound; // even moments.end() can be a bound

            if (lowerBound->first <= key && key <= lowerBound->second) {
                return static_cast<int>(std::distance(moments.begin(), lowerBound));
            }
        }
        return -1;
    }

    int quadraticSearch(double key, const std::vector<Moment> &moments) {

        if (!moments.empty()) {
            int p1 = 0;
            int p2 = 0;
            int mid = 0;
            int first = 0;
            int last = moments.size() - 1;

            while (first <= last) {
                mid = (first + last) * 0.5;
                p1 = first + (last - first) * 0.25;
                p2 = first + (last - first) * 0.75;

                // TODO: Out of the box needs 6 comparisons to check for a hit, should be
                // reducable by nesting
                /*if(key >= moments[p1].first){ // found
                    if(key <= moments[p1].second)
                        return p1;
                    else if(key >= moments[mid].first){ // found
                        if(key <= moments[mid].second)
                            return mid;
                        else if(key >= moments[p2].first && key <= moments[p2].second){ //
                found return p2;
                        }
                    }
                }*/
                if (p1 > mid || p2 < mid) {
                    fmt::print("Indices messed up {} < {} < {}\n", p1, mid, p2);
                }
                if (key >= moments[mid].first && key <= moments[mid].second) { // found
                    return mid;
                } else if (key >= moments[p1].first &&
                           key <= moments[p1].second) { // found
                    return p1;
                } else if (key >= moments[p2].first &&
                           key <= moments[p2].second) { // found
                    return p2;
                } else if (key < moments[mid].second && key < moments[p1].first) {
                    last = p1 - 1;
                } else if (key < moments[mid].second && key > moments[p1].first) {
                    first = p1 + 1;
                    last = mid - 1;
                } else if (key > moments[mid].first && key > moments[p2].second) {
                    first = p2 + 1;
                } else if (key > moments[mid].first && key < moments[p2].second) {
                    first = mid + 1;
                    last = p2 - 1;
                }
            }
            return -1;
        }
        return -1;
    }

    int interpolationSearch(double key, const std::vector<Moment> &moments) {

        // Find indexes of two corners
        int lo = 0, hi = (moments.size() - 1);

        // Since array is sorted, an element present
        // in array must be in range defined by corner
        while (lo <= hi && key >= moments[lo].first && key <= moments[hi].second) {
            if (lo == hi) {
                if (moments[lo].first <= key && key <= moments[lo].second)
                    return lo;
                return -1;
            }
            // Probing the position with keeping
            // uniform distribution in mind.
            int pos =
                    lo + (((double) (hi - lo) / (moments[hi].second - moments[lo].first)) *
                          (key - moments[lo].first));
            // fmt::print("{} -> {} <- {}\n", lo, pos, hi);
            // Condition of target found
            if (moments[pos].first <= key && key <= moments[pos].second)
                return pos;

            // If x is larger, x is in upper part
            if (moments[pos].second < key)
                lo = pos + 1;

                // If x is smaller, x is in the lower part
            else
                hi = pos - 1;
        }
        return -1;
    }


    int interpolationSearch(double key, const std::vector<Moment> &moments, int startIndex) {

        // Find indexes of two corners
        int lo = startIndex, hi = (moments.size() - 1);

        // Since array is sorted, an element present
        // in array must be in range defined by corner
        while (lo <= hi && key >= moments[lo].first && key <= moments[hi].second) {
            if (lo == hi) {
                if (moments[lo].first <= key && key <= moments[lo].second)
                    return lo;
                return -1;
            }
            // Probing the position with keeping
            // uniform distribution in mind.
            int pos =
                    lo + (((double) (hi - lo) / (moments[hi].second - moments[lo].first)) *
                          (key - moments[lo].first));
            // fmt::print("{} -> {} <- {}\n", lo, pos, hi);
            // Condition of target found
            if (moments[pos].first <= key && key <= moments[pos].second)
                return pos;

            // If x is larger, x is in upper part
            if (moments[pos].second < key)
                lo = pos + 1;

                // If x is smaller, x is in the lower part
            else
                hi = pos - 1;
        }
        return -1;
    }

    int jumpSearchProg(const std::vector<Moment> &arr, double noToSearch,
                       int ArrayLim) {
        int previous = 0;
        const int stepSize = std::sqrt(ArrayLim);
        int step = std::min(stepSize, ArrayLim - 1);
        // Step to skip elements for jumping

#ifdef DEBUG
        int controlIndex = LookupMomentIndex(noToSearch, arr);
#endif
        while (arr[step].first <= noToSearch && step < ArrayLim) {
            previous = step;
            step += stepSize;

            if (step > ArrayLim - 1) {
                step = ArrayLim;
/*#ifdef DEBUG
                if (-1 != controlIndex) {
                    fmt::print(fg(fmt::color::pale_violet_red), "[{:e}] pre [{:e} , {:e}] -- ", noToSearch,
                               arr[step - 1].first, arr[step - 1].second);
                }
#endif*/
                //return -1;
            }
        }

        /*Applying linear Search and starting from the previous elements*/
        while (arr[previous].second < noToSearch) {
            previous++;
            /*If element has not found yet then it means element is not present in the
             * array*/
            if (previous == std::min(step, ArrayLim)) {
#ifdef DEBUG
                if (-1 != controlIndex) {
                    fmt::print(fg(fmt::color::pale_violet_red), "[{:e}] pre [{:e} , {:e}] -- ", noToSearch,
                               arr[step - 1].first, arr[step - 1].second);
                }
#endif
                return -1;
            }
        }
        // if we found the element then
        if (arr[previous].first <= noToSearch && noToSearch <= arr[previous].second) {
            return previous;
        }

#ifdef DEBUG
        if (-1 != controlIndex) {
            fmt::print(fg(fmt::color::pale_violet_red), "[{:e}] pre [{:e} , {:e}] -- ", noToSearch,
                       arr[previous].first, arr[previous].second);
        }
#endif
        return -1;
    }

    int jumpSearchProg(const std::vector<Moment> &arr, double noToSearch,
                       int ArrayLim, int startIndex) {
        int previous = startIndex;
        const int stepSize = std::sqrt(ArrayLim - startIndex);
        int step = std::min(startIndex + stepSize, ArrayLim - 1);
        // Step to skip elements for jumping

        while (arr[step].first <= noToSearch && step < ArrayLim) {
            previous = step;
            step += stepSize;

            if (step > ArrayLim - 1) {
                step = ArrayLim;
            }
        }

        /*Applying linear Search and starting from the previous elements*/
        while (arr[previous].second < noToSearch) {
            previous++;
            /*If element has not found yet then it means element is not present in the
             * array*/
            if (previous == std::min(step, ArrayLim)) {
                return -1;
            }
        }
        // if we found the element then
        if (arr[previous].first <= noToSearch && noToSearch <= arr[previous].second) {
            return previous;
        }

        return -1;
    }

    int calcSearch(double key, const std::vector<Moment> &moments,
                   const std::vector<MomentInterval> &userMoments) {

        int start = -1; // TODO: ....
        // size_t indexOffset = 0;

#ifdef DEBUG
        int controlIndex = LookupMomentIndex(key, moments);
#endif
        for (auto &uMom: userMoments) {
            const double halfTimeWindow = uMom.timeWindow * 0.5;
            double calced_end = uMom.start + std::ceil((uMom.end - uMom.start) / uMom.interval) * uMom.interval + halfTimeWindow;
            if (key <= calced_end/*uMom.end + halfTimeWindow*/) {
                // fmt::print("Found {:e} <= {:e} (??? {:e} < ???)\n", key, uMom.end +
                // uMom.interval * 0.5, uMom.start - uMom.timeWindow * 0.5);

                if (key >= uMom.start - halfTimeWindow) {
                    // found it
                    const double nbMoments =
                            std::floor((uMom.end + 0.5 * uMom.timeWindow - uMom.start) / uMom.interval)/* + 1*/;
                            //std::ceil((uMom.end - uMom.start + uMom.timeWindow) / uMom.interval)/* + 1.0*/;
                    start = (int) std::min(((key - uMom.start + halfTimeWindow) / uMom.interval),
                                           nbMoments); // can go above limits on edge values
                    // fmt::print("[{:e}] Potential find at start {} ({}) [{:e} , {:e}]
                    // [[{:e} , {:e}]] [[[{:e} , {:e}]]]\n", key, start, nbMoments,
                    // moments[start].first, moments[start].second, uMom.start, uMom.end,
                    // uMom.start - uMom.timeWindow * 0.5 , uMom.end + uMom.timeWindow *
                    // 0.5);
                    if (key <= uMom.start + start * uMom.interval - halfTimeWindow ||
                        key >= uMom.start + start * uMom.interval + halfTimeWindow) {
                        // fmt::print("[{:e}] Calc not in window [{:e} , {:e}] , {} * {}
                        // ({})\n", key, uMom.start + start * uMom.interval - uMom.timeWindow,
                        // uMom.start + start * uMom.interval + uMom.timeWindow, start,
                        // nbMoments, uMom.interval); return -1;
                        start = -1;
                    } else {
                        if(uMom.startIndex) {
                            //start += 1;
                            start += uMom.startIndex;
                        }
                    }

                    if(start != -1 && start >= moments.size()) {
                        //fmt::print(fg(fmt::color::red), "{} / {} out of range for {}\n", start, moments.size(), key);
                        return -1;
                    }

#if 0
                    if (start != -1 &&
                        !(key >= moments[start].first && key <= moments[start].second)) {
                        // fmt::print("[{:e}] Calc passed [{:e} , {:e}] , {} * {} ({})\n",
                        // key, uMom.start + (start-uMom.startIndex) * uMom.interval -
                        // halfTimeWindow, uMom.start + (start-uMom.startIndex) * uMom.interval
                        // + halfTimeWindow, (start-uMom.startIndex), nbMoments,
                        // uMom.interval); fmt::print("[{:e}] Moments not in window [{:e} ,
                        // {:e}] [[{:e} , {:e}]]\n", key, moments[start].first,
                        // moments[start].second, uMom.start - halfTimeWindow, uMom.end +
                        // halfTimeWindow);
                        //start = -1;
                    }
#endif
#ifdef DEBUG
                    if (start != controlIndex) {
                        fmt::print(fg(fmt::color::pale_violet_red),
                                   "[{:e}] {} ({} - {}) vs {} [{:e} , {:e}] [[{:e} , {:e}]]\n", key,
                                start,
                                std::floor((key - uMom.start + halfTimeWindow) / uMom.interval),
                                (int) ((key - uMom.start + halfTimeWindow) / uMom.interval),
                                controlIndex, moments[start].first, moments[start].second,
                                moments[controlIndex].first, moments[controlIndex].second);
                        fmt::print(fg(fmt::color::pale_violet_red), "[{:e}] pre [{:e} , {:e}] -- ", key,
                                   moments[start - 1].first, moments[start - 1].second);
                        fmt::print(fg(fmt::color::pale_violet_red), "to [{:e} , {:e}]\n",
                                   moments[start + 1].first, moments[start + 1].second);
                    }
#endif
                    return start;
                }
                else {
#ifdef DEBUG
                    if (controlIndex != -1) {
                        fmt::print(fg(fmt::color::orange_red),
                                   "[{:e}] start but not end {} [{:e} , {:e}]\n", key,
                                   controlIndex, moments[controlIndex].first,
                                   moments[controlIndex].second);
                    }
#endif

                    //int ctrl = LookupMomentIndex(key, moments);
                    /*if(controlIndex != -1) {
                        fmt::print(fg(fmt::color::orange_red),
                                   "[{}] start but not end {} - {}\n", key,
                                   uMom.start, halfTimeWindow);
                    }*/
                    return -2;
                }
            }
            // indexOffset += ((uMom.end - uMom.start) / uMom.interval + 1);
            // fmt::print("Indexoffset {}\n", indexOffset);
        }
        /*if(controlIndex != -1) {
            fmt::print("[{:e}] no start no end {} [{:e} , {:e}] [[{:e} , {:e}]]\n",
        key, controlIndex, moments[controlIndex].first, moments[controlIndex].second);
            fmt::print("[{:e}] first [{:e} , {:e}]\n", key, moments.front().first,
        moments.front().second); fmt::print("[{:e}] last [{:e} , {:e}]\n", key,
        moments.back().first, moments.back().second);
        }*/
        return -3;
    }

} // namespace MFTD