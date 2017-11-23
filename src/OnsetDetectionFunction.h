//=======================================================================
/** @file OnsetDetectionFunction.h
 *  @brief A class for calculating onset detection functions
 *  @author Zach Banks, Adam Stark
 *  @copyright Copyright (C) 2015 Zach Banks
 *  @copyright Copyright (C) 2008-2014  Queen Mary University of London
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
//=======================================================================

#ifndef __ONSETDETECTIONFUNCTION_H
#define __ONSETDETECTIONFUNCTION_H

#include "fftw3.h"
#include "common.h"

#include <memory>
#include <utility>
#include <functional>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iterator>
#include "Plan.hpp"

//=======================================================================
/** The type of onset detection function to calculate */
enum OnsetDetectionFunctionType
{
    EnergyEnvelope,
    EnergyDifference,
    SpectralDifference,
    SpectralDifferenceHWR,
    PhaseDeviation,
    ComplexSpectralDifference,
    ComplexSpectralDifferenceSq,
    ComplexSpectralDifferenceHWR,
    ComplexSpectralDifferenceSqHWR,
    HighFrequencyContent,
    HighFrequencySpectralDifference,
    HighFrequencySpectralDifferenceHWR
};

//=======================================================================
/** The type of window to use when calculating onset detection function samples */
enum WindowType
{
    RectangularWindow,
    HanningWindow,
    HammingWindow,
    BlackmanWindow,
    TukeyWindow
};

//=======================================================================
/** A class for calculating onset detection functions. */
struct odf {
    using Plan = retrack::Plan;
    using float_type = Plan::float_type;
	int frameSize;						/**< audio framesize */
	int hopSize;						/**< audio hopsize */
	OnsetDetectionFunctionType type;		/**< type of detection function */
    WindowType windowType;                     /**< type of window used in calculations */
    Plan p;
    retrack::detail::fftwf_ptr<float_type> complexIn;
    retrack::detail::fftwf_ptr<float_type> complexOut;

	int initialised;					/**< flag indicating whether buffers and FFT plans are initialised */

    std::unique_ptr<float[]> frame;
    std::unique_ptr<float[]> window;                     /**< audio frame */

	float prevEnergySum;				/**< to hold the previous energy sum value */

    std::unique_ptr<float[]> magSpec;                   /**< magnitude spectrum */
    std::unique_ptr<float[]> prevMagSpec;               /**< previous magnitude spectrum */

    std::unique_ptr<float[]> phase;                     /**< FFT phase values */
    std::unique_ptr<float[]> prevPhase;                 /**< previous phase values */
    std::unique_ptr<float[]> prevPhase2;                 /**< second order previous phase values */

    odf() = default;;
    odf(int hop_size, int frame_size, OnsetDetectionFunctionType odf_type, WindowType window_type);
  ~odf();
    float process_frame(const btrack_chunk_t * buffer);
    float process_fft_frame(const btrack_chunk_t * fft_buffer); // XXX: This is not implemented yet
    void set_type(enum OnsetDetectionFunctionType type);

};

/** Constructor
* @param hopSize_ the hop size in audio samples
* @param frameSize_ the frame size in audio samples
* @param onsetDetectionFunctionType_ the type of onset detection function to use - (see OnsetDetectionFunctionType)
* @param windowType the type of window to use (see WindowType)
*/


#endif
