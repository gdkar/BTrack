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
#include "sse_mathfun.h"
#include <string.h>
#include "common.h"
#include <fftw3.h>
#include <math.h>

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
    ComplexSpectralDifferenceHWR,
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
	int frameSize;						/**< audio framesize */
	int hopSize;						/**< audio hopsize */
	enum OnsetDetectionFunctionType type;		/**< type of detection function */
    enum WindowType windowType;                     /**< type of window used in calculations */
	
	fftwf_plan p;						/**< fftw plan */
	float  *realIn;			/**< to hold complex fft values for input */
	float  *imagIn;			/**< to hold complex fft values for input */
  float  *realOut;			/**< to hold complex fft values for input */
	float  *imagOut;			/**< to hold complex fft values for input */

	int initialised;					/**< flag indicating whether buffers and FFT plans are initialised */

    float * frame;                     /**< audio frame */
    float * window;                    /**< window */
	
	float prevEnergySum;				/**< to hold the previous energy sum value */
	
    float * magSpec;                   /**< magnitude spectrum */
    float * prevMagSpec;               /**< previous magnitude spectrum */
	
    float * phase;                     /**< FFT phase values */
    float * prevPhase;                 /**< previous phase values */
    float * prevPhase2;                 /**< second order previous phase values */

};

/** Constructor 
* @param hopSize_ the hop size in audio samples
* @param frameSize_ the frame size in audio samples
* @param onsetDetectionFunctionType_ the type of onset detection function to use - (see OnsetDetectionFunctionType)
* @param windowType the type of window to use (see WindowType)
*/
int odf_init(struct odf * odf, int hop_size, int frame_size, enum OnsetDetectionFunctionType odf_type, enum WindowType window_type);
void odf_del(struct odf * odf);

float odf_calculate_sample(struct odf * odf, float * buffer);
void odf_set_type(struct odf * odf, enum OnsetDetectionFunctionType type);

#endif
