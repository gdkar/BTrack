//=======================================================================
/** @file OnsetDetectionFunction.h
 *  @brief A class for calculating onset detection functions
 *  @author Adam Stark
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

#include <cmath>
#include <vector>
#include <fftw3.h>
#include <ffts/ffts.h>
#include "sse_mathfun.h"

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
class OnsetDetectionFunction
{
public:
    /** Constructor 
     * @param hopSize_ the hop size in audio samples
     * @param frameSize_ the frame size in audio samples
     * @param onsetDetectionFunctionType_ the type of onset detection function to use - (see OnsetDetectionFunctionType)
     * @param windowType the type of window to use (see WindowType)
     */
	OnsetDetectionFunction(int hopSize_
      ,int frameSize_
      ,int onsetDetectionFunctionType_ = ComplexSpectralDifferenceHWR
      ,int windowType_ = HanningWindow);
    
    /** Destructor */
	~OnsetDetectionFunction();
	
    /** Process input frame and calculate detection function sample 
     * @param buffer a pointer to an array containing the audio samples to be processed
     * @returns the onset detection function sample
     */
	float calculateOnsetDetectionFunctionSample(float *buffer);
    
    /** Set the detection function type 
     * @param onsetDetectionFunctionType_ the type of onset detection function to use - (see OnsetDetectionFunctionType)
     */
	void setOnsetDetectionFunctionType(int onsetDetectionFunctionType_);
	
private:
	
    /** Perform the FFT on the data in 'frame' */
	void performFFT();
  void toPolar();

    //=======================================================================
    /** Calculate energy envelope detection function sample */
	float energyEnvelope();
    
    /** Calculate energy difference detection function sample */
	float energyDifference();
    
    /** Calculate spectral difference detection function sample */
	float spectralDifference();
    
    /** Calculate spectral difference (half wave rectified) detection function sample */
	float spectralDifferenceHWR();
    
    /** Calculate phase deviation detection function sample */
	float phaseDeviation();
    
    /** Calculate complex spectral difference detection function sample */
	float complexSpectralDifference();
    
    /** Calculate complex spectral difference detection function sample (half-wave rectified) */
	float complexSpectralDifferenceHWR();
    
    /** Calculate high frequency content detection function sample */
	float highFrequencyContent();
    
    /** Calculate high frequency spectral difference detection function sample */
	float highFrequencySpectralDifference();
    
    /** Calculate high frequency spectral difference detection function sample (half-wave rectified) */
	float highFrequencySpectralDifferenceHWR();

    //=======================================================================
    /** Calculate a Rectangular window */
	void calculateRectangularWindow();
    
    /** Calculate a Hanning window */
	void calculateHanningWindow();
    
    /** Calculate a Hamming window */
	void calclulateHammingWindow();
    
    /** Calculate a Blackman window */
	void calculateBlackmanWindow();
    
    /** Calculate a Tukey window */
	void calculateTukeyWindow();

    //=======================================================================
	/** Set phase values between [-pi, pi] 
     * @param phaseVal the phase value to process
     * @returns the wrapped phase value
     */
	float princarg(float phaseVal);
  v4sf _mm_princarg_ps(v4sf);	
	
	int frameSize;						/**< audio framesize */
	int hopSize;						/**< audio hopsize */
	int onsetDetectionFunctionType;		/**< type of detection function */
  int windowType;                     /**< type of window used in calculations */
	
	fftwf_plan p;						/**< fftw plan */
  float         *realIn;
  float         *imagIn;
  float         *realOut;
  float         *imagOut;

  std::vector<float> frame;          /**< audio frame */
  std::vector<float> window;         /**< window */
	
	v4sf prevEnergySum;				/**< to hold the previous energy sum value */
	
  std::vector<float> magSpec;        /**< magnitude spectrum */
  std::vector<float> prevMagSpec;    /**< previous magnitude spectrum */
	
  std::vector<float> phase;          /**< FFT phase values */
  std::vector<float> prevPhase;      /**< previous phase values */
  std::vector<float> prevPhase2;     /**< second order previous phase values */

};


#endif
