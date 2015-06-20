//=======================================================================
/** @file OnsetDetectionFunction.c
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
#define _ISOC11_SOURCE
#define _GNU_SOURCE
#define _XOPEN_SOURCE
#include <math.h>
#include "OnsetDetectionFunction.h"
#include "common.h"

static void performFFT(struct odf * odf);
static double princarg(double phaseVal);

static double energyEnvelope(struct odf * odf);
static double energyDifference(struct odf * odf);
static double spectralDifference(struct odf * odf);
static double spectralDifferenceHWR(struct odf * odf);
static double phaseDeviation(struct odf * odf);
static double complexSpectralDifference(struct odf * odf);
static double complexSpectralDifferenceHWR(struct odf * odf);
static double highFrequencyContent(struct odf * odf);
static double highFrequencySpectralDifference(struct odf * odf);
static double highFrequencySpectralDifferenceHWR(struct odf * odf);

static void calculateRectangularWindow(double * window, int frameSize);	
static void calculateHanningWindow(double * window, int frameSize);		
static void calclulateHammingWindow(double * window, int frameSize);		
static void calculateBlackmanWindow(double * window, int frameSize);		
static void calculateTukeyWindow(double * window, int frameSize);          

int odf_init(struct odf * odf, int hop_size, int frame_size, enum OnsetDetectionFunctionType odf_type, enum WindowType window_type){
    // if we have already initialised FFT plan
    odf_del(odf);
	
	odf->hopSize = hop_size; // set hopsize
	odf->frameSize = frame_size; // set framesize
	
    odf->type = odf_type;
    odf->windowType = window_type;
		
	// initialise buffers
    odf->frame = malloc(sizeof(double) * frame_size);
    odf->window = malloc(sizeof(double) * frame_size);
    odf->magSpec = malloc(sizeof(double) * frame_size);
    odf->prevMagSpec = malloc(sizeof(double) * frame_size);
    odf->phase = malloc(sizeof(double) * frame_size);
    odf->prevPhase = malloc(sizeof(double) * frame_size);
    odf->prevPhase2 = malloc(sizeof(double) * frame_size);
    ASSERT(odf->frame && odf->window && odf->magSpec && odf->prevMagSpec && odf->phase && odf->prevPhase && odf->prevPhase2);
	
	// set the window to the specified type
	switch (window_type){
		case RectangularWindow:
			calculateRectangularWindow(odf->window, frame_size);		// Rectangular window
			break;	
		case HanningWindow:
			calculateHanningWindow(odf->window, frame_size);			// Hanning Window
			break;
		case HammingWindow:
			calclulateHammingWindow(odf->window, frame_size);			// Hamming Window
			break;
		case BlackmanWindow:
			calculateBlackmanWindow(odf->window, frame_size);			// Blackman Window
			break;
		case TukeyWindow:
			calculateTukeyWindow(odf->window, frame_size);             // Tukey Window
			break;
		default:
			calculateHanningWindow(odf->window, frame_size);			// DEFAULT: Hanning Window
	}
	
	// initialise previous magnitude spectrum to zero
	for (int i = 0;i < frame_size;i++) {
		odf->prevMagSpec[i] = 0.0;
		odf->prevPhase[i] = 0.0;
		odf->prevPhase2[i] = 0.0;
		odf->frame[i] = 0.0;
	}
	
	odf->prevEnergySum = 0.0;	// initialise previous energy sum value to zero
	
	/*  Init fft */
	odf->complexIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * frame_size);		// complex array to hold fft data
	odf->complexOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * frame_size);	// complex array to hold fft data
	odf->p = fftw_plan_dft_1d(frame_size, odf->complexIn, odf->complexOut, FFTW_FORWARD, FFTW_ESTIMATE);	// FFT plan initialisation
	
	odf->initialised = true;
    return 0;
}

void odf_del(struct odf * odf){
    if (odf->initialised) {
        // destroy fft plan
        fftw_destroy_plan(odf->p);
        fftw_free(odf->complexIn);
        fftw_free(odf->complexOut);
    }
}

void odf_set_type(struct odf * odf, enum OnsetDetectionFunctionType type){
    odf->type = type;
}

double odf_calculate_sample(struct odf * odf, double * buffer){
	double odfSample;
	
  memmove(odf->frame,odf->frame+odf->hopSize,(odf->frameSize-odf->hopSize)*sizeof(odf->frame[0]));
	
	// add new samples to frame from input buffer
  memmove(odf->frame + (odf->frameSize-odf->hopSize),buffer,odf->hopSize*sizeof(buffer[0]));
		
	switch (odf->type){
		case EnergyEnvelope:
            // calculate energy envelope detection function sample
			odfSample = energyEnvelope(odf);
			break;
		case EnergyDifference:
            // calculate half-wave rectified energy difference detection function sample
			odfSample = energyDifference(odf);
			break;
		case SpectralDifference:
            // calculate spectral difference detection function sample
			odfSample = spectralDifference(odf);
			break;
		case SpectralDifferenceHWR:
            // calculate spectral difference detection function sample (half wave rectified)
			odfSample = spectralDifferenceHWR(odf);
			break;
		case PhaseDeviation:
            // calculate phase deviation detection function sample (half wave rectified)
			odfSample = phaseDeviation(odf);
			break;
		case ComplexSpectralDifference:
            // calcualte complex spectral difference detection function sample
			odfSample = complexSpectralDifference(odf);
			break;
		case ComplexSpectralDifferenceHWR:
            // calcualte complex spectral difference detection function sample (half-wave rectified)
			odfSample = complexSpectralDifferenceHWR(odf);
			break;
		case HighFrequencyContent:
            // calculate high frequency content detection function sample
			odfSample = highFrequencyContent(odf);
			break;
		case HighFrequencySpectralDifference:
            // calculate high frequency spectral difference detection function sample
			odfSample = highFrequencySpectralDifference(odf);
			break;
		case HighFrequencySpectralDifferenceHWR:
            // calculate high frequency spectral difference detection function (half-wave rectified)
			odfSample = highFrequencySpectralDifferenceHWR(odf);
			break;
		default:
			odfSample = 1.0;
	}
		
	return odfSample;
}

static void performFFT(struct odf * odf) {
	int fsize2 = (odf->frameSize/2);
	
	// window frame and copy to complex array, swapping the first and second half of the signal
	for (int i = 0;i < fsize2;i++)
	{
		odf->complexIn[i][0] = odf->frame[i+fsize2] * odf->window[i+fsize2];
		odf->complexIn[i][1] = 0.0;
		odf->complexIn[i+fsize2][0] = odf->frame[i] * odf->window[i];
		odf->complexIn[i+fsize2][1] = 0.0;
	}
	
	// perform the fft
	fftw_execute(odf->p);
}

////////////////////////////// Methods for Detection Functions /////////////////////////////////

static double energyEnvelope(struct odf * odf){
	double sum = 0;	
	
	// sum the squares of the samples
	for (int i = 0;i < odf->frameSize;i++) {
		sum = sum + (odf->frame[i]*odf->frame[i]);
	}
	
	return sum;		// return sum
}
static void genMagnitudeSpectrum(struct odf * odf ){
	performFFT(odf);
	for (int i = 0;i < (odf->frameSize/2)+1;i++) {
		double thisSpec = hypot(odf->complexOut[i][0],odf->complexOut[i][1]);
    odf->magSpec[i] = thisSpec;
    odf->magSpec[odf->frameSize-1-i] = thisSpec;
	}
}
static double energyDifference(struct odf * odf){
	double sum = 0;
	double sample;
	
	// sum the squares of the samples
	for (int i = 0;i < odf->frameSize;i++) {
		sum = sum + (odf->frame[i]*odf->frame[i]);
	}
	sample = sum - odf->prevEnergySum;	// sample is first order difference in energy
  odf->prevEnergySum = sum;

  return MIN(sample, 0);
}

static double spectralDifference(struct odf * odf){
	double sum = 0;
	
	// perform the FFT
  genMagnitudeSpectrum(odf);
	
	// compute first (N/2)+1 mag values
	for (int i = 0;i < odf->frameSize;i++) {
    sum += 2*(odf->prevMagSpec[i]-odf->magSpec[i]);
    odf->prevMagSpec[i] = odf->magSpec[i];;
	}
	return sum;		
}

static double spectralDifferenceHWR(struct odf * odf) {
	double sum = 0;

	// perform the FFT
  genMagnitudeSpectrum(odf);
	
	sum = 0;	// initialise sum to zero
	for (int i = 0;i < odf->frameSize;i++) {
    double diff = odf->magSpec[i]-odf->prevMagSpec[i];
    sum += MAX(diff,0);
    odf->prevMagSpec[i]=odf->magSpec[i];
  }
	return sum;		
}

static double phaseDeviation(struct odf * odf){
	double dev;
	double sum = 0;
	
	// perform the FFT
  genMagnitudeSpectrum(odf);	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// calculate phase value
		odf->phase[i] = atan2(odf->complexOut[i][1],odf->complexOut[i][0]);
		
		// if bin is not just a low energy bin then examine phase deviation
		if (odf->magSpec[i] > 0.1) {
			dev = odf->phase[i] - (2*odf->prevPhase[i]) + odf->prevPhase2[i];	// phase deviation
			sum += fabs(princarg(dev));	// wrap into [-pi,pi] range
		}
				
		// store values for next calculation
		odf->prevPhase2[i] = odf->prevPhase[i];
		odf->prevPhase[i] = odf->phase[i];
	}
	
	return sum;		
}

double complexSpectralDifference(struct odf * odf){
	double dev,pdev;
	double sum = 0;
	double mag_diff,phase_diff;
	double value;
	
	// perform the FFT
	genMagnitudeSpectrum(odf);
	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// calculate phase value
		odf->phase[i] = atan2(odf->complexOut[i][1],odf->complexOut[i][0]);
		
				// phase deviation
		dev = odf->phase[i] - (2*odf->prevPhase[i]) + odf->prevPhase2[i];
		
		// wrap into [-pi,pi] range
		pdev = princarg(dev);	
		
		// calculate magnitude difference (real part of Euclidean distance between complex frames)
		mag_diff = odf->magSpec[i] - odf->prevMagSpec[i];
		
		// calculate phase difference (imaginary part of Euclidean distance between complex frames)
		phase_diff = -odf->magSpec[i]*sin(pdev);
		
		// square real and imaginary parts, sum and take square root
		sum += value = hypot(mag_diff,phase_diff);
	
		// store values for next calculation
		odf->prevPhase2[i] = odf->prevPhase[i];
		odf->prevPhase[i] = odf->phase[i];
		odf->prevMagSpec[i] = odf->magSpec[i];
	}
	
	return sum;		
}

double complexSpectralDifferenceHWR(struct odf * odf) {
	double dev,pdev;
	double sum = 0;
	double mag_diff,phase_diff;
	
	genMagnitudeSpectrum(odf);
	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// calculate phase value
		odf->phase[i] = atan2(odf->complexOut[i][1],odf->complexOut[i][0]);
		
		// phase deviation
		dev = odf->phase[i] - (2*odf->prevPhase[i]) + odf->prevPhase2[i];
		
		// wrap into [-pi,pi] range
		pdev = princarg(dev);	
		
		// calculate magnitude difference (real part of Euclidean distance between complex frames)
		mag_diff = odf->magSpec[i] - odf->prevMagSpec[i];
		
		// if we have a positive change in magnitude, then include in sum, otherwise ignore (half-wave rectification)
		if (mag_diff > 0) {
			// calculate phase difference (imaginary part of Euclidean distance between complex frames)
			phase_diff = -odf->magSpec[i]*sin(pdev);

			// square real and imaginary parts, sum and take square root
			sum += hypot(mag_diff,phase_diff);
		
		}
		
		// store values for next calculation
		odf->prevPhase2[i] = odf->prevPhase[i];
		odf->prevPhase[i] = odf->phase[i];
		odf->prevMagSpec[i] = odf->magSpec[i];
	}
	
	return sum;		
}

double highFrequencyContent(struct odf * odf) {
	double sum = 0;
	
	// perform the FFT
	genMagnitudeSpectrum(odf);
	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {		
		
		sum = sum + (odf->magSpec[i]*((double) (i+1)));
		
		// store values for next calculation
		odf->prevMagSpec[i] = odf->magSpec[i];
	}
	
	return sum;		
}

static double highFrequencySpectralDifference(struct odf * odf) {
	double sum = 0;
	
	// perform the FFT
	genMagnitudeSpectrum(odf);
	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {		
		// calculate difference
		sum += (i+1)*fabs(odf->magSpec[i] - odf->prevMagSpec[i]);
		
		// store values for next calculation
		odf->prevMagSpec[i] = odf->magSpec[i];
	}
	
	return sum;		
}

static double highFrequencySpectralDifferenceHWR(struct odf * odf) {
	double sum = 0;
	
	// perform the FFT
	genMagnitudeSpectrum(odf);
	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {		
		// calculate difference
		sum += (i+1)*MAX(odf->magSpec[i] - odf->prevMagSpec[i],0);
		
		// store values for next calculation
		odf->prevMagSpec[i] = odf->magSpec[i];
	}
	
	return sum;		
}

////////////////////////////// Methods to Calculate Windows ////////////////////////////////////

static void calculateHanningWindow(double * window, int frameSize){
	double N = (double) (frameSize-1);	// framesize minus 1
	
	// Hanning window calculation
	for (int n = 0;n < frameSize;n++) {
		window[n] = 0.5*(1-cos(2*M_PI*(n/N)));
	}
}

static void calclulateHammingWindow(double * window, int frameSize) {
	double N = (double) (frameSize-1);	// framesize minus 1
	double n_val = 0;// double version of index 'n'
	
	// Hamming window calculation
	for (int n = 0;n < frameSize;n++) {
		window[n] = 0.54 - (0.46*cos(2*M_PI*(n_val/N)));
		n_val = n_val+1;
	}
}

static void calculateBlackmanWindow(double * window, int frameSize) {
	double N = (double) (frameSize-1);	// framesize minus 1
	double n_val = 0;	// double version of index 'n'
	
	// Blackman window calculation
	for (int n = 0;n < frameSize;n++) {
		window[n] = 0.42 - (0.5*cos(2*M_PI*(n_val/N))) + (0.08*cos(4*M_PI*(n_val/N)));
		n_val = n_val+1;
	}
}

static void calculateTukeyWindow(double * window, int frameSize) {
	double alpha = 0.5;	// alpha [default value = 0.5];
	double N = (double) (frameSize-1);	// framesize minus 1
	double n_val = (double) (-1*((frameSize/2)))+1;
		
	// Tukey window calculation
    // left taper
	for (int n = 0;n < frameSize;n++) {
		if ((n_val >= 0) && (n_val <= (alpha*(N/2)))) {
			window[n] = 1.0;
		} else if ((n_val <= 0) && (n_val >= (-1*alpha*(N/2)))) {
			window[n] = 1.0;
		} else {
			window[n] = 0.5*(1+cos(M_PI*(((2*n_val)/(alpha*N))-1)));
		}

		n_val = n_val+1;			 
	}
}

static void calculateRectangularWindow(double * window, int frameSize) {
	for (int n = 0;n < frameSize;n++) {
		window[n] = 1.0;
	}
}


///////////////////////////////// Other Handy Methods //////////////////////////////////////////

static inline double princarg(double phaseVal) {	
	// if phase value is less than or equal to -pi then add 2*pi
  return phaseVal - (2*M_PI)*floor(phaseVal/(2*M_PI));
}
