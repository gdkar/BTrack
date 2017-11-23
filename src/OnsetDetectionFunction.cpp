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

#include "OnsetDetectionFunction.h"
#include "common.h"
#include "Math.hpp"
#include "SysUtils.hpp"

using namespace retrack;

static void performFFT(struct odf * odf);

static float energyEnvelope(struct odf * odf);
static float energyDifference(struct odf * odf);
static float spectralDifference(struct odf * odf);
static float spectralDifferenceHWR(struct odf * odf);
static float phaseDeviation(struct odf * odf);
static float complexSpectralDifference(struct odf * odf);
static float complexSpectralDifferenceSq(struct odf * odf);
static float complexSpectralDifferenceHWR(struct odf * odf);
static float complexSpectralDifferenceSqHWR(struct odf * odf);
static float highFrequencyContent(struct odf * odf);
static float highFrequencySpectralDifference(struct odf * odf);
static float highFrequencySpectralDifferenceHWR(struct odf * odf);

static void calculateRectangularWindow(float * window, int frameSize);
static void calculateHanningWindow(float * window, int frameSize);
static void calclulateHammingWindow(float * window, int frameSize);
static void calculateBlackmanWindow(float * window, int frameSize);
static void calculateTukeyWindow(float * window, int frameSize);

odf::odf(int hop_size, int frame_size, float rate, OnsetDetectionFunctionType odf_type, WindowType window_type)
: sampleRate{rate}{
    // if we have already initialised FFT plan
    //odf_del(odf);

    BTRACK_ASSERT(frame_size >= hop_size);
	hopSize = hop_size; // set hopsize
	frameSize = frame_size; // set framesize

    type = odf_type;
    windowType = window_type;

	// initialise buffers
    frame = SlideBuffer<float>(frame_size);
//    frame = std::make_unique<float[]>(frame_size);
    window = std::make_unique<float[]>(frame_size);
    magSpec= std::make_unique<float[]>(frame_size);
    prevMagSpec= std::make_unique<float[]>(frame_size);
    phase= std::make_unique<float[]>(frame_size);
    prevPhase= std::make_unique<float[]>(frame_size);
    prevPhase2= std::make_unique<float[]>(frame_size);
	// set the window to the specified type
	switch (window_type){
		case RectangularWindow:
			calculateRectangularWindow(window.get(), frame_size);		// Rectangular window
			break;
		case HanningWindow:
			calculateHanningWindow(window.get(), frame_size);			// Hanning Window
			break;
		case HammingWindow:
			calclulateHammingWindow(window.get(), frame_size);			// Hamming Window
			break;
		case BlackmanWindow:
			calculateBlackmanWindow(window.get(), frame_size);			// Blackman Window
			break;
		case TukeyWindow:
			calculateTukeyWindow(window.get(), frame_size);             // Tukey Window
			break;
		default:
			calculateHanningWindow(window.get(), frame_size);			// DEFAULT: Hanning Window
	}

	prevEnergySum = 0.0;	// initialise previous energy sum value to zero

	complexIn = detail::make_fftwf<float_type>(frame_size );
	complexOut = detail::make_fftwf<float_type>(frame_size * 2);
    /*  Init fft */
//	complexIn = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * frame_size);		// complex array to hold fft data
//	complexOut = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * frame_size);	// complex array to hold fft data
    p = Plan::dft_1d_r2c(frame_size,complexIn.get(),complexOut.get(),complexOut.get() + frame_size,  FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
//	p = fftwf_plan_dft_1d(frame_size, complexIn, complexOut, FFTW_FORWARD, FFTW_ESTIMATE);	// FFT plan initialisation

}

odf::~odf() = default;

void odf::set_type(enum OnsetDetectionFunctionType _type){
    type = _type;
}

float odf::process_frame(const btrack_chunk_t * buffer){
	float odfSample;

	// shift audio samples back in frame by hop size
    frame.append(buffer,buffer + hopSize);
/*    auto req = frame.request(
	for (int i = 0; i < (frameSize-hopSize);i++) {
		frame[i] = frame[i+hopSize];
	}

	// add new samples to frame from input buffer
	int j = 0;
	for (int i = (frameSize-hopSize);i < frameSize;i++) {
		frame[i] = buffer[j];
		j++;
	}*/

	switch (type){
		case EnergyEnvelope:
            // calculate energy envelope detection function sample
			odfSample = energyEnvelope(this);
			break;
		case EnergyDifference:
            // calculate half-wave rectified energy difference detection function sample
			odfSample = energyDifference(this);
			break;
		case SpectralDifference:
            // calculate spectral difference detection function sample
			odfSample = spectralDifference(this);
			break;
		case SpectralDifferenceHWR:
            // calculate spectral difference detection function sample (half wave rectified)
			odfSample = spectralDifferenceHWR(this);
			break;
		case PhaseDeviation:
            // calculate phase deviation detection function sample (half wave rectified)
			odfSample = phaseDeviation(this);
			break;
		case ComplexSpectralDifference:
            // calcualte complex spectral difference detection function sample
			odfSample = complexSpectralDifference(this);
			break;
		case ComplexSpectralDifferenceSq:
            // calcualte complex spectral difference squared detection function sample
			odfSample = complexSpectralDifferenceSq(this);
			break;
		case ComplexSpectralDifferenceHWR:
            // calcualte complex spectral difference detection function sample (half-wave rectified)
			odfSample = complexSpectralDifferenceHWR(this);
			break;
		case ComplexSpectralDifferenceSqHWR:
            // calcualte complex spectral difference detection squared function sample (half-wave rectified)
			odfSample = complexSpectralDifferenceSqHWR(this);
			break;
		case HighFrequencyContent:
            // calculate high frequency content detection function sample
			odfSample = highFrequencyContent(this);
			break;
		case HighFrequencySpectralDifference:
            // calculate high frequency spectral difference detection function sample
			odfSample = highFrequencySpectralDifference(this);
			break;
		case HighFrequencySpectralDifferenceHWR:
            // calculate high frequency spectral difference detection function (half-wave rectified)
			odfSample = highFrequencySpectralDifferenceHWR(this);
			break;
		default:
			odfSample = 1.0;
	}

	return odfSample;
}

float odf::process_fft_frame(const btrack_chunk_t * buffer) {
    // XXX TODO
    return 0;
}

static void performFFT(struct odf * odf) {
    auto fsize = odf->frameSize;
    auto fsize2 = fsize / 2;
    auto fsize2p = fsize2 + 1;

    auto r_out = odf->complexOut.get();
    auto i_out = r_out + fsize;

	// window frame and copy to complex array, swapping the first and second half of the signal
    auto rd = odf->frame.kept_range();
    std::transform(
        std::begin(rd)
       ,std::end(rd)
       ,&odf->window[0]
       ,&odf->complexIn[0]
       ,[](auto x,auto y){return x*y;}
        );

    std::rotate(&odf->complexIn[0],&odf->complexIn[0] + fsize2, &odf->complexIn[0] + fsize);
	// perform the fft
    odf->p.execute(odf->complexIn.get(),r_out);
	// compute first (N/2)+1 mag values
	for (int i = 0;i < fsize2p;i++) {
		// calculate phase value
        auto _i = i_out[i];
        auto _r = r_out[i];
        auto _a = std::atan2(_i,_r);
        auto _m = std::hypot(_i,_r);
        odf->phase[i] = _a;
        odf->magSpec[i] = _m;
        if(i) {
            odf->phase[fsize-i] = -_a;
            odf->magSpec[fsize-i] = _m;
        }
	}
}

////////////////////////////// Methods for Detection Functions /////////////////////////////////

static float energyEnvelope(struct odf * odf){
    auto rd = odf->frame.kept_range();
    return std::inner_product(rd.begin(),rd.end(),rd.begin(),0.f);
}

static float energyDifference(struct odf * odf){
    auto rd = odf->frame.kept_range();
    auto sum = std::inner_product(rd.begin(),rd.end(),rd.begin(),0.f);
    using std::swap;
    swap(sum,odf->prevEnergySum);
    return std::min(0.f,odf->prevEnergySum - sum);	// sample is first order difference in energy
}

static float spectralDifference(struct odf * odf){
	float diff;
	float sum;

	// perform the FFT
	performFFT(odf);
	sum = 0;	// initialise sum to zero
	for (int i = 0;i < odf->frameSize;i++) {
		// calculate difference
		diff = odf->magSpec[i] - odf->prevMagSpec[i];

		// ensure all difference values are positive
		if (diff < 0) {
			diff = diff*-1;
		}

		// add difference to sum
		sum = sum+diff;

		// store magnitude spectrum bin for next detection function sample calculation
		odf->prevMagSpec[i] = odf->magSpec[i];
	}

	return sum;
}

static float spectralDifferenceHWR(struct odf * odf) {
	float diff;
	float sum;

	// perform the FFT
	performFFT(odf);
	sum = 0;	// initialise sum to zero
	for (int i = 0;i < odf->frameSize;i++) {
		// calculate difference
		diff = odf->magSpec[i] - odf->prevMagSpec[i];

		// only add up positive differences
		if (diff > 0)
		{
			// add difference to sum
			sum = sum+diff;
		}

		// store magnitude spectrum bin for next detection function sample calculation
		odf->prevMagSpec[i] = odf->magSpec[i];
	}

	return sum;
}

static float phaseDeviation(struct odf * odf){
	float dev;
    float pdev;
	float sum = 0;
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// calculate phase value

		// if bin is not just a low energy bin then examine phase deviation
		if (odf->magSpec[i] > 0.1) {
			dev = odf->phase[i] - (2*odf->prevPhase[i]) + odf->prevPhase2[i];	// phase deviation
			pdev = princarg(dev);	// wrap into [-pi,pi] range

			// make all values positive
			if (pdev < 0){
				pdev = pdev*-1;
			}

			// add to sum
			sum = sum + pdev;
		}

		// store values for next calculation
		odf->prevPhase2[i] = odf->prevPhase[i];
		odf->prevPhase[i] = odf->phase[i];
	}

	return sum;
}

float complexSpectralDifference(struct odf * odf){
	float dev,pdev;
	float sum = 0;
	float value;

	// perform the FFT
	performFFT(odf);

	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// phase deviation
		dev = odf->phase[i] - (2*odf->prevPhase[i]) + odf->prevPhase2[i];

		// wrap into [-pi,pi] range
		pdev = princarg(dev);

        // Calculate the euclidean distance between previous frame & expected current frame
        value = std::sqrt(sqr(odf->magSpec[i]) + sqr(odf->prevMagSpec[i]) - 2 * odf->magSpec[i] * odf->prevMagSpec[i] * std::cos(pdev));

		// add to sum
		sum = sum + value;

		// store values for next calculation
		odf->prevPhase2[i] = odf->prevPhase[i];
		odf->prevPhase[i] = odf->phase[i];
		odf->prevMagSpec[i] = odf->magSpec[i];
	}

	return sum;
}

float complexSpectralDifferenceSq(struct odf * odf){
	float dev,pdev;
	float sum = 0;
	float value;

	// perform the FFT
	performFFT(odf);

	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// phase deviation
		dev = odf->phase[i] - (2*odf->prevPhase[i]) + odf->prevPhase2[i];

		// wrap into [-pi,pi] range
		pdev = princarg(dev);

        // Calculate the euclidean distance squared between previous frame & expected current frame
        value = (sqr(odf->magSpec[i]) + sqr(odf->prevMagSpec[i]) - 2 * odf->magSpec[i] * odf->prevMagSpec[i] * std::cos(pdev));

		// add to sum
		sum = sum + value;

		// store values for next calculation
		odf->prevPhase2[i] = odf->prevPhase[i];
		odf->prevPhase[i] = odf->phase[i];
		odf->prevMagSpec[i] = odf->magSpec[i];
	}

	return sum;
}

float complexSpectralDifferenceHWR(struct odf * odf) {
	float dev,pdev;
	float sum = 0;
	float mag_diff;
	float value;

	// perform the FFT
	performFFT(odf);

	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		dev = odf->phase[i] - (2*odf->prevPhase[i]) + odf->prevPhase2[i];

		// wrap into [-pi,pi] range
		pdev = princarg(dev);

		// calculate magnitude difference (real part of Euclidean distance between complex frames)
		mag_diff = odf->magSpec[i] - odf->prevMagSpec[i];

		// if we have a positive change in magnitude, then include in sum, otherwise ignore (half-wave rectification)
		if (mag_diff > 0) {
            // Calculate the euclidean distance between previous frame & expected current frame
            value = std::sqrt(sqr(odf->magSpec[i]) + sqr(odf->prevMagSpec[i]) - 2 * odf->magSpec[i] * odf->prevMagSpec[i] * std::cos(pdev));

			// add to sum
			sum = sum + value;
		}

		// store values for next calculation
		odf->prevPhase2[i] = odf->prevPhase[i];
		odf->prevPhase[i] = odf->phase[i];
		odf->prevMagSpec[i] = odf->magSpec[i];
	}

	return sum;
}

float complexSpectralDifferenceSqHWR(struct odf * odf) {
	float dev,pdev;
	float sum = 0;
	float mag_diff;
	float value;

	// perform the FFT
	performFFT(odf);

	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// phase deviation
		dev = odf->phase[i] - (2*odf->prevPhase[i]) + odf->prevPhase2[i];
		// wrap into [-pi,pi] range
		pdev = princarg(dev);
		// calculate magnitude difference (real part of Euclidean distance between complex frames)
		mag_diff = odf->magSpec[i] - odf->prevMagSpec[i];
		// if we have a positive change in magnitude, then include in sum, otherwise ignore (half-wave rectification)
		if (mag_diff > 0) {
            // Calculate the euclidean distance squared between previous frame & expected current frame
            value = (sqr(odf->magSpec[i]) + sqr(odf->prevMagSpec[i]) - 2 * odf->magSpec[i] * odf->prevMagSpec[i] * std::cos(pdev));
			// add to sum
			sum = sum + value;
		}
		// store values for next calculation
		odf->prevPhase2[i] = odf->prevPhase[i];
		odf->prevPhase[i] = odf->phase[i];
		odf->prevMagSpec[i] = odf->magSpec[i];
	}

	return sum;
}

float highFrequencyContent(struct odf * odf) {
	float sum = 0;

	// perform the FFT
	performFFT(odf);

	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		sum = sum + (odf->magSpec[i]*((float) (i+1)));
		// store values for next calculation
		odf->prevMagSpec[i] = odf->magSpec[i];
	}

	return sum;
}

static float highFrequencySpectralDifference(struct odf * odf) {
	float sum = 0;
	float mag_diff;

	// perform the FFT
	performFFT(odf);

	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// calculate difference
		mag_diff = std::abs(odf->magSpec[i] - odf->prevMagSpec[i]);
		sum = sum + (mag_diff*((float) (i+1)));
		// store values for next calculation
		odf->prevMagSpec[i] = odf->magSpec[i];
	}

	return sum;
}

static float highFrequencySpectralDifferenceHWR(struct odf * odf) {
	float sum = 0;
	float mag_diff;

	// perform the FFT
	performFFT(odf);

	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// calculate difference
		mag_diff = odf->magSpec[i] - odf->prevMagSpec[i];
		if (mag_diff > 0) {
            sum = sum + (mag_diff*((float ) (i+1)));
		}

		// store values for next calculation
		odf->prevMagSpec[i] = odf->magSpec[i];
	}

	return sum;
}

////////////////////////////// Methods to Calculate Windows ////////////////////////////////////

static void calculateHanningWindow(float * window, int frameSize){
	float N = (float) (frameSize-1);	// framesize minus 1

	// Hanning window calculation
	for (int n = 0;n < frameSize;n++) {
		window[n] = 0.5*(1-cos(2*M_PI*(n/N)));
	}
}

static void calclulateHammingWindow(float * window, int frameSize) {
	float N = (float) (frameSize-1);	// framesize minus 1
	float n_val = 0;// float version of index 'n'

	// Hamming window calculation
	for (int n = 0;n < frameSize;n++) {
		window[n] = 0.54 - (0.46*cos(2*M_PI*(n_val/N)));
		n_val = n_val+1;
	}
}

static void calculateBlackmanWindow(float * window, int frameSize) {
	float N = (float) (frameSize-1);	// framesize minus 1
	float n_val = 0;	// float version of index 'n'

	// Blackman window calculation
	for (int n = 0;n < frameSize;n++) {
		window[n] = 0.42 - (0.5*cos(2*M_PI*(n_val/N))) + (0.08*cos(4*M_PI*(n_val/N)));
		n_val = n_val+1;
	}
}

static void calculateTukeyWindow(float * window, int frameSize) {
	float alpha = 0.5;	// alpha [default value = 0.5];
	float N = (float) (frameSize-1);	// framesize minus 1
	float n_val = (float) (-1*((frameSize/2)))+1;

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
static void calculateRectangularWindow(float * window, int frameSize) {
	for (int n = 0;n < frameSize;n++) {
		window[n] = 1.0;
	}
}
///////////////////////////////// Other Handy Methods //////////////////////////////////////////
