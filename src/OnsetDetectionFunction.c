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
#include <math.h>

static void performFFT(struct odf * odf);

static float energyEnvelope(struct odf * odf);
static float energyDifference(struct odf * odf);
static float spectralDifference(struct odf * odf);
static float spectralDifferenceHWR(struct odf * odf);
static float phaseDeviation(struct odf * odf);
static float complexSpectralDifference(struct odf * odf);
static float complexSpectralDifferenceHWR(struct odf * odf);
static float highFrequencyContent(struct odf * odf);
static float highFrequencySpectralDifference(struct odf * odf);
static float highFrequencySpectralDifferenceHWR(struct odf * odf);
static inline float __attribute__((unused)) princarg(float arg)
{
  static const float one_over_two_pi = 1/(2*M_PI);
  return arg - (2*M_PI)*floorf(arg*one_over_two_pi);
}
static inline __m128 princarg_ps(const __m128);
static void calculateRectangularWindow(float * window, int frameSize);	
static void calculateHanningWindow(float * window, int frameSize);		
static void calclulateHammingWindow(float * window, int frameSize);		
static void calculateBlackmanWindow(float * window, int frameSize);		
static void calculateTukeyWindow(float * window, int frameSize);          

int odf_init(struct odf * odf, int hop_size, int frame_size, enum OnsetDetectionFunctionType odf_type, enum WindowType window_type){
    // if we have already initialised FFT plan
    odf_del(odf);
	
	odf->hopSize = hop_size; // set hopsize
	odf->frameSize = frame_size; // set framesize
	
    odf->type = odf_type;
    odf->windowType = window_type;
		
	// initialise buffers
    odf->frame = malloc(sizeof(float) * frame_size);
    odf->window = malloc(sizeof(float) * frame_size);
    odf->magSpec = malloc(sizeof(float) * frame_size);
    odf->prevMagSpec = malloc(sizeof(float) * frame_size);
    odf->phase = malloc(sizeof(float) * frame_size);
    odf->prevPhase = malloc(sizeof(float) * frame_size);
    odf->prevPhase2 = malloc(sizeof(float) * frame_size);
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
  memset(odf->prevMagSpec,0,frame_size*sizeof(float));
  memset(odf->prevPhase,0,frame_size*sizeof(float));
  memset(odf->prevPhase2,0,frame_size*sizeof(float));
  memset(odf->frame,0,frame_size*sizeof(float));
	odf->prevEnergySum = 0.0;	// initialise previous energy sum value to zero
	/*  Init fft */
	odf->realIn= (float*) fftwf_malloc(sizeof(float) * frame_size);		// complex array to hold fft data
	odf->imagIn= (float*) fftwf_malloc(sizeof(float) * frame_size);		// complex array to hold fft data
	odf->realOut = (float*) fftwf_malloc(sizeof(float) * frame_size);	// complex array to hold fft data
	odf->imagOut = (float*) fftwf_malloc(sizeof(float) * frame_size);	// complex array to hold fft data
  fftw_iodim transform_dims[] = {{ .n = frame_size, .is = 1, .os=1}};
  odf->p = fftwf_plan_guru_split_dft(1, transform_dims, 0, NULL,
  odf->realIn,odf->imagIn,odf->realOut,odf->imagOut,FFTW_ESTIMATE);
	
	odf->initialised = true;
    return 0;
}

void odf_del(struct odf * odf){
    if (odf->initialised) {
        // destroy fft plan
        fftwf_destroy_plan(odf->p);
        fftwf_free(odf->realIn);
        fftwf_free(odf->imagIn);
        fftwf_free(odf->realOut);
        fftwf_free(odf->imagOut);
    }
}

void odf_set_type(struct odf * odf, enum OnsetDetectionFunctionType type){
    odf->type = type;
}

float odf_calculate_sample(struct odf * odf, float * buffer){
	float odfSample;
	
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
  //
  for(int i = 0; i < fsize2; i+=4)
  {
    *(__m128*)(odf->realIn+fsize2+i) = 
      _mm_mul_ps(*(__m128*)(odf->frame+i),*(__m128*)(odf->window+i));
    *(__m128*)(odf->realIn+i) = 
      _mm_mul_ps(*(__m128*)(odf->frame+fsize2+i),*(__m128*)(odf->window+fsize2+i));

  }
  memset(odf->imagIn,0,odf->frameSize*sizeof(float));
	// perform the fft
	fftwf_execute(odf->p);
}

////////////////////////////// Methods for Detection Functions /////////////////////////////////

static float energyEnvelope(struct odf * odf){
	// sum the squares of the samples
  __m128 accum = _mm_setzero_ps();
  for(int i = 0; i < odf->frameSize; i+=4)
  {
    __m128 these_4 = *(__m128 *)(odf->frame+i);
    accum = _mm_add_ps(accum,_mm_mul_ps(these_4,these_4));
  }
  const __m128 accum1 = _mm_add_ps(accum,_mm_movehl_ps(accum,accum));
  const __m128 accum2 = _mm_add_ss(accum1,_mm_shuffle_ps(accum1,accum1,1));
	return accum2[0];		// return sum
}
static void genMagnitudeSpectrum(struct odf * odf ){
	performFFT(odf);
  for ( int i = 0; i < odf->frameSize; i+=4)
  {
    v4sf mag,phase;
    *(__m128*)(odf->prevMagSpec+i)=*(__m128*)(odf->magSpec+i);
    *(__m128*)(odf->prevPhase2+i)=*(__m128*)(odf->prevPhase+i);
    *(__m128*)(odf->prevPhase+i)=*(__m128*)(odf->phase+i);
    _approx_magphase_ps(&mag,&phase,*(__m128*)(odf->realOut+i),*(__m128*)(odf->imagOut+i));
    *(__m128*)(odf->magSpec+i)=mag;
    *(__m128*)(odf->phase+i)=phase;
  }
}
static float energyDifference(struct odf * odf){
  __m128 sum = _mm_setzero_ps();	
	// sum the squares of the samples
	for (int i = 0;i < odf->frameSize;i+=4) {
    __m128 these_4 = *(__m128*)(odf->frame+i);
    sum = _mm_add_ps(sum,_mm_mul_ps(these_4,these_4));
	}
  const __m128 accum1 = _mm_add_ps(sum,_mm_movehl_ps(sum,sum));
  const __m128 accum2 = _mm_add_ss(accum1,_mm_shuffle_ps(accum1,accum1,1));
	

	float sample = accum2[0] - odf->prevEnergySum;	// sample is first order difference in energy
  odf->prevEnergySum = accum2[0];
  return MAX(sample, 0);
}

static float spectralDifference(struct odf * odf){
	__m128 sum =_mm_setzero_ps();
	// perform the FFT
  genMagnitudeSpectrum(odf);
	// compute first (N/2)+1 mag values
	for (int i = 0;i < odf->frameSize;i+=4) {
    __m128 these_4 = *(__m128 *)(odf->magSpec+i);
    sum = _mm_add_ps(sum,_mm_abs_ps(_mm_sub_ps(*(__m128*)(odf->prevMagSpec+i),these_4)));
    *(__m128 *)(odf->prevMagSpec+i) = these_4;
	}
  const __m128 accum1 = _mm_add_ps(sum,_mm_movehl_ps(sum,sum));
  const __m128 accum2 = _mm_add_ss(accum1,_mm_shuffle_ps(accum1,accum1,1));
	return accum2[0];		// return sum

}

static float spectralDifferenceHWR(struct odf * odf) {
	__m128 sum = _mm_setzero_ps();

	// perform the FFT
  genMagnitudeSpectrum(odf);
	for (int i = 0;i < odf->frameSize;i+=4) {
    __m128 these_4 = *(__m128*)(odf->magSpec+i);
    sum = _mm_add_ps(sum,_mm_max_ps(_mm_sub_ps(these_4,*(__m128*)(odf->prevMagSpec+i)),_mm_setzero_ps()));
    *(__m128*)(odf->prevMagSpec+i)=these_4;
  }
	return sum[0];
}

static float phaseDeviation(struct odf * odf){
	__m128 sum = _mm_setzero_ps();
	
	// perform the FFT
  genMagnitudeSpectrum(odf);	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// calculate phase value
		
    __m128 mag_diff = _mm_sub_ps(*(__m128*)(odf->magSpec+i),*(__m128*)(odf->prevMagSpec+i));
    __m128 mask = _mm_cmpgt_ps(mag_diff,_mm_set1_ps(0.0));
      __m128 dev = 
      _mm_sub_ps(_mm_add_ps(*(__m128*)(odf->phase+i),*(__m128*)(odf->prevPhase2+i)),
      _mm_add_ps(*(__m128*)(odf->prevPhase+i),*(__m128*)(odf->prevPhase+i)));
		// wrap into [-pi,pi] range
		__m128 pdev = princarg_ps(dev);	
    sum = _mm_add_ps(sum,_mm_abs_ps(_mm_and_ps(mask,pdev)));
  }
  const __m128 accum1 = _mm_add_ps(sum,_mm_movehl_ps(sum,sum));
  const __m128 accum2 = _mm_add_ss(accum1,_mm_shuffle_ps(accum1,accum1,1));
	return accum2[0];		// return sum

}

float complexSpectralDifference(struct odf * odf){
	__m128 sum = _mm_setzero_ps();
	// perform the FFT
	genMagnitudeSpectrum(odf);
	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {
		// calculate phase value
				// phase deviation
		 __m128 dev = 
      _mm_sub_ps(_mm_add_ps(*(__m128*)(odf->phase+i),*(__m128*)(odf->prevPhase2+i)),
      _mm_add_ps(*(__m128*)(odf->prevPhase+i),*(__m128*)(odf->prevPhase+i)));
		// wrap into [-pi,pi] range
		__m128 pdev = princarg_ps(dev);	
    __m128 mag_diff = _mm_sub_ps(*(__m128*)(odf->magSpec+i),*(__m128*)(odf->prevMagSpec+i));
		
		// if we have a positive change in magnitude, then include in sum, otherwise ignore (half-wave rectification)
//    __m128 mag_diff_mask = _mm_cmpgt_ps(mag_diff,_mm_setzero_ps());
    __m128 phase_diff = _mm_mul_ps(*(__m128*)(odf->magSpec+i),sin_ps(pdev));
		
		// square real and imaginary parts, sum and take square root
		sum =_mm_add_ps(sum, _mm_hypot_ps(mag_diff,phase_diff));
	
		// store values for next calculation
	}
	
	const __m128 accum1 = _mm_add_ps(sum,_mm_movehl_ps(sum,sum));
  const __m128 accum2 = _mm_add_ss(accum1,_mm_shuffle_ps(accum1,accum1,1));
	return accum2[0];		// return sum
}

float complexSpectralDifferenceHWR(struct odf * odf) {
	__m128 sum = _mm_setzero_ps();
	genMagnitudeSpectrum(odf);
	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i+=4) {
    __m128 dev = 
      _mm_sub_ps(_mm_add_ps(*(__m128*)(odf->phase+i),*(__m128*)(odf->prevPhase2+i)),
      _mm_add_ps(*(__m128*)(odf->prevPhase+i),*(__m128*)(odf->prevPhase+i)));
		// wrap into [-pi,pi] range
		__m128 pdev = princarg_ps(dev);	
		
		// calculate magnitude difference (real part of Euclidean distance between complex frames)
		__m128 mag_diff = _mm_sub_ps(*(__m128*)(odf->magSpec+i),*(__m128*)(odf->prevMagSpec+i));
		
		// if we have a positive change in magnitude, then include in sum, otherwise ignore (half-wave rectification)
    __m128 mag_diff_mask = _mm_cmpgt_ps(mag_diff,_mm_setzero_ps());
    __m128 phase_diff = _mm_mul_ps(*(__m128*)(odf->magSpec+i),
        sin_ps(pdev));
    __m128 _hypot = _mm_hypot_ps(mag_diff,phase_diff);
    sum = _mm_add_ps(sum,_mm_and_ps(mag_diff_mask,_hypot));
		
		// store values for next calculation
		odf->prevPhase2[i] = odf->prevPhase[i];
		odf->prevPhase[i] = odf->phase[i];
		odf->prevMagSpec[i] = odf->magSpec[i];
	}
	  const __m128 accum1 = _mm_add_ps(sum,_mm_movehl_ps(sum,sum));
  const __m128 accum2 = _mm_add_ss(accum1,_mm_shuffle_ps(accum1,accum1,1));
	return accum2[0];		// return sum
}

float highFrequencyContent(struct odf * odf) {
	// perform the FFT
	genMagnitudeSpectrum(odf);
	__m128 sum = _mm_setzero_ps();
  __m128 counter = _mm_set_ps(1,2,3,4);
  const __m128 count = _mm_set1_ps(1);
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {		
		sum = _mm_add_ps(sum,_mm_mul_ps(counter,*(__m128*)(odf->magSpec+i)));
    counter = _mm_add_ps(counter,count);
	}
	const __m128 accum1 = _mm_add_ps(sum,_mm_movehl_ps(sum,sum));
  const __m128 accum2 = _mm_add_ss(accum1,_mm_shuffle_ps(accum1,accum1,1));
	return accum2[0];		// return sum
}

static float highFrequencySpectralDifference(struct odf * odf) {
	// perform the FFT
	genMagnitudeSpectrum(odf);

  __m128 sum = _mm_setzero_ps();
  __m128 counter = _mm_set_ps(1,2,3,4);
  const __m128 count = _mm_set1_ps(1);
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {		
    __m128 diff = _mm_abs_ps(_mm_sub_ps(*(__m128*)(odf->magSpec+i),
          *(__m128*)(odf->prevMagSpec+i)));
		sum = _mm_add_ps(sum,_mm_mul_ps(counter,diff));
    counter = _mm_add_ps(counter,count);
	}
	const __m128 accum1 = _mm_add_ps(sum,_mm_movehl_ps(sum,sum));
  const __m128 accum2 = _mm_add_ss(accum1,_mm_shuffle_ps(accum1,accum1,1));
	return accum2[0];		// return sum
}

static float highFrequencySpectralDifferenceHWR(struct odf * odf) {
	// perform the FFT
	genMagnitudeSpectrum(odf);
	
  __m128 sum = _mm_setzero_ps();
  __m128 counter = _mm_set_ps(1,2,3,4);
  const __m128 count = _mm_set1_ps(1);
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < odf->frameSize;i++) {		
    __m128 diff = _mm_max_ps(_mm_sub_ps(*(__m128*)(odf->magSpec+i),
          *(__m128*)(odf->prevMagSpec+i)),_mm_setzero_ps());
		sum = _mm_add_ps(sum,_mm_mul_ps(counter,diff));
    counter = _mm_add_ps(counter,count);
	}
	const __m128 accum1 = _mm_add_ps(sum,_mm_movehl_ps(sum,sum));
  const __m128 accum2 = _mm_add_ss(accum1,_mm_shuffle_ps(accum1,accum1,1));
	return accum2[0];		// return sum
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

_PS_CONST(cephes_2PI,2*M_PI);
_PS_CONST(cephes_minus_2PI,-2*M_PI);
static __m128 princarg_ps(__m128 phaseVal) {	
	// if phase value is less than or equal to -pi then add 2*pi
  __m128 less_than_minus_pi_mask = _mm_cmplt_ps(
      phaseVal,
      *(const __m128*)_ps_cephes_minus_PI);
  __m128 greater_than_pi_mask = _mm_cmpgt_ps(
      phaseVal,
      *(const __m128*)_ps_cephes_PI);
  __m128 addend = _mm_or_ps(
              _mm_and_ps(
                greater_than_pi_mask,
                *(const __m128*)_ps_cephes_minus_2PI),
              _mm_and_ps(
                less_than_minus_pi_mask,
                *(const __m128*)_ps_cephes_2PI));
  return _mm_add_ps(phaseVal,addend);
}
