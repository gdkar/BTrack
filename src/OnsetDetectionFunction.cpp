//=======================================================================
/** @file OnsetDetectionFunction.cpp
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

#include "OnsetDetectionFunction.h"
#include <cstring>

//=======================================================================
//=======================================================================
OnsetDetectionFunction::OnsetDetectionFunction(
    int hopSize_
    ,int frameSize_
    ,int onsetDetectionFunctionType_
    ,int windowType_)
: frameSize(frameSize_)
, hopSize(hopSize_)
, onsetDetectionFunctionType(ComplexSpectralDifferenceHWR)
, windowType(HanningWindow)
, p(nullptr)
, realIn((float*)fftwf_malloc(frameSize_*sizeof(float)))
, imagIn((float*)fftwf_malloc(frameSize_*sizeof(float)))
, realOut((float*)fftwf_malloc(frameSize_*sizeof(float)))
, imagOut((float*)fftwf_malloc(frameSize_*sizeof(float)))
, frame(frameSize_)
, window(frameSize_)
, prevEnergySum(_mm_setzero_ps())
, magSpec(frameSize_)
, prevMagSpec(frameSize_)
, phase(frameSize_)
, prevPhase(frameSize_)
, prevPhase2(frameSize_)
{	
  fftwf_iodim dimensions{ frameSize_, 1, 1 };
	p = fftwf_plan_guru_split_dft(
      1               /* rank */
      ,&dimensions    /* strides */
      , 0             /* howmany_rank */
      , 0             /* howmany_dims */
      ,realIn        /* ri */
      ,imagIn        /* ii */
      ,realOut        /* ro */
      ,imagOut        /* io */
      ,FFTW_ESTIMATE);

	switch (windowType){
		case RectangularWindow:
			calculateRectangularWindow();		// Rectangular window
			break;	
		case HanningWindow:
			calculateHanningWindow();			// Hanning Window
			break;
		case HammingWindow:
			calclulateHammingWindow();			// Hamming Window
			break;
		case BlackmanWindow:
			calculateBlackmanWindow();			// Blackman Window
			break;
		case TukeyWindow:
			calculateTukeyWindow();             // Tukey Window
			break;
		default:
			calculateHanningWindow();			// DEFAULT: Hanning Window
	}
}


//=======================================================================
OnsetDetectionFunction::~OnsetDetectionFunction()
{
    // destroy fft plan
    if(realIn)fftwf_free(realIn);
    if(imagIn)fftwf_free(imagIn);
    if(realOut)fftwf_free(realOut);
    if(imagOut)fftwf_free(imagOut);
    if(p)fftwf_destroy_plan(p);
}


//=======================================================================
void OnsetDetectionFunction :: setOnsetDetectionFunctionType(int odt){onsetDetectionFunctionType = odt;}
//=======================================================================
float OnsetDetectionFunction :: calculateOnsetDetectionFunctionSample(float *buffer){	
	// shift audio samples back in frame by hop size
  memmove(&frame[0],&frame[hopSize],(frameSize-hopSize)*sizeof(float));
  memmove(&frame[frameSize-hopSize],buffer,hopSize);
		
	switch (onsetDetectionFunctionType){
		case EnergyEnvelope:
        {
            // calculate energy envelope detection function sample
			return energyEnvelope();
			break;
        }
		case EnergyDifference:
        {
            // calculate half-wave rectified energy difference detection function sample
			return energyDifference();
			break;
        }
		case SpectralDifference:
        {
            // calculate spectral difference detection function sample
			return spectralDifference();
			break;
        }
		case SpectralDifferenceHWR:
        {
            // calculate spectral difference detection function sample (half wave rectified)
			return spectralDifferenceHWR();
			break;
        }
		case PhaseDeviation:
        {
            // calculate phase deviation detection function sample (half wave rectified)
			return phaseDeviation();
			break;
        }
		case ComplexSpectralDifference:
        {
            // calcualte complex spectral difference detection function sample
			return complexSpectralDifference();
			break;
        }
		case ComplexSpectralDifferenceHWR:
        {
            // calcualte complex spectral difference detection function sample (half-wave rectified)
			return complexSpectralDifferenceHWR();
			break;
        }
		case HighFrequencyContent:
        {
            // calculate high frequency content detection function sample
			return highFrequencyContent();
			break;
        }
		case HighFrequencySpectralDifference:
        {
            // calculate high frequency spectral difference detection function sample
			return highFrequencySpectralDifference();
			break;
        }
		case HighFrequencySpectralDifferenceHWR:
        {
            // calculate high frequency spectral difference detection function (half-wave rectified)
			return highFrequencySpectralDifferenceHWR();
			break;
        }
		default:
        {
			return 1.0;
        }
	}
		
}


//=======================================================================
void OnsetDetectionFunction :: performFFT()
{
	const int fsize2 = (frameSize/2);
	for (int i = 0;i < fsize2;i+=4){
    *(v4sf*)(realIn+i)        = _mm_mul_ps(*(v4sf*)(&frame[fsize2+i]),*(v4sf*)(&window[fsize2+i]));
    *(v4sf*)(realIn+i+fsize2) = _mm_mul_ps(*(v4sf*)(&frame[i])       ,*(v4sf*)(&window[i]));
	}
  std::memset(imagIn,0,frameSize*sizeof(float));
	// perform the fft
	fftwf_execute(p);
}

void OnsetDetectionFunction :: toPolar()
{
  memmove(&prevPhase2[0],&prevPhase[0],frameSize*sizeof(float));
  memmove(&prevPhase[0], &phase[0],frameSize*sizeof(float));
  memmove(&prevMagSpec[0],&magSpec[0],frameSize*sizeof(float));
  for(int i = 0; i < frameSize; i+= 4)
  {
    _approx_magphase_ps ((v4sf*) &magSpec[i], (v4sf*)&phase[i], *(v4sf*)(realOut+i),*(v4sf*)(imagOut+i));
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Methods for Detection Functions /////////////////////////////////

//=======================================================================
float OnsetDetectionFunction :: energyEnvelope()
{
	v4sf sum = _mm_setzero_ps();
	// sum the squares of the samples
	for (int i = 0;i < frameSize;i+=4)
	{
    sum = _mm_add_ps(sum,_mm_mul_ps(*(v4sf*)&frame[i],*(v4sf*)&frame[i]));
	}
  sum = _mm_hadd_ps(sum,sum);	
  sum = _mm_hadd_ps(sum,sum);	
	return sum[0];		// return sum
}

//=======================================================================
float OnsetDetectionFunction :: energyDifference()
{
  v4sf sum = _mm_setzero_ps();
  for ( int i = 0; i < frameSize; i += 4 )
  {sum = _mm_add_ps(sum,_mm_mul_ps(*(v4sf*)&frame[i],*(v4sf*)&frame[i]));}
  sum = _mm_hadd_ps(sum,sum);	
  sum = _mm_hadd_ps(sum,sum);	
	v4sf sample = _mm_sub_ps(sum,prevEnergySum);
  prevEnergySum = sum;
  return _mm_max_ss(sample,_mm_setzero_ps())[0];		// return difference
}

//=======================================================================
float OnsetDetectionFunction :: spectralDifference()
{
	// perform the FFT
	performFFT();
  toPolar();	
	v4sf sum = _mm_setzero_ps();
  for ( int i = 0; i < frameSize; i += 4 )
  {sum = _mm_add_ps(sum,_mm_abs_ps(_mm_sub_ps(*(v4sf*)&magSpec[i],*(v4sf*)&prevMagSpec[i])));}
  sum = _mm_hadd_ps(sum,sum);
  sum = _mm_hadd_ps(sum,sum);
	return sum[0];
}
//=======================================================================
float OnsetDetectionFunction :: spectralDifferenceHWR()
{
	
	// perform the FFT
	performFFT();
  toPolar();	
	v4sf sum = _mm_setzero_ps();
  for ( int i = 0; i < frameSize; i+=4 )
  {
    sum = _mm_add_ps(sum,_mm_max_ps(_mm_sub_ps(*(v4sf*)&magSpec[i],*(v4sf*)&prevMagSpec[i]),_mm_setzero_ps()));
  }
  sum = _mm_hadd_ps(sum,sum);
  sum = _mm_hadd_ps(sum,sum);
	return sum[0];
}


//=======================================================================
float OnsetDetectionFunction :: phaseDeviation()
{
	// perform the FFT
	performFFT();
  toPolar();	
	v4sf sum = _mm_setzero_ps();
	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < frameSize;i+=4){
    v4sf _pp = *(v4sf*)&prevPhase[i];
    v4sf dev = _mm_add_ps(_mm_sub_ps(*(v4sf*)&phase[i],_pp),
                     _mm_sub_ps(*(v4sf*)&prevPhase2[i],_pp));
    v4sf mask = _mm_cmpgt_ps(*(v4sf*)&magSpec[i],_mm_set1_ps(.1));
    sum = _mm_add_ps(sum,_mm_and_ps(mask,_mm_abs_ps(dev)));
	}
  sum = _mm_hadd_ps(sum,sum);sum=_mm_hadd_ps(sum,sum);	
	return sum[0];
}

//=======================================================================
float OnsetDetectionFunction :: complexSpectralDifference()
{
	// perform the FFT
	performFFT();
  toPolar();	
  v4sf sum = _mm_setzero_ps();
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < frameSize;i+=4)
	{

    v4sf _pp = *(v4sf*)&prevPhase[i];
    v4sf dev = _mm_princarg_ps(_mm_add_ps(_mm_sub_ps(*(v4sf*)&phase[i],_pp),
                     _mm_sub_ps(*(v4sf*)&prevPhase2[i],_pp)));

	  v4sf mag_diff = _mm_sub_ps(*(v4sf*)&magSpec[i],*(v4sf*)&prevMagSpec[i]);
	  	
		// calculate phase difference (imaginary part of Euclidean distance between complex frames)
		v4sf phase_diff = _mm_mul_ps(*(v4sf*)&magSpec[i],sin_ps(dev));
	  sum = _mm_add_ps(sum,_mm_hypot_ps(phase_diff,mag_diff));	
	}
	
	sum = _mm_hadd_ps(sum,sum);sum=_mm_hadd_ps(sum,sum);	
	return sum[0];
}

//=======================================================================
float OnsetDetectionFunction :: complexSpectralDifferenceHWR()
{
	// perform the FFT
	performFFT();
  toPolar();
  v4sf sum = _mm_setzero_ps();	
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < frameSize;i+=4){
		v4sf _pp = *(v4sf*)&prevPhase[i];
    v4sf dev = _mm_princarg_ps(_mm_add_ps(_mm_sub_ps(*(v4sf*)&phase[i],_pp),
                     _mm_sub_ps(*(v4sf*)&prevPhase2[i],_pp)));

	  v4sf mag_diff = _mm_sub_ps(*(v4sf*)&magSpec[i],*(v4sf*)&prevMagSpec[i]);
    v4sf phase_diff = _mm_mul_ps(*(v4sf*)&magSpec[i],sin_ps(dev));
    v4sf mask = _mm_cmpgt_ps(mag_diff,_mm_setzero_ps());
	  sum = _mm_add_ps(sum,_mm_and_ps(mask,_mm_hypot_ps(phase_diff,mag_diff)));
	}
	
	sum = _mm_hadd_ps(sum,sum);sum=_mm_hadd_ps(sum,sum);	
	return sum[0];
}


//=======================================================================
float OnsetDetectionFunction :: highFrequencyContent()
{
	performFFT();
  toPolar();	
  v4sf sum = _mm_setzero_ps();	
  v4sf inc = _mm_set_ps(4,3,2,0);
  const v4sf adv = _mm_set1_ps(1);
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < frameSize;i+=4){		
		sum = _mm_add_ps(sum,_mm_mul_ps(inc,*(v4sf*)&magSpec[i]));
    inc = _mm_add_ps(inc,adv);
	}
  sum = _mm_hadd_ps(sum,sum);sum=_mm_hadd_ps(sum,sum);	
	return sum[0];
}

//=======================================================================
float OnsetDetectionFunction :: highFrequencySpectralDifference()
{
	performFFT();
  toPolar();	
  v4sf sum = _mm_setzero_ps();	
  v4sf inc = _mm_set_ps(4,3,2,0);
  const v4sf adv = _mm_set1_ps(1);
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < frameSize;i+=4)
	{		
    sum = _mm_add_ps(sum,_mm_abs_ps(_mm_mul_ps(inc,_mm_sub_ps(*(v4sf*)&magSpec[i],
                                                          *(v4sf*)&prevMagSpec[i]))));
    inc = _mm_add_ps(inc,adv);
	}
  sum = _mm_hadd_ps(sum,sum);sum=_mm_hadd_ps(sum,sum);	
	return sum[0];
}

//=======================================================================
float OnsetDetectionFunction :: highFrequencySpectralDifferenceHWR()
{
	performFFT();
  toPolar();	
  v4sf sum = _mm_setzero_ps();	
  v4sf inc = _mm_set_ps(4,3,2,0);
  const v4sf adv = _mm_set1_ps(1);
	// compute phase values from fft output and sum deviations
	for (int i = 0;i < frameSize;i+=4)
	{		
    sum = _mm_add_ps(sum,_mm_max_ps(_mm_setzero_ps(),_mm_mul_ps(inc,_mm_sub_ps(*(v4sf*)&magSpec[i],
                                                                           *(v4sf*)&prevMagSpec[i]))));
    inc = _mm_add_ps(inc,adv);
	}
  sum = _mm_hadd_ps(sum,sum);sum=_mm_hadd_ps(sum,sum);	
	return sum[0];
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Methods to Calculate Windows ////////////////////////////////////

//=======================================================================
void OnsetDetectionFunction :: calculateHanningWindow()
{
	float N;		// variable to store framesize minus 1
	
	N = (float) (frameSize-1);	// framesize minus 1
	
	// Hanning window calculation
	for (int n = 0;n < frameSize;n++)
	{
		window[n] = 0.5*(1-cos(2*M_PI*(n/N)));
	}
}

//=======================================================================
void OnsetDetectionFunction :: calclulateHammingWindow()
{
	float N;		// variable to store framesize minus 1
	float n_val;	// float version of index 'n'
	
	N = (float) (frameSize-1);	// framesize minus 1
	n_val = 0;
	
	// Hamming window calculation
	for (int n = 0;n < frameSize;n++)
	{
		window[n] = 0.54 - (0.46*cos(2*M_PI*(n_val/N)));
		n_val = n_val+1;
	}
}

//=======================================================================
void OnsetDetectionFunction :: calculateBlackmanWindow()
{
	float N;		// variable to store framesize minus 1
	float n_val;	// float version of index 'n'
	
	N = (float) (frameSize-1);	// framesize minus 1
	n_val = 0;
	
	// Blackman window calculation
	for (int n = 0;n < frameSize;n++)
	{
		window[n] = 0.42 - (0.5*cos(2*M_PI*(n_val/N))) + (0.08*cos(4*M_PI*(n_val/N)));
		n_val = n_val+1;
	}
}

//=======================================================================
void OnsetDetectionFunction :: calculateTukeyWindow()
{
	float N;		// variable to store framesize minus 1
	float n_val;	// float version of index 'n'
	float alpha;	// alpha [default value = 0.5];
	
	alpha = 0.5;
	
	N = (float) (frameSize-1);	// framesize minus 1
		
	// Tukey window calculation
	
	n_val = (float) (-1*((frameSize/2)))+1;

	for (int n = 0;n < frameSize;n++)	// left taper
	{
		if ((n_val >= 0) && (n_val <= (alpha*(N/2))))
		{
			window[n] = 1.0;
		}
		else if ((n_val <= 0) && (n_val >= (-1*alpha*(N/2))))
		{
			window[n] = 1.0;
		}else{
			window[n] = 0.5*(1+cos(M_PI*(((2*n_val)/(alpha*N))-1)));
		}
		n_val = n_val+1;			 
	}

}

//=======================================================================
void OnsetDetectionFunction :: calculateRectangularWindow()
{
	// Rectangular window calculation
	for (int n = 0;n < frameSize;n++){window[n] = 1.0;}
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Other Handy Methods //////////////////////////////////////////

//=======================================================================
float OnsetDetectionFunction :: princarg(float phaseVal)
{	return phaseVal - (float)(2*M_PI) * floor(phaseVal/(float)(2*M_PI));}
v4sf OnsetDetectionFunction::_mm_princarg_ps(v4sf ph)
{
  v4sf mask_neg = _mm_and_ps(_mm_cmpgt_ps(ph,*(v4sf*)_ps_cephes_PI),*(v4sf*)_ps_cephes_TWOPI);
  v4sf mask_pos = _mm_and_ps(_mm_cmplt_ps(ph,*(v4sf*)_ps_cephes_minus_PI),*(v4sf*)_ps_cephes_TWOPI);
  return _mm_add_ps(mask_pos,_mm_sub_ps(ph,mask_neg));
}
