//=======================================================================
/** @file BTrack.c
 *  @brief BTrack - a real-time beat tracker
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

#include "BTrack.h"
#include "samplerate.h"
#include "alloca.h"
#include "common.h"
#include "Plan.hpp"
#include "ACF.hpp"

namespace retrack {
namespace detail {
template<class T>
constexpr T recip(T x)
{
    return T{1}/x;
}
template<class T>
constexpr T gaussian(T x, T sigma, T mu)
{
    return std::exp(-T(0.5) * sqr((x-mu)* recip(sigma)))* (recip(std::sqrt(2*T(M_PI)) * sigma));
}
template<class T>
constexpr T gaussian_unorm(T x, T sigma, T mu)
{
    return std::exp(-T(0.5) * sqr((x-mu)* recip(sigma)));
}

template<class T>
constexpr T gaussian_warp(T x, T sigma, T mu)
{
    return std::exp(-T(0.5) * sqr(std::log(std::abs(x) * recip(mu)) * recip(sigma)));
}
template<class T>
constexpr T rayleigh(T x, T sigma, T mu = T{})
{
    return x * std::exp(-T(0.5) * sqr((x-mu) * recip(sigma))) * recip(sqr(sigma));
}
template<class I>
void normalize_array(I ibeg, I iend)
{
    using T = typename std::iterator_traits<I>::value_type;
    auto acc = std::accumulate(ibeg,iend,T{});
    if(acc) {
        std::transform(ibeg,iend,ibeg,[n=recip(acc)](auto x){return x * n;});
    }
}
template<class I>
void adaptive_threshold(I ibeg, I iend, size_t _winlen = 15)
{
    auto ibeg_orig = ibeg, iend_orig = iend;
    using T = typename std::iterator_traits<I>::value_type;
    auto winlen  = ptrdiff_t(_winlen);
    auto winpre  = ptrdiff_t((winlen-1) >>1);
    auto buflen  = roundup(winlen);
    auto bufmask = buflen - 1ul;
    auto buf     = (T*)alloca(buflen * sizeof(T));
    std::fill_n(buf, buflen, T{});
    auto acc    = T{};
    auto wptr   = ptrdiff_t{};
    auto rptr   = ptrdiff_t{};
    auto extra  = ptrdiff_t{};
    auto _thresh= [](auto &v, auto t) {
        v = std::max<T>(T{},v-t);
    };
    for ( ; wptr < winpre && ibeg != iend;) {
        acc += (buf[wptr++&bufmask] = *ibeg++);
        extra++;
    }
    if ( ibeg == iend ) {
        std::for_each(ibeg_orig,iend_orig,
            [&_thresh,t=T(acc/wptr)](auto &v) {_thresh(v,t);}
            );
        return;
    }
    auto obeg = ibeg_orig;
    for ( ; ibeg != iend && wptr < winlen;) {
        acc += (buf[wptr++&bufmask] = *ibeg++);
        _thresh(*obeg++,T(acc/T(wptr)));
    }
    if ( ibeg == iend ) {
        for ( ; extra--;) {
            acc -= buf[rptr++ & bufmask];
            _thresh(*obeg++, T(acc /(wptr-rptr)));
        }
        return;
    }
    auto mult = T(1)/T(winlen);
    for ( ; ibeg != iend;) {
        acc -= (buf[rptr++&bufmask]);
        acc += (buf[wptr++&bufmask] = *ibeg++);
        _thresh(*obeg++, acc*mult);
    }
    for ( ; extra--; ) {
        acc -= (buf[rptr++&bufmask]);
        _thresh(*obeg++,T(acc/T(wptr-rptr)));
    }
}
}
}
using namespace retrack;
static void resampleOnsetDetectionFunction(struct btrack * bt);
static void calculateTempo(struct btrack * bt);
//static void adaptiveThreshold(float * x, int N, float * aux_buffer);
static void calculateOutputOfCombFilterBank(struct btrack * bt);
//static void calculateBalancedACF(float* onsetDetectionFunction, float * acf);
//static float calculateMeanOfArray(float * array, int startIndex, int endIndex);
static void updateCumulativeScore(struct btrack * bt, float odfSample);
static void predictBeat(struct btrack * bt);

// Ex. hop_size = 512; frame_size = 1024
btrack::btrack(int hop_size, int frame_size, int sample_rate)
: m_odf(hop_size, frame_size, ComplexSpectralDifferenceSq, HanningWindow)
, m_acf(512)
{
    hopSize = hop_size;
    frameSize = frame_size;
    invSampleRate = 1.0 / (double) sample_rate;

    auto rayparam = 43.;

	// initialise parameters
	tightness = 5;
	alpha = 0.9;
	tempo = 120;
	estimatedTempo = 120.0;
	tempoToLagFactor = 60.* (double) sample_rate / (double) hop_size;

	m0 = 10;
	beatCounter = -1;

	beatDueInFrame = false;

	// create rayleigh weighting vector
	for (int n = 0;n < 128;n++) {
		weightingVector[n] = detail::rayleigh(float(n),float(rayparam),float{});
	}

	// initialise prev_delta
	for (int i = 0;i < 41;i++) {
		prevDelta[i] = 1;
	}

	// create tempo transition matrix
	auto m_sig = 41/8.f;
	for (int i = 0;i < 41;i++) {
		for (int j = 0;j < 41;j++)
			tempoTransitionMatrix[i][j] = detail::gaussian(float(j+1),float(m_sig),float(i+1));
	}

	// tempo is not fixed
	tempoFixed = false;

    // initialise latest cumulative score value
    // in case it is requested before any processing takes place
    latestCumulativeScoreValue = 0;
    latestODF = 0;

    btrack_set_hop_size(this, hop_size);
}

int btrack_beat_due_in_current_frame(const struct btrack * bt){
    return bt->beatDueInFrame;
}

float btrack_get_bpm(const struct btrack * bt){
    return bt->estimatedTempo;
}

float btrack_get_latest_score(const struct btrack * bt){
    return bt->latestCumulativeScoreValue;
}

float btrack_get_latest_odf(const struct btrack * bt){
    return bt->latestODF;
}

float btrack_get_latest_confidence(const struct btrack * bt){
    return bt->latestConfidence;
}

int btrack_get_frames_until_beat(const struct btrack * bt) {
    return bt->beatCounter;
}

float btrack_get_time_until_beat(const struct btrack * bt) {
    auto n_frames = btrack_get_frames_until_beat(bt);
    return n_frames * bt->hopSize * bt->invSampleRate;
}

void btrack_process_audio_frame(struct btrack * bt, const btrack_chunk_t * frame){
    auto sample = bt->m_odf.process_frame(frame);
    btrack_process_odf_sample(bt, sample);
}

void btrack_process_fft_frame(struct btrack * bt, const btrack_chunk_t * frame){
    auto sample = bt->m_odf.process_fft_frame( frame);
    btrack_process_odf_sample(bt, sample);
}

void btrack_process_odf_sample(struct btrack * bt, float odf_sample){
    // we need to ensure that the onset
    // detection function sample is positive
    // add a tiny constant to the sample to stop it from ever going
    // to zero. this is to avoid problems further down the line
    odf_sample = fabs(odf_sample) + 0.0001;

	bt->m0--;
	bt->beatCounter--;
	bt->beatDueInFrame = false;

	// move all samples back one step
	for (int i=0;i < (bt->onsetDFBufferSize-1);i++)
	{
		bt->onsetDF[i] = bt->onsetDF[i+1];
	}

	// add new sample at the end
	bt->onsetDF[bt->onsetDFBufferSize-1] = odf_sample;
    bt->latestODF = odf_sample;

	// update cumulative score
	updateCumulativeScore(bt, odf_sample);

	// if we are halfway between beats
	if (bt->m0 == 0) {
		predictBeat(bt);
	}

	// if we are at a beat
	if (bt->beatCounter == 0) {
		bt->beatDueInFrame = true;	// indicate a beat should be output

		// recalculate the tempo
		resampleOnsetDetectionFunction(bt);
		calculateTempo(bt);
	}
}

void btrack_set_bpm(struct btrack * bt, float bpm){
	/////////// TEMPO INDICATION RESET //////////////////

	// firstly make sure tempo is between 80 and 160 bpm..
	while (bpm > 160) {
		bpm = bpm/2;
	}
	while (bpm < 80) {
		bpm = bpm * 2;
	}

	// convert tempo from bpm value to integer index of tempo probability
	int tempo_index = (int) round((bpm - 80)/2);

	// now set previous tempo observations to zero
	for (int i=0;i < 41;i++) {
		bt->prevDelta[i] = 0;
	}

	// set desired tempo index to 1
	bt->prevDelta[tempo_index] = 1;

	/////////// CUMULATIVE SCORE ARTIFICAL TEMPO UPDATE //////////////////

	// calculate new beat period
	int new_bperiod = (int) std::round(60/((((float) bt->hopSize)/44100)*bpm));

	int bcounter = 1;
	// initialise df_buffer to zeros
	for (int i = (bt->onsetDFBufferSize-1);i >= 0;i--) {
		if (bcounter == 1) {
			bt->cumulativeScore[i] = 150;
			bt->onsetDF[i] = 150;
		} else {
			bt->cumulativeScore[i] = 10;
			bt->onsetDF[i] = 10;
		}

		bcounter++;

		if (bcounter > new_bperiod) {
			bcounter = 1;
		}
	}

	/////////// INDICATE THAT THIS IS A BEAT //////////////////

	// beat is now
	bt->beatCounter = 0;

	// offbeat is half of new beat period away
	bt->m0 = (int) round(((float) new_bperiod)/2);
}

void btrack_fix_bpm(struct btrack * bt, float bpm){
	// firstly make sure tempo is between 80 and 160 bpm..
	while (bpm > 160) {
		bpm = bpm/2;
	}
	while (bpm < 80) {
		bpm = bpm * 2;
	}

	// convert tempo from bpm value to integer index of tempo probability
	int tempo_index = (int) round((bpm - 80)/2);

	// now set previous fixed previous tempo observation values to zero
	for (int i=0;i < 41;i++) {
		bt->prevDeltaFixed[i] = 0;
	}

	// set desired tempo index to 1
	bt->prevDeltaFixed[tempo_index] = 1;

	// set the tempo fix flag
	bt->tempoFixed = true;
}

void btrack_nofix_bpm(struct btrack * bt){
	// set the tempo fix flag
	bt->tempoFixed = false;
}

void btrack_set_hop_size(struct btrack * bt, int hop_size){
    bt->hopSize = hop_size;
	bt->onsetDFBufferSize = (512*512)/bt->hopSize;		// calculate df buffer size

    bt->onsetDF = std::make_unique<float[]>(bt->onsetDFBufferSize);
    BTRACK_ASSERT(bt->onsetDF);

    bt->m_acf = ACF(bt->onsetDFBufferSize);
    bt->cumulativeScore = std::make_unique<float[]>(bt->onsetDFBufferSize);
    BTRACK_ASSERT(bt->cumulativeScore);

    bt->w1 = std::make_unique<float[]>(bt->onsetDFBufferSize);
    BTRACK_ASSERT(bt->w1);

    bt->beatPeriod = std::round(60/((((float) bt->hopSize)/44100)*bt->tempo));
    // initialise df_buffer to zeros
    for (int i = 0;i < bt->onsetDFBufferSize;i++) {
        bt->onsetDF[i] = 0;
        bt->cumulativeScore[i] = 0;
        if ((i %  ((int) std::round(bt->beatPeriod))) == 0) {
            bt->onsetDF[i] = 1;
        }
    }
}

//=======================================================================

static void resampleOnsetDetectionFunction(struct btrack * bt) {
	float output[512];
    float * input = static_cast<float*>(alloca(bt->onsetDFBufferSize * sizeof(float)));
    BTRACK_ASSERT(input);

    for (int i = 0;i < bt->onsetDFBufferSize;i++) {
        input[i] = (float) bt->onsetDF[i];
    }

	auto src_ratio = 512.0/((double) bt->onsetDFBufferSize);
	int BUFFER_LEN = bt->onsetDFBufferSize;
	int output_len;
	SRC_DATA	src_data ;

	//output_len = (int) floor (((double) BUFFER_LEN) * src_ratio) ;
	output_len = 512;

	src_data.data_in = input;
	src_data.input_frames = BUFFER_LEN;

	src_data.src_ratio = src_ratio;

	src_data.data_out = output;
	src_data.output_frames = output_len;

	src_simple (&src_data, SRC_SINC_BEST_QUALITY, 1);

	for (int i = 0;i < output_len;i++) {
		bt->resampledOnsetDF[i] = (float) src_data.data_out[i];
	}
}

static void calculateTempo(struct btrack * bt){
	// adaptive threshold on input
	retrack::detail::adaptive_threshold(bt->resampledOnsetDF,bt->resampledOnsetDF + 512, 15);

	// calculate auto-correlation function of detection function
    bt->m_acf.process(bt->resampledOnsetDF, bt->acf);
	// calculate output of comb filterbank
	calculateOutputOfCombFilterBank(bt);
	// adaptive threshold on rcf
	retrack::detail::adaptive_threshold(bt->combFilterBankOutput,bt->combFilterBankOutput+128, 15);

	int t_index;
	int t_index2;
	// calculate tempo observation vector from beat period observation vector
	for (int i = 0;i < 41;i++) {
		t_index = (int) std::round(bt->tempoToLagFactor / ((float) ((2*i)+80)));
		t_index2 = (int) std::round(bt->tempoToLagFactor / ((float) ((4*i)+160)));

		bt->tempoObservationVector[i] = bt->combFilterBankOutput[t_index-1] + bt->combFilterBankOutput[t_index2-1];
	}

	float maxval;
	float maxind;
	float curval;

	// if tempo is fixed then always use a fixed set of tempi as the previous observation probability function
	if(bt->tempoFixed) {
		for (int k = 0;k < 41;k++) {
			bt->prevDelta[k] = bt->prevDeltaFixed[k];
		}
	}

	for(int j=0;j < 41;j++) {
		maxval = -1;
		for (int i = 0;i < 41;i++) {
			curval = bt->prevDelta[i]*bt->tempoTransitionMatrix[i][j];

			if (curval > maxval) {
				maxval = curval;
			}
		}

		bt->delta[j] = maxval*bt->tempoObservationVector[j];
	}

    retrack::detail::normalize_array(bt->delta,bt->delta+41);

	maxind = -1;
	maxval = -1;

	for (int j=0;j < 41;j++) {
		if (bt->delta[j] > maxval) {
			maxval = bt->delta[j];
			maxind = j;
		}
		bt->prevDelta[j] = bt->delta[j];
	}

	bt->beatPeriod = std::round((60.0*44100.0)/(((2*maxind)+80)*((double) bt->hopSize)));
	if (bt->beatPeriod > 0) {
		bt->estimatedTempo = 60.0/((((double) bt->hopSize) / 44100.0)*bt->beatPeriod);
	}
}
static void calculateOutputOfCombFilterBank(struct btrack * bt) {
	int numelem;

	for (int i = 0;i < 128;i++) {
		bt->combFilterBankOutput[i] = 0;
	}

	numelem = 4;

    // max beat period
	for (int i = 2;i <= 127;i++) {
        // number of comb elements
		for (int a = 1;a <= numelem;a++) {
            // general state using normalisation of comb elements
			for (int b = 1-a;b <= a-1;b++) {
				bt->combFilterBankOutput[i-1] = bt->combFilterBankOutput[i-1] + (bt->acf[(a*i+b)-1]*bt->weightingVector[i-1])/(2*a-1);	// calculate value for comb filter row
			}
		}
	}
}
static void updateCumulativeScore(struct btrack * bt, float odfSample) {
	float max, maxw, min;


	float wcumscore;


	// create window
    if(bt->beatPeriod != bt->w1size) {
        auto start = int(bt->onsetDFBufferSize - round(2*bt->beatPeriod));
        auto end = int(bt->onsetDFBufferSize - round(bt->beatPeriod/2));
        auto winsize = end-start+1;
        auto v = -2*bt->beatPeriod;
        auto sigma = detail::recip(bt->tightness);
        for (int i = 0;i < winsize;i++) {
            bt->w1[i] = detail::gaussian_warp(v, sigma, bt->beatPeriod);
    //        exp((-1*pow(bt->tightness*log(-v/bt->beatPeriod),2))/2);
            v = v+1;
        }
        bt->w1size = bt->beatPeriod;
    }

	// calculate new cumulative score value
	maxw = 0;
    max = 0;
    min = bt->cumulativeScore[start];
	int n = 0;
	for (int i=start;i <= end;i++) {
        wcumscore = bt->cumulativeScore[i]*bt->w1[n];
        if (wcumscore > maxw) {
            maxw = wcumscore;
        }
        if (bt->cumulativeScore[i] < min) {
            min = bt->cumulativeScore[i];
        }
        if (bt->cumulativeScore[i] > max) {
            max = bt->cumulativeScore[i];
        }
		n++;
	}


	// shift cumulative score back one
	for (int i = 0;i < (bt->onsetDFBufferSize-1);i++) {
		bt->cumulativeScore[i] = bt->cumulativeScore[i+1];
	}

	// add new value to cumulative score
	bt->cumulativeScore[bt->onsetDFBufferSize-1] = ((1-bt->alpha)*odfSample) + (bt->alpha*maxw);
	bt->latestCumulativeScoreValue = bt->cumulativeScore[bt->onsetDFBufferSize-1];
    bt->latestConfidence = 1.0 - (min / max);
}

static void predictBeat(struct btrack * bt){
	int windowSize = (int) bt->beatPeriod;
	float futureCumulativeScore[bt->onsetDFBufferSize + windowSize + 1];
    float w2[windowSize + 1];

	// copy cumscore to first part of fcumscore
	for (int i = 0;i < bt->onsetDFBufferSize;i++) {
		futureCumulativeScore[i] = bt->cumulativeScore[i];
	}

	// create future window
	for (int i = 0;i < windowSize;i++) {
		w2[i] = detail::gaussian_unorm(float(i+1), bt->beatPeriod/2.f, bt->beatPeriod/2.f);
	}
	// create window
    if(bt->beatPeriod != bt->w1size) {
        auto start = int(bt->onsetDFBufferSize - round(2*bt->beatPeriod));
        auto end = int(bt->onsetDFBufferSize - round(bt->beatPeriod/2));
        auto winsize = end-start+1;
        auto v = -2*bt->beatPeriod;
        auto sigma = detail::recip(bt->tightness);
        for (auto i = 0;i < winsize;) {
            bt->w1[i++] = detail::gaussian_warp(v++, sigma, bt->beatPeriod);
        }
        bt->w1size = bt->beatPeriod;
    }
	// calculate future cumulative score
	float max;
	int n;
	float wcumscore;
	for (int i = bt->onsetDFBufferSize;i < (bt->onsetDFBufferSize + windowSize); i++) {
		start = i - round(2*bt->beatPeriod);
		end = i - round(bt->beatPeriod/2);

		max = 0;
		n = 0;
        wcumscore = 0.f;
		for (int k=start;k <= end;k++) {
			wcumscore = std::max(wcumscore,futureCumulativeScore[k]*bt->w1[n]);
			n++;
		}

		futureCumulativeScore[i] = wcumscore;
	}


	// predict beat
	max = 0;
	n = 0;

	for (auto i = bt->onsetDFBufferSize;i < (bt->onsetDFBufferSize+windowSize);i++) {
		wcumscore = futureCumulativeScore[i]*w2[n];
		if (wcumscore > max) {
			max = wcumscore;
			bt->beatCounter = n;
		}
		n++;
	}
	// set next prediction time
	bt->m0 = bt->beatCounter+round(bt->beatPeriod/2);
}
