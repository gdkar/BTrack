//=======================================================================
/** @file BTrack.h
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

#ifndef __BTRACK_H
#define __BTRACK_H

#include <memory>
#include <utility>
#include <functional>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iterator>

#include "SlideBuffer.hpp"
#include "ACF.hpp"
#include "OnsetDetectionFunction.h"


//=======================================================================
/** The main beat tracking class and the interface to the BTrack
 * beat tracking algorithm. The algorithm can process either
 * audio frames or onset detection function samples and also
 * contains some static functions for calculating beat times in seconds
 */
struct btrack {
    btrack(int hop_size, int frame_size, float sample_rate);
    btrack() = default;
   ~btrack() = default;

    void configure();
    bool    m_configured{false};
    odf     m_odf;
    float   m_sample_rate;
    int     m_frame_size;
    int     m_hop_size;
    float   m_odf_rate;
    float   m_acf_rate;

    float m_sample_time;
    float m_odf_time;
    float m_acf_time;

    float m_acf_duration;
    int   m_odf_buf_size;
    int   m_acf_buf_size;
    int   m_acf_size;

    retrack::ACF m_acf{};

    std::unique_ptr<float[]> onsetDF;
           /**< to hold onset detection function */
    std::unique_ptr<float[]> cumulativeScore;               /**< to hold cumulative score */
    std::unique_ptr<float[]> w1;
    float                    w1period{};
    int                      w1start{};
    int                      w1end{};
    float                    w1threshold{};

    struct tempo_type {
        float   tempo{};
        float   period{};
        std::array<int, 2> lags{};
    };
    std::vector<tempo_type> m_tempos{};

    struct weight_type{
        int src_skip{};
        int dst_skip{};
        std::vector<float> weights{};
    };
    std::vector<weight_type> m_weights{};
    float                    m_weights_threshold{};

    using viterbi_type = std::pair<float, int>;

    retrack::SlideBuffer<float>        m_odf_buf{};
    retrack::SlideBuffer<viterbi_type> m_cum_buf{};
    std::vector<float>                 m_acf_buf{};
    std::vector<float>                 m_comb_buf{};

    std::vector<viterbi_type>          m_delta{};
    std::vector<viterbi_type>          m_delta_prev{};
    std::vector<viterbi_type>          m_delta_fixed{};

    float resampledOnsetDF[512];           /**< to hold resampled detection function */
    float acf[512];                        /**<  to hold autocorrelation function */
    float weightingVector[128];            /**<  to hold weighting vector */
    float combFilterBankOutput[128];       /**<  to hold comb filter output */
    float tempoObservationVector[41];      /**<  to hold tempo version of comb filter output */

    float delta[41];                       /**<  to hold final tempo candidate array */
    float prevDelta[41];                   /**<  previous delta */
    float prevDeltaFixed[41];              /**<  fixed tempo version of previous delta */

    float tempoTransitionMatrix[41][41];   /**<  tempo transition matrix */

	//=======================================================================
    // parameters

    float tightness;                       /**< the tightness of the weighting used to calculate cumulative score */
    float alpha;                           /**< the mix between the current detection function sample and the cumulative score's "momentum" */
    float beatPeriod;                      /**< the beat period, in detection function samples */
    float tempo;                           /**< the tempo in beats per minute */
    float estimatedTempo;                  /**< the current tempo estimation being used by the algorithm */
    float latestCumulativeScoreValue;      /**< holds the latest value of the cumulative score function */
    float latestODF;                       /**< holds the latest value of the onset detection function*/
    float latestConfidence;                /**< holds the latest confidence value, the ratio between max score & min score in the last beat */
    float tempoToLagFactor;                /**< factor for converting between lag and tempo */
    int m0;                                 /**< indicates when the next point to predict the next beat is */
    int beatCounter;                        /**< keeps track of when the next beat is - will be zero when the beat is due, and is set elsewhere in the algorithm to be positive once a beat prediction is made */
    int onsetDFBufferSize;                  /**< the onset detection function buffer size */
    int tempoFixed;                         /**< indicates whether the tempo should be fixed or not */
    int beatDueInFrame;                     /**< indicates whether a beat is due in the current frame */

};

/** Constructor taking both hopSize and frameSize
* @param hop_size the hop size in audio samples
* @param frame_size the frame size in audio samples
*/
void btrack_process_audio_frame(struct btrack * bt, const btrack_chunk_t * frame);
void btrack_process_fft_frame(struct btrack * bt, const btrack_chunk_t * fft_frame);
void btrack_process_odf_sample(struct btrack * bt, float odf_sample);

int btrack_beat_due_in_current_frame(const struct btrack * bt);
float btrack_get_bpm(const struct btrack * bt);
float btrack_get_latest_score(const struct btrack * bt);
float btrack_get_latest_odf(const struct btrack * bt);
float btrack_get_latest_confidence(const struct btrack * bt);
int btrack_get_frames_until_beat(const struct btrack * bt);
float btrack_get_time_until_beat(const struct btrack * bt);

void btrack_set_bpm(struct btrack * bt, float bpm);
void btrack_fix_bpm(struct btrack * bt, float bpm);
void btrack_nofix_bpm(struct btrack * bt);

void btrack_set_hop_size(struct btrack * bt, int hop_size);

#endif
