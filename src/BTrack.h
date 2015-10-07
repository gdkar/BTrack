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

_Pragma("once")
#ifndef _ISOC11_SOURCE
#define _ISOC11_SOURCE
#endif
#include "OnsetDetectionFunction.h"

//=======================================================================
/** The main beat tracking class and the interface to the BTrack
 * beat tracking algorithm. The algorithm can process either
 * audio frames or onset detection function samples and also
 * contains some static functions for calculating beat times in seconds
 */
struct btrack {
    struct odf odf;
    int frameSize;
    float * onsetDF;                       /**< to hold onset detection function */
    float * cumulativeScore;               /**< to hold cumulative score */
    float * w1;
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
    float tempoToLagFactor;                /**< factor for converting between lag and tempo */
    int m0;                                 /**< indicates when the next point to predict the next beat is */
    int beatCounter;                        /**< keeps track of when the next beat is - will be zero when the beat is due, and is set elsewhere in the algorithm to be positive once a beat prediction is made */
    int hopSize;                            /**< the hop size being used by the algorithm */
    int onsetDFBufferSize;                  /**< the onset detection function buffer size */
    int tempoFixed;                         /**< indicates whether the tempo should be fixed or not */
    int beatDueInFrame;                     /**< indicates whether a beat is due in the current frame */
};
/** Constructor taking both hopSize and frameSize
* @param hop_size the hop size in audio samples
* @param frame_size the frame size in audio samples
*/
int btrack_init(struct btrack * bt, int hop_size, int frame_size);
void btrack_del(struct btrack * bt);
void btrack_process_audio_frame(struct btrack * bt, float * frame);
void btrack_process_odf_sample(struct btrack * bt, float odf_sample);
int btrack_beat_due_in_current_frame(struct btrack * bt);
float btrack_get_bpm(struct btrack * bt);
float btrack_get_latest_score(struct btrack * bt);
float btrack_get_latest_odf(struct btrack * bt);
void btrack_set_bpm(struct btrack * bt, float bpm);
void btrack_fix_bpm(struct btrack * bt, float bpm);
void btrack_nofix_bpm(struct btrack * bt);
void btrack_set_hop_size(struct btrack * bt, int hop_size);
