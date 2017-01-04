//
//  AllpassCrossover.hpp
//
//  Created by Raphael Pistachio on 22.12.16.
//  Copyright (c) 2016 Raphael Pistachio. All rights reserved.
//

#pragma once

#include <cmath>

#include "math_constants.h"

namespace afx {
namespace xover {

// State of delays of a 2nd order allpass, transposed direct form II
template <typename T>
struct AllpassState {
    T q_nm1;
    T p_nm1;
    
    //------------------------------
    AllpassState():q_nm1(0), p_nm1(0) {}
};


// stateless (does not own delays) 2nd order allpass, transposed direct form II
// usage: init, tune, process
template <typename T>
class AllpassFilterStateless {
public:
    // coefficients
    T a1;
    T a2;
    
    // beta (squared magnitude of the purely imaginary pole
    // of the startup half-band filter)
    T beta;
    
    // switch between 1st and 2nd order allpass
    bool isBiquad;
    
    //------------------------------
    
    // setBeta, setBiquad
    // init() ?
    
    AllpassFilterStateless():a1(0),a2(0),beta(0),isBiquad(true) {}
    
    inline void init(bool isBiquad, T beta = T(0)) {
        // reset coefficients
        a1 = T(0);
        a2 = T(0);
        
        this->beta = beta;
        this->isBiquad = isBiquad;
        // IDEA: tune here to half-band
    }
    
    inline void tuneCrossoverFrequency(T alpha_1, T alpha) {
        if (isBiquad) {
            T beta_i = (beta + alpha_1*alpha_1) / (beta*alpha_1*alpha_1 + T(1));
            a1 = alpha * (T(1) + beta_i);
            a2 = beta_i;
        } else {
            a1 = alpha_1;
        }
    }
    
    // Tune crossover frequency (0=DC, 1=Nyquist)
    inline void tuneCrossoverFrequency(T frequency) {
        T f_3dB = (T(1)-frequency)/T(2); // mirror and scale
        
        T alpha;
        
        T tang = std::tan(M_PI*f_3dB);
        T alpha_1 = (T(1)-tang) / (T(1)+tang);
        if (isBiquad) {
            alpha = (T(1)-(tang*tang)) / (T(1)+(tang*tang));
        }
        
        tuneCrossoverFrequency(alpha_1, alpha);
    }
    
    
    inline void process(T* input, T* output, int numSamples, AllpassState<T>* state) {
        T q_nm1 = state->q_nm1;
        T p_nm1 = state->p_nm1;
        
        for (int i=0; i<numSamples; i++) {
            T x = input[i];
            
            T y;
            if (isBiquad) {
                // for allpass case: b0 = a2; b1 = a1; b2 = a0 = 1
                y = q_nm1 + a2*x;
                q_nm1 = p_nm1 + a1*x - a1*y;
                p_nm1 = x - a2*y;
            } else {
                // 1st order
                // for allpass case: b0 = a1; b1 = a0 = 1;
                y = q_nm1 + a1*x;
                q_nm1 = x - a1*y;
            }
            
            output[i] = y;
        }
        
        state->q_nm1 = q_nm1;
        state->p_nm1 = p_nm1;
    }
    
    inline void processHB(T* input, T* output, int numSamples, AllpassState<T>* state) {
        T q_nm1 = state->q_nm1;
        T p_nm1 = state->p_nm1;
        
        for (int i=0; i<numSamples; i++) {
            T x = input[i];
            
            T y;
            if (isBiquad) {
                // for allpass/HB case: b0 = a2 = beta; b1 = a1 = 0; b2 = a0 = 1
                y = q_nm1 + beta*x;
                q_nm1 = p_nm1;
                p_nm1 = x - beta*y;
            } else {
                // pure Delay
                y = q_nm1;
                q_nm1 = x;
            }

            output[i] = y;
        }
        
        state->q_nm1 = q_nm1;
        state->p_nm1 = p_nm1;
    }
};

// crossover owns state, which is allocated outside
// setFilters(AllpassFilterStateless*, int numFilterPerChannel)
// setChannels(AllpassState*, int numChannels)

// setEMQFHalfbandFilter(T* betas, int numPairs) // distribute poles and reset states
// tuneCrossoverFrequency(T f) // iter and tune

// process(int channel, int numSamples, in, outLow, outHigh)
// processHB(int channel, int numSamples, in, outLow, outHigh)
// processPhase(int channel, int numSamples, in, out)
    
// no bounds checking is done by this class, make shure you allocated engough states & filters
template <typename T>
class Crossover {
    AllpassFilterStateless<T>* filters;
    AllpassState<T>* states;
    
    int numFiltersPerChannel;
    int numChannels;
    int offset; // offset to the second filter within the channel
public:
    Crossover() {
        reset();
    }
    
    inline void init(AllpassFilterStateless<T>* filters, AllpassState<T>* states, int numChannels = 1) {
        reset();
        this->filters = filters;
        this->states = states;
    }
    
    inline void reset() {
        numFiltersPerChannel = numChannels = offset = 0;
    }
    
    inline void setEMQFHalfbandFilter(T* betas, int numPairs, int numChannels) {
        numFiltersPerChannel = numPairs+1;
        
        bool firstSection;
        if (numPairs%2) {
            // section containing the real pole has lower order
            offset = (numPairs+1)/2;
            firstSection = true;
        } else {
            // section without real pole has lower order
            offset = (numPairs/2);
            firstSection = false;
        }
        
        AllpassFilterStateless<T>* firstSectionFilter  = filters+offset-1;
        AllpassFilterStateless<T>* secondSectionFilter = filters+numFiltersPerChannel-1;
        
        for (int i = 0; i < numFiltersPerChannel; i++) {
            T beta = 0;
            bool isBiquad = true;
            
            if (i == 0) { // real pole
                isBiquad = false;
            } else {
                isBiquad = true;
                beta = betas[i-1];
            }
            
            if (firstSection) {
                firstSectionFilter->init(isBiquad, beta);
                firstSectionFilter--;
            } else {
                secondSectionFilter->init(isBiquad, beta);
                secondSectionFilter--;
            }
            
            firstSection = (!firstSection); // switch section
        }
        
        for (int i = 0; i < (numChannels*numFiltersPerChannel); i++) {
            states[i] = AllpassState<T>();
        }
    }
    
    
};


    
} // namespace xover
} // namespace afx