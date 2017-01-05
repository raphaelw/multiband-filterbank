//
//  AllpassCrossover.hpp
//
//  Created by Raphael Pistachio on 22.12.16.
//  Copyright (c) 2016 Raphael Pistachio. All rights reserved.
//

#pragma once

#include <vector>
#include <algorithm>

#include <cmath>

#include "math_constants.h"

namespace afx {
namespace emqf {

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
//template <typename T>
//class AllpassFilterStateless {
//public:
//    // coefficients
//    T a1;
//    T a2;
//    
//    // beta (squared magnitude of the purely imaginary pole
//    // of the startup half-band filter)
//    T beta;
//    
//    // switch between 1st and 2nd order allpass
//    bool isBiquad;
//    
//    //------------------------------
//    
//    // setBeta, setBiquad
//    // init() ?
//    
//    AllpassFilterStateless():a1(0),a2(0),beta(0),isBiquad(true) {}
//    
//    inline void init(bool isBiquad, T beta = T(0)) {
//        // reset coefficients
//        a1 = T(0);
//        a2 = T(0);
//        
//        this->beta = beta;
//        this->isBiquad = isBiquad;
//        // IDEA: tune here to half-band
//    }
//    
//    inline void __tuneCrossoverFrequency(T alpha_1, T alpha) {
//        if (isBiquad) {
//            T beta_i = (beta + alpha_1*alpha_1) / (beta*alpha_1*alpha_1 + T(1));
//            a1 = alpha * (T(1) + beta_i);
//            a2 = beta_i;
//        } else {
//            a1 = alpha_1;
//        }
//    }
//    
//    inline void _tuneCrossoverFrequency(T frequency, T& _alpha_1, T& _alpha) {
//        T f_3dB = (T(1)-frequency)/T(2); // mirror and scale
//        
//        T alpha;
//        
//        T tang = std::tan(M_PI*f_3dB);
//        T alpha_1 = (T(1)-tang) / (T(1)+tang);
//        if (isBiquad) {
//            alpha = (T(1)-(tang*tang)) / (T(1)+(tang*tang));
//        }
//        
//        __tuneCrossoverFrequency(alpha_1, alpha);
//        
//        _alpha_1 = alpha_1;
//        _alpha   = alpha;
//    }
//    
//    // Tune crossover frequency (0=DC, 1=Nyquist)
//    inline void tuneCrossoverFrequency(T frequency) {
//        T alpha_1;
//        T alpha;
//        _tuneCrossoverFrequency(frequency, alpha_1, alpha);
//    }
//
//    
//    inline void process(int numSamples, T* input, T* output, AllpassState<T>* state) {
//        T q_nm1 = state->q_nm1;
//        T p_nm1 = state->p_nm1;
//        
//        for (int i=0; i<numSamples; i++) {
//            T x = input[i];
//            
//            T y;
//            if (isBiquad) {
//                // for allpass case: b0 = a2; b1 = a1; b2 = a0 = 1
//                y = q_nm1 + a2*x;
//                q_nm1 = p_nm1 + a1*x - a1*y;
//                p_nm1 = x - a2*y;
//            } else {
//                // 1st order
//                // for allpass case: b0 = a1; b1 = a0 = 1;
//                y = q_nm1 + a1*x;
//                q_nm1 = x - a1*y;
//            }
//            
//            output[i] = y;
//        }
//        
//        state->q_nm1 = q_nm1;
//        state->p_nm1 = p_nm1;
//    }
//    
//    inline void processHB(T* input, T* output, int numSamples, AllpassState<T>* state) {
//        T q_nm1 = state->q_nm1;
//        T p_nm1 = state->p_nm1;
//        
//        for (int i=0; i<numSamples; i++) {
//            T x = input[i];
//            
//            T y;
//            if (isBiquad) {
//                // for allpass/HB case: b0 = a2 = beta; b1 = a1 = 0; b2 = a0 = 1
//                y = q_nm1 + beta*x;
//                q_nm1 = p_nm1;
//                p_nm1 = x - beta*y;
//            } else {
//                // pure Delay
//                y = q_nm1;
//                q_nm1 = x;
//            }
//
//            output[i] = y;
//        }
//        
//        state->q_nm1 = q_nm1;
//        state->p_nm1 = p_nm1;
//    }
//};
    

// multichannel (statefull) 2nd order allpass, transposed direct form II
// usage: init, tune, process
// IDEA: use custom allocator for channelStates to achieve an contiguous arena of filters
template <typename T>
class AllpassFilterMultichannel {
public:
    // state
    std::vector< AllpassState<T> > channelStates;
    
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
    
    AllpassFilterMultichannel(int numChannels)
        : channelStates(numChannels)
        , a1(0)
        , a2(0)
        , beta(0)
        , isBiquad(true)
    {}
    
    inline void resetStates() {
        std::fill(channelStates.begin(), channelStates.end(), AllpassState<T>());
    }
    
    // note that you should call also call resetStates(),
    //   especially if the filter order changed
    inline void setFilter(bool isBiquad, T beta = T(0)) {
        // reset coefficients
        a1 = T(0);
        a2 = T(0);
        
        this->beta = beta;
        this->isBiquad = isBiquad;
        // IDEA: tune here to half-band
    }
    
    inline void __tuneCrossoverFrequency(T alpha_1, T alpha) {
        if (isBiquad) {
            T beta_i = (beta + alpha_1*alpha_1) / (beta*alpha_1*alpha_1 + T(1));
            a1 = alpha * (T(1) + beta_i);
            a2 = beta_i;
        } else {
            a1 = alpha_1;
        }
    }
    
    inline void _tuneCrossoverFrequency(T frequency, T& _alpha_1, T& _alpha) {
        T f_3dB = (T(1)-frequency)/T(2); // mirror and scale
        
        T alpha;
        
        T tang = std::tan(M_PI*f_3dB);
        T alpha_1 = (T(1)-tang) / (T(1)+tang);
        if (isBiquad) {
            alpha = (T(1)-(tang*tang)) / (T(1)+(tang*tang));
        }
        
        __tuneCrossoverFrequency(alpha_1, alpha);
        
        _alpha_1 = alpha_1;
        _alpha   = alpha;
    }
    
    // Tune crossover frequency (0=DC, 1=Nyquist)
    inline void tuneCrossoverFrequency(T frequency) {
        T alpha_1;
        T alpha;
        _tuneCrossoverFrequency(frequency, alpha_1, alpha);
    }
    
    
    inline void process(int channel, int numSamples, T* input, T* output) {
        AllpassState<T> state = channelStates.at(channel);
        T q_nm1 = state.q_nm1;
        T p_nm1 = state.p_nm1;
        
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
        
        state.q_nm1 = q_nm1;
        state.p_nm1 = p_nm1;
        channelStates[channel] = state;
    }
    
    inline void processHB(int channel, int numSamples, T* input, T* output) {
        AllpassState<T> state = channelStates.at(channel);
        T q_nm1 = state.q_nm1;
        T p_nm1 = state.p_nm1;
        
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
        
        state.q_nm1 = q_nm1;
        state.p_nm1 = p_nm1;
        channelStates[channel] = state;
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

// TODO:
    // verify pole distribution (both cases, numPair odd&even) in debugging mode
    // ? verify tuning formula ?
    
    // implement rest of methods: minimum processPhase()
    
// crossover owns state, which is allocated outside
// no bounds checking is done by this class, make shure you allocated engough states & filters

// the xover needs:
// - filter instances: (numPair+1)=numFiltersPerChannel
// - states: numChannels*numFiltersPerChannel
template <typename T>
class Crossover {
    // member types -------------------
    typedef std::vector< AllpassFilterMultichannel<T> > FilterList;
    typedef typename FilterList::iterator FilterIterator;
    
    // members ------------------------
    FilterList filters;
    
    
    int numFiltersPerChannel; // number of ACTIVE filters
    //int numChannels; // filters[0].size()
    int offset; // offset to the second filter within the channel
public:
    Crossover(int numChannels, int numMaxFiltersPerChannel)
        : filters(numMaxFiltersPerChannel, AllpassFilterMultichannel<T>(numChannels) )
        , offset(0)
        , numFiltersPerChannel(0)
    {
        //reset();
    }
    
    
    inline void tuneCrossoverFrequency(T f) {
        T alpha_1;
        T alpha;
        for (int i = 0; i < numFiltersPerChannel; i++) {
            if (i==0) {
                filters[i]._tuneCrossoverFrequency(f, alpha_1, alpha);
            } else {
                filters[i].__tuneCrossoverFrequency(alpha_1, alpha);
            }
        }
        // no smoothing yet -> reset delay states too? would need numChannels
    }
    
    inline void setEMQFHalfbandFilter(T* betas, int numPairs) {
        numFiltersPerChannel = numPairs+1;
        
        bool firstSection;
        if (numPairs%2) {
            // section containing the real pole has lower order
            offset = (numPairs+1)/2;
            firstSection = true;
        } else {
            // section without real pole has lower order
            offset = numPairs/2;
            firstSection = false;
        }
        
        FilterIterator firstSectionFilter  = filters.begin()+offset-1;
        FilterIterator secondSectionFilter = filters.begin()+numFiltersPerChannel-1;
        
        //AllpassFilterStateless<T>* firstSectionFilter  = filters+offset-1;
        //AllpassFilterStateless<T>* secondSectionFilter = filters+numFiltersPerChannel-1;
        
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
                firstSectionFilter->setFilter(isBiquad, beta);
                firstSectionFilter->resetStates();
                
                firstSectionFilter--;
            } else {
                secondSectionFilter->setFilter(isBiquad, beta);
                secondSectionFilter->resetStates();
                
                secondSectionFilter--;
            }
            
            firstSection = (!firstSection); // switch section
        }
    }
    
    // first channel is zero
    // case numPairs==0: only works if one of the output buffers is also input
    inline void process(int channel, int numSamples
                        , T* input, T* outputLowpass, T* outputHighpass
                        , bool onlyLowpass = false) {
        T* buf1 = outputLowpass;
        T* buf2 = outputHighpass;
        if (input == outputLowpass) {
            // swap to prevent overwriting of input in 1st section
            buf1 = outputHighpass;
            buf2 = outputLowpass;
        }
        if (offset == 0) {
            // numPairs = 0; one-pole case
            // swap buffers
            T* tmp = buf1;
            buf1 = buf2;
            buf2 = tmp;
        }
        
        
        // 1st section
        bool firstPass = true;
        for (int i=0; i<offset; i++) {
            T* in;
            if (firstPass) {
                firstPass = false;
                in = input;
            } else {
                in = buf1;
            }
            
            filters[i].process(channel, numSamples, in, buf1);
        }
        // 2nd section
        firstPass = true;
        for (int i=offset; i<numFiltersPerChannel; i++) {
            T* in;
            if (firstPass) {
                firstPass = false;
                in = input;
            } else {
                in = buf2;
            }
            
            filters[i].process(channel, numSamples, in, buf2);
        }
        
        // criss-cross
        for (int i=0; i<numSamples; i++) {
            T low = (buf1[i] + buf2[i]) / T(2);
            
            if (!onlyLowpass) {
                T hi  = (buf1[i] - buf2[i]) / T(2);
                outputHighpass[i] = hi;
            }
            
            outputLowpass[i] = low;
        }
    }
    
    
};


    
} // namespace emqf
} // namespace afx