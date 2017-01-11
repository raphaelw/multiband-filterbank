//
//  AllpassFilterbank.hpp
//
//  Created by Raphael Pistachio
//  Copyright (c) 2017 Raphael Pistachio. All rights reserved.
//

#pragma once

#include <vector>

#include "AllpassCrossover.hpp"

namespace afx {
namespace emqf {

// Filterbank
    // buffers? -> maxBlockSize -> outside
    // numChannels -> numStages !
    // set EMQF startup-filter: betas or attenuation?
    
    // tune!
    
    // process
        // numSamples
        // inputBuffer
        // array of outputBuffers? T**
    
template <typename T>
class Filterbank {
public:
    // member types -------------------
    typedef Crossover<T> TCrossover;
    typedef std::vector<TCrossover> TCrossoverList;
    typedef typename TCrossoverList::iterator TCrossoverIterator;
    
    // members ------------------------
    TCrossoverList crossovers;
    
    // methods ------------------------
    Filterbank(int numStages, int numChannels, int maxNumPairs) {
        crossovers.reserve(numStages);
        for (int i = 0; i < numStages; i++) {
            crossovers.push_back(TCrossover( (i+1)*numChannels , maxNumPairs ));
        }
    }
    
    inline void setEMQFHalfbandFilter(T* betas, int numPairs) {
        for (int i = 0; i < crossovers.size(); i++) {
            crossovers[i].setEMQFHalfbandFilter(betas, numPairs);
        }
    }
    
    // frequency from 0 to 1, where 1 is nyquist
    inline void tuneCrossoverFrequency(int stage, T frequency) {
        crossovers[stage].tuneCrossoverFrequency(frequency);
    }
    
    // the last stage's highpass will be split up
    // first channel is zero
    inline void process(int channel, int numSamples, T* input, T** outputs) {
        int numStages = crossovers.size();
        T* bufferToSplit = input;

        for (int i = 0 ; i < numStages; i++) {
            int chOffset = (i+1)*channel;
            T* outLP = outputs[i];
            T* outHP = outputs[i+1];
            
            crossovers[i].process(chOffset+0, numSamples, bufferToSplit, outLP, outHP);
            
            bufferToSplit = outHP;
            
            for (int k = 0; k < i; k++) {
                crossovers[i].processPhase(chOffset+(k+1), numSamples, outputs[k], outputs[k]);
            }
        }
    }
    
};
    
    
    
} // namespace emqf
} // namespace afx