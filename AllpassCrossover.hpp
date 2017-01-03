//
//  AllpassCrossover.hpp
//
//  Created by Raphael Pistachio on 22.12.16.
//  Copyright (c) 2016 Raphael Pistachio. All rights reserved.
//

namespace afx {
namespace emqf {

// State of delays of a 2nd order allpass, transposed direct form II
template <typename T>
struct AllpassState {
    T q_nm1;
    T p_nm1;
};


// stateless (does not own delays) 2nd order allpass, transposed direct form II
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
    
    // tuneCrossoverFrequency(T frequency);
    
    inline void process(T *input, T *output, int numSamples, AllpassState<T> *state) {
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
    
    inline void processHB(T *input, T *output, int numSamples, AllpassState<T> *state) {
        T q_nm1 = state->q_nm1;
        T p_nm1 = state->p_nm1;
        
        for (int i=0; i<numSamples; i++) {
            T x = input[i];
            
            // for allpass/HB case: b0 = a2 = beta; b1 = a1 = 0; b2 = a0 = 1
            T y = q_nm1 + beta*x;
            q_nm1 = p_nm1;
            p_nm1 = x - beta*y;

            output[i] = y;
        }
        
        state->q_nm1 = q_nm1;
        state->p_nm1 = p_nm1;
    }
};




    
} // namespace xover
} // namespace afx