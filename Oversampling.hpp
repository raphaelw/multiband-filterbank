//
//  Oversampling.hpp
//
//  Created by Raphael Pistachio
//  Copyright (c) 2017 Raphael Pistachio. All rights reserved.
//

#pragma once

#include "AllpassCrossover.hpp"


namespace afx {
namespace emqf {

template <typename T>
class Resampler {
    // members -------------------------
    Crossover<T> xover;
    int numStages;

public:
    Resampler():xover,numStages(0){}
};


} // namespace emqf
} // namespace afx