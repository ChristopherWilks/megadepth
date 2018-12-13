//
// Created by Ben Langmead on 2018-12-12.
//

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main() {
    std::cout << CharString("Hello SeqAn!") << std::endl;
    return 0;
}
