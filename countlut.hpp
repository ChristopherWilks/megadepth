/** 
Copyright (C) 2014 Milo Yip

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
**/
#include "countdecimaldigit.h"
#include "digitslut.h"


// Additional count number of digit pass
// Use lookup table of two gDigitsLut

uint32_t u32toa_countlut(uint32_t value, char* buffer, char delim) {
    unsigned digit = CountDecimalDigit32(value);
    buffer += digit;
    *buffer = delim;
    //char* buffer_end = buffer;

    while (value >= 100) {
        const unsigned i = (value % 100) << 1;
        value /= 100;
        *--buffer = gDigitsLut[i + 1];
        *--buffer = gDigitsLut[i];
    }

    if (value < 10) {
        *--buffer = char(value) + '0';
    }
    else {
        const unsigned i = value << 1;
        *--buffer = gDigitsLut[i + 1];
        *--buffer = gDigitsLut[i];
    }
    //return buffer_end;
    return digit;
}

void i32toa_countlut(int32_t value, char* buffer) {
    uint32_t u = static_cast<uint32_t>(value);
    if (value < 0) {
        *buffer++ = '-';
        u = ~u + 1;
    }
    u32toa_countlut(u, buffer, '\0');
}

void u64toa_countlut(uint64_t value, char* buffer) {
    unsigned digit = CountDecimalDigit64(value);
    buffer += digit;
    *buffer = '\0';

    while (value >= 100000000) {
        const uint32_t a = static_cast<uint32_t>(value % 100000000);
        value /= 100000000;

        const uint32_t b = a / 10000;
        const uint32_t c = a % 10000;

        const uint32_t b1 = (b / 100) << 1;
        const uint32_t b2 = (b % 100) << 1;
        const uint32_t c1 = (c / 100) << 1;
        const uint32_t c2 = (c % 100) << 1;

        buffer -= 8;

        buffer[0] = gDigitsLut[b1];
        buffer[1] = gDigitsLut[b1 + 1];
        buffer[2] = gDigitsLut[b2];
        buffer[3] = gDigitsLut[b2 + 1];
        buffer[4] = gDigitsLut[c1];
        buffer[5] = gDigitsLut[c1 + 1];
        buffer[6] = gDigitsLut[c2];
        buffer[7] = gDigitsLut[c2 + 1];
    }

    uint32_t value32 = static_cast<uint32_t>(value);
    while (value32 >= 100) {
        const unsigned i = static_cast<unsigned>(value32 % 100) << 1;
        value32 /= 100;
        *--buffer = gDigitsLut[i + 1];
        *--buffer = gDigitsLut[i];
    }

    if (value32 < 10) {
        *--buffer = char(value32) + '0';
    }
    else {
        const unsigned i = static_cast<unsigned>(value32) << 1;
        *--buffer = gDigitsLut[i + 1];
        *--buffer = gDigitsLut[i];
    }
}

void i64toa_countlut(int64_t value, char* buffer) {
    uint64_t u = static_cast<uint64_t>(value);
    if (value < 0) {
        *buffer++ = '-';
        u = ~u + 1;
    }
    u64toa_countlut(u, buffer);
}
