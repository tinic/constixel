/*
MIT License

Copyright (c) 2025 Tinic Uro

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// API tests run by ctest. Prints nothing on success; prints the failed checks and returns 1 otherwise.

#include "constixel.hpp"

#include <cstdio>

namespace {

int failures = 0;

void check(bool ok, const char *what, int line) {
    if (!ok) {
        printf("FAILED line %d: %s\n", line, what);
        failures++;
    }
}

#define CHECK(cond) check((cond), #cond, __LINE__)

// ----- image basics ---------------------------------------------------------------------

consteval bool consteval_image_basics() {
    constixel::image<constixel::format_1bit, 16, 8> image;
    image.plot(3, 2, 1);
    return image.data_ref()[2 * 2] == 0x10 && image.size() == 16;
}
static_assert(consteval_image_basics());

// formats without a runtime palette are pixel storage only
static_assert(sizeof(constixel::image<constixel::format_1bit, 192, 64>) == 1536);

void test_image_basics() {
    constixel::image<constixel::format_8bit, 16, 8> image;
    image.fill_rect(1, 1, 2, 2, constixel::color::WHITE);
    CHECK(image.data_ref()[1 * 16 + 1] == constixel::color::WHITE);
    CHECK(image.data_ref()[3 * 16 + 3] == 0);
}

}  // namespace

int main() {
    test_image_basics();
    return failures == 0 ? 0 : 1;
}
