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

#include "fonts/ibmplexsans_bold_18_aa.hpp"
#include "fonts/ibmplexsans_bold_18_mono.hpp"

#include <array>
#include <cstdio>
#include <string>
#include <string_view>

namespace {

int failures = 0;

void check(bool ok, int line) {
    if (!ok) {
        printf("FAILED: check at line %d\n", line);
        failures++;
    }
}


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
    check(image.data_ref()[1 * 16 + 1] == constixel::color::WHITE, __LINE__);
    check(image.data_ref()[3 * 16 + 3] == 0, __LINE__);
}

// ----- text takes std::string_view ------------------------------------------------------

using mono = constixel::ibmplexsans_bold_18_mono;
using aa = constixel::ibmplexsans_bold_18_aa;
using image8 = constixel::image<constixel::format_8bit, 160, 48>;
using image1 = constixel::image<constixel::format_1bit, 64, 24>;

template <typename I>
[[nodiscard]] constexpr bool same_pixels(I &a, I &b) {
    return a.data_ref() == b.data_ref();
}

template <typename I>
[[nodiscard]] constexpr bool blank(I &a) {
    for (auto v : a.data_ref()) {
        if (v != 0) {
            return false;
        }
    }
    return true;
}

// A view which is not NUL-terminated: reading one byte past it is not a constant expression.
consteval bool consteval_string_view() {
    constexpr std::array<char, 2> raw{'A', 'V'};
    // ends in the middle of the three byte sequence E2 9A A0
    constexpr std::array<char, 3> cut{'A', '\xE2', '\x9A'};
    image1 a;
    image1 b;
    const int32_t end_a = a.draw_string_mono<mono, true>(1, 1, "AV", 1);
    const int32_t end_b = b.draw_string_mono<mono, true>(1, 1, std::string_view{raw.data(), raw.size()}, 1);
    bool ok = end_a == end_b && same_pixels(a, b) && !blank(a);
    ok = ok && a.string_width<mono, true>("AV") == b.string_width<mono, true>(std::string_view{raw.data(), raw.size()});
    a.clear();
    b.clear();
    a.draw_string_mono<mono, true>(1, 1, "A\xE2\x9A", 1);
    b.draw_string_mono<mono, true>(1, 1, std::string_view{cut.data(), cut.size()}, 1);
    return ok && same_pixels(a, b) && !blank(a);
}
static_assert(consteval_string_view());

void test_string_view() {
    static constexpr char text[] = "Hello W\xC3\xB6rld";         // 12 bytes
    static constexpr char buffer[] = "xxxHello W\xC3\xB6rldyyy";  // the same text inside a larger buffer
    const std::string_view inner{buffer + 3, 12};
    image8 a;
    image8 b;

    // literal, view into the middle of a buffer, std::string
    const int32_t end = a.draw_string_mono<mono, true>(2, 4, text, 1);
    check(!blank(a), __LINE__);
    check(b.draw_string_mono<mono, true>(2, 4, inner, 1) == end, __LINE__);
    check(same_pixels(a, b), __LINE__);
    b.clear();
    check(b.draw_string_mono<mono, true>(2, 4, std::string{text}, 1) == end, __LINE__);
    check(same_pixels(a, b), __LINE__);

    // the kerning lookahead stops at the end of the view: "AV" out of "AVA"
    a.clear();
    b.clear();
    a.draw_string_mono<mono, true>(2, 4, "AV", 1);
    b.draw_string_mono<mono, true>(2, 4, std::string_view{"AVA", 2}, 1);
    check(same_pixels(a, b), __LINE__);
    check(a.string_width<mono, true>("AV") == a.string_width<mono, true>(std::string_view{"AVA", 2}), __LINE__);

    // a view which ends in the middle of a UTF-8 sequence
    static constexpr char warning[] = "A\xE2\x9A\xA0" "B";
    a.clear();
    b.clear();
    a.draw_string_mono<mono>(2, 4, "A\xE2\x9A", 1);
    b.draw_string_mono<mono>(2, 4, std::string_view{warning, 3}, 1);
    check(same_pixels(a, b), __LINE__);
    check(a.string_width<mono>("A\xE2\x9A") == a.string_width<mono>(std::string_view{warning, 3}), __LINE__);

    // width and character counts
    size_t actual_a = 0;
    size_t actual_b = 0;
    check(a.string_width<mono, true>(text, 5, &actual_a) == a.string_width<mono, true>(inner, 5, &actual_b), __LINE__);
    check(actual_a == actual_b && actual_a != 0, __LINE__);
    check(a.string_width<aa>(text) == a.string_width<aa>(inner), __LINE__);

    // antialiased, centered and rotated
    a.clear();
    b.clear();
    check(a.draw_string_aa<aa, true>(2, 4, text, 3) == b.draw_string_aa<aa, true>(2, 4, inner, 3), __LINE__);
    check(same_pixels(a, b) && !blank(a), __LINE__);
    a.clear();
    b.clear();
    a.draw_string_centered_mono<mono, true>(80, 4, text, 1);
    b.draw_string_centered_mono<mono, true>(80, 4, inner, 1);
    check(same_pixels(a, b) && !blank(a), __LINE__);
    a.clear();
    b.clear();
    a.draw_string_centered_aa<aa, false, constixel::DEGREE_180>(80, 30, text, 1);
    b.draw_string_centered_aa<aa, false, constixel::DEGREE_180>(80, 30, inner, 1);
    check(same_pixels(a, b) && !blank(a), __LINE__);

    // fluent shapes
    a.clear();
    b.clear();
    a.text_mono<mono, true>(2, 4, text).color(1);
    b.text_mono<mono, true>(2, 4, inner).color(1);
    check(same_pixels(a, b) && !blank(a), __LINE__);
    a.clear();
    b.clear();
    check(a.text_aa<aa>(2, 4, text).draw(2) == b.text_aa<aa>(2, 4, inner).draw(2), __LINE__);
    a.text_centered_mono<mono>(80, 20, text).color(1);
    b.text_centered_mono<mono>(80, 20, inner).color(1);
    a.text_centered_aa<aa>(80, 4, text).color(4);
    b.text_centered_aa<aa>(80, 4, inner).color(4);
    check(same_pixels(a, b) && !blank(a), __LINE__);

    // the draw_string struct and a null pointer keep working
    a.clear();
    b.clear();
    a.draw_string_mono<mono>({.x = 2, .y = 4, .str = "Hello", .col = 1});
    b.draw_string_mono<mono>(2, 4, std::string_view{"Hello"}, 1);
    check(same_pixels(a, b) && !blank(a), __LINE__);
    a.clear();
    const char *null_string = nullptr;
    check(a.draw_string_mono<mono>(2, 4, null_string, 1) == 0, __LINE__);
    check(a.draw_string_aa<aa>(2, 4, null_string, 1) == 0, __LINE__);
    check(blank(a), __LINE__);
}

}  // namespace

int main() {
    test_image_basics();
    test_string_view();
    return failures == 0 ? 0 : 1;
}
