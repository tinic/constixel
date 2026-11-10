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
#include <vector>

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
    constexpr std::array<char, 2> raw{{'A', 'V'}};
    // ends in the middle of the three byte sequence E2 9A A0
    constexpr std::array<char, 3> cut{{'A', '\xE2', '\x9A'}};
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

// ----- draw_string: one entry point with alignment --------------------------------------

consteval bool consteval_draw_string() {
    image1 a;
    image1 b;
    // a mono font on a 1-bit image, right aligned
    const int32_t end_a = a.draw_string_mono<mono, true>(60 - a.string_width<mono, true>("AV"), 1, "AV", 1);
    const int32_t end_b = b.draw_string<mono, true>(60, 1, "AV", 1, constixel::ALIGN_RIGHT);
    return end_a == end_b && same_pixels(a, b) && !blank(a);
}
static_assert(consteval_draw_string());

void test_draw_string() {
    static constexpr char text[] = "Hello W\xC3\xB6rld";
    image8 a;
    image8 b;

    // dispatch on the kind of font
    check(a.draw_string_mono<mono, true>(2, 4, text, 1) == b.draw_string<mono, true>(2, 4, text, 1), __LINE__);
    check(same_pixels(a, b) && !blank(a), __LINE__);
    a.clear();
    b.clear();
    check(a.draw_string_aa<aa, true>(2, 4, text, 3) == b.draw_string<aa, true>(2, 4, std::string_view{text}, 3),
          __LINE__);
    check(same_pixels(a, b) && !blank(a), __LINE__);

    // right and center are the manual string_width arithmetic
    a.clear();
    b.clear();
    const int32_t width = a.string_width<mono, true>(text);
    a.draw_string_mono<mono, true>(150 - width, 4, text, 1);
    b.draw_string<mono, true>(150, 4, text, 1, constixel::ALIGN_RIGHT);
    check(same_pixels(a, b) && !blank(a), __LINE__);
    a.clear();
    b.clear();
    a.draw_string_mono<mono, true>(80 - width / 2, 4, text, 1);
    b.draw_string<mono, true>(80, 4, text, 1, constixel::ALIGN_CENTER);
    check(same_pixels(a, b) && !blank(a), __LINE__);
    a.clear();
    a.draw_string_centered_mono<mono, true>(80, 4, text, 1);
    check(same_pixels(a, b), __LINE__);

    // antialiased, right aligned, with a character limit
    a.clear();
    b.clear();
    const int32_t width5 = a.string_width<aa>(text, 5);
    a.draw_string_aa<aa>(150 - width5, 4, text, 2, 5);
    b.draw_string<aa>(150, 4, text, 2, constixel::ALIGN_RIGHT, 5);
    check(same_pixels(a, b) && !blank(a), __LINE__);

    // rotated: the offset follows the direction of the text
    a.clear();
    b.clear();
    a.draw_string_mono<mono, false, constixel::DEGREE_180>(10 + a.string_width<mono>("Hi"), 40, "Hi", 1);
    b.draw_string<mono, false, constixel::DEGREE_180>(10, 40, "Hi", 1, constixel::ALIGN_RIGHT);
    check(same_pixels(a, b) && !blank(a), __LINE__);

    // the param struct carries the alignment, also through the _mono and _aa names
    a.clear();
    b.clear();
    a.draw_string<mono, true>(150, 4, text, 1, constixel::ALIGN_RIGHT);
    b.draw_string<mono, true>({.x = 150, .y = 4, .str = text, .col = 1, .align = constixel::ALIGN_RIGHT});
    check(same_pixels(a, b) && !blank(a), __LINE__);
    b.clear();
    b.draw_string_mono<mono, true>({.x = 150, .y = 4, .str = text, .col = 1, .align = constixel::ALIGN_RIGHT});
    check(same_pixels(a, b), __LINE__);
    a.clear();
    b.clear();
    a.draw_string_aa<aa>(2, 4, text, 1);
    b.draw_string_aa<aa>({.x = 2, .y = 4, .str = text, .col = 1});
    check(same_pixels(a, b) && !blank(a), __LINE__);

    // dynamic_image and a bare surface have the same entry point
    constixel::dynamic_image<constixel::pixel_8bit<false>> d(160, 48);
    a.clear();
    a.draw_string<aa, true>(80, 4, text, 1, constixel::ALIGN_CENTER);
    d.draw_string<aa, true>(80, 4, text, 1, constixel::ALIGN_CENTER);
    b.clear();
    b.copy<160 * 48>(d.data_ref().data());
    check(same_pixels(a, b) && !blank(a), __LINE__);
    b.clear();
    b.view().draw_string<aa, true>(80, 4, text, 1, constixel::ALIGN_CENTER);
    check(same_pixels(a, b), __LINE__);
}

// ----- convert: rectangles, chunks, 1-bit display layouts -------------------------------

template <constixel::device_format FORMAT, typename I>
[[nodiscard]] std::vector<uint8_t> convert_all(const I &image) {
    std::vector<uint8_t> out;
    image.template convert<FORMAT>([&out](uint8_t byte) {
        out.push_back(byte);
    });
    return out;
}

template <constixel::device_format FORMAT, typename I>
[[nodiscard]] std::vector<uint8_t> convert_rect(const I &image, const constixel::rect<int32_t> &r) {
    std::vector<uint8_t> out;
    image.template convert<FORMAT>(
        [&out](uint8_t byte) {
            out.push_back(byte);
        },
        r);
    return out;
}

template <constixel::device_format FORMAT, typename I>
[[nodiscard]] std::vector<uint8_t> convert_chunks(const I &image, const constixel::rect<int32_t> &r, size_t chunk) {
    std::vector<uint8_t> out;
    std::vector<uint8_t> buffer(chunk);
    constixel::convert_cursor cursor;
    size_t calls = 0;
    while (size_t n = image.template convert_chunk<FORMAT>(buffer, cursor, r)) {
        out.insert(out.end(), buffer.begin(), buffer.begin() + static_cast<std::ptrdiff_t>(n));
        if (++calls > 100000) {
            break;
        }
    }
    // finished: further calls write nothing
    check(image.template convert_chunk<FORMAT>(buffer, cursor, r) == 0, __LINE__);
    return out;
}

// rows y .. y + h - 1, bytes x_bytes .. x_bytes + w_bytes - 1 of each line of a full stream
[[nodiscard]] std::vector<uint8_t> cut(const std::vector<uint8_t> &full, size_t line_bytes, size_t x_bytes, size_t y,
                                       size_t w_bytes, size_t h) {
    std::vector<uint8_t> out;
    for (size_t row = y; row < y + h; row++) {
        for (size_t i = 0; i < w_bytes; i++) {
            out.push_back(full[row * line_bytes + x_bytes + i]);
        }
    }
    return out;
}

template <typename I>
void draw_pattern(I &image) {
    for (int32_t y = 0; y < image.height(); y++) {
        for (int32_t x = 0; x < image.width(); x++) {
            image.plot(x, y, static_cast<uint8_t>((x * 7 + y * 13 + (x * y) / 3) & 0xFF));
        }
    }
}

template <constixel::device_format FORMAT, typename I>
void check_chunks(const I &image, const constixel::rect<int32_t> &r) {
    const std::vector<uint8_t> whole = convert_rect<FORMAT>(image, r);
    check(whole.size() == image.template convert_size<FORMAT>(r), __LINE__);
    for (size_t chunk : {size_t{1}, size_t{3}, size_t{7}, size_t{64}, size_t{1000}}) {
        check(convert_chunks<FORMAT>(image, r, chunk) == whole, __LINE__);
    }
}

// SSD1306 page layout and a chunked conversion, as constant expressions
consteval bool consteval_convert() {
    constixel::image<constixel::format_1bit, 16, 10> image;
    image.plot(0, 0, 1);
    image.plot(0, 7, 1);
    image.plot(1, 1, 1);
    image.plot(15, 8, 1);
    image.plot(3, 9, 1);
    std::array<uint8_t, 32> pages{};
    size_t count = 0;
    image.convert<constixel::MONO_1BIT_VERTICAL_PAGES>([&pages, &count](uint8_t byte) {
        pages[count++] = byte;
    });
    bool ok = count == 32 && image.convert_size<constixel::MONO_1BIT_VERTICAL_PAGES>() == 32;
    // page 0: column 0 has rows 0 and 7, column 1 has row 1. page 1 holds rows 8 and 9 only.
    ok = ok && pages[0] == 0x81 && pages[1] == 0x02 && pages[2] == 0x00 && pages[15] == 0x00;
    ok = ok && pages[16 + 15] == 0x01 && pages[16 + 3] == 0x02 && pages[16] == 0x00;

    std::array<uint8_t, 32> chunked{};
    std::array<uint8_t, 5> buffer{};
    constixel::convert_cursor cursor;
    size_t total = 0;
    while (size_t n = image.convert_chunk<constixel::MONO_1BIT_VERTICAL_PAGES>(buffer, cursor)) {
        for (size_t i = 0; i < n; i++) {
            chunked[total++] = buffer[i];
        }
    }
    return ok && total == 32 && chunked == pages;
}
static_assert(consteval_convert());

void test_convert() {
    // ---- 8-bit, RGB565: hand-computed colors, 0xRRRRRGGG 0xGGGBBBBB
    constixel::image<constixel::format_8bit, 32, 16> image;
    image.plot(0, 0, constixel::color::RED);
    image.plot(1, 0, constixel::color::LIME);
    image.plot(2, 0, constixel::color::BLUE);
    image.plot(3, 0, constixel::color::WHITE);
    const std::vector<uint8_t> first = convert_all<constixel::RGB565_8BIT_SERIAL>(image);
    check(first.size() == 32 * 16 * 2, __LINE__);
    check(first[0] == 0xF8 && first[1] == 0x00, __LINE__);
    check(first[2] == 0x07 && first[3] == 0xE0, __LINE__);
    check(first[4] == 0x00 && first[5] == 0x1F, __LINE__);
    check(first[6] == 0xFF && first[7] == 0xFF, __LINE__);
    check(first[8] == 0x00 && first[9] == 0x00, __LINE__);
    const std::vector<uint8_t> rgb666 = convert_all<constixel::RGB666_8BIT_SERIAL_1>(image);
    check(rgb666[0] == 0x3F && rgb666[1] == 0x00 && rgb666[2] == 0x00 && rgb666[5] == 0x00 && rgb666[8] == 0x3F,
          __LINE__);
    check(convert_all<constixel::RGB666_8BIT_SERIAL_2>(image)[0] == 0xFC, __LINE__);

    // ---- a rectangle is the corresponding bytes of the full stream
    draw_pattern(image);
    {
        const std::vector<uint8_t> full = convert_all<constixel::RGB565_8BIT_SERIAL>(image);
        check(convert_rect<constixel::RGB565_8BIT_SERIAL>(image, {5, 3, 10, 6}) == cut(full, 64, 10, 3, 20, 6),
              __LINE__);
        // clipped to the surface, and empty outside of it
        check(convert_rect<constixel::RGB565_8BIT_SERIAL>(image, {28, 12, 10, 10}) == cut(full, 64, 56, 12, 8, 4),
              __LINE__);
        check(convert_rect<constixel::RGB565_8BIT_SERIAL>(image, {-4, -2, 8, 4}) == cut(full, 64, 0, 0, 8, 2),
              __LINE__);
        check(convert_rect<constixel::RGB565_8BIT_SERIAL>(image, {40, 0, 4, 4}).empty(), __LINE__);
        check(image.convert_size<constixel::RGB565_8BIT_SERIAL>({40, 0, 4, 4}) == 0, __LINE__);
        check(image.convert_size<constixel::RGB565_8BIT_SERIAL>() == full.size(), __LINE__);
        const std::vector<uint8_t> full666 = convert_all<constixel::RGB666_8BIT_SERIAL_2>(image);
        check(convert_rect<constixel::RGB666_8BIT_SERIAL_2>(image, {5, 3, 10, 6}) == cut(full666, 96, 15, 3, 30, 6),
              __LINE__);
        const std::vector<uint8_t> raw = convert_all<constixel::STRAIGHT_THROUGH>(image);
        check(raw.size() == image.size() && raw[5 * 32 + 7] == image.data_ref()[5 * 32 + 7], __LINE__);
        check(convert_rect<constixel::STRAIGHT_THROUGH>(image, {5, 3, 10, 6}) == cut(raw, 32, 5, 3, 10, 6), __LINE__);
    }

    // ---- chunks concatenate to the one-shot stream
    check_chunks<constixel::RGB565_8BIT_SERIAL>(image, {0, 0, 32, 16});
    check_chunks<constixel::RGB565_8BIT_SERIAL>(image, {5, 3, 11, 7});
    check_chunks<constixel::RGB666_8BIT_SERIAL_1>(image, {5, 3, 11, 7});
    check_chunks<constixel::STRAIGHT_THROUGH>(image, {1, 1, 30, 14});
    {
        // without a rectangle
        std::vector<uint8_t> out;
        std::array<uint8_t, 13> buffer{};
        constixel::convert_cursor cursor;
        while (size_t n = image.convert_chunk<constixel::RGB565_8BIT_SERIAL>(buffer, cursor)) {
            out.insert(out.end(), buffer.begin(), buffer.begin() + static_cast<std::ptrdiff_t>(n));
        }
        check(out == convert_all<constixel::RGB565_8BIT_SERIAL>(image), __LINE__);
    }

    // ---- 24-bit and 4-bit sources give the colors of their pixels
    {
        constixel::image<constixel::format_24bit, 8, 4> rgb;
        rgb.plot(2, 1, constixel::color::RED);
        const std::vector<uint8_t> out = convert_rect<constixel::RGB565_8BIT_SERIAL>(rgb, {2, 1, 1, 1});
        check(out.size() == 2 && out[0] == 0xF8 && out[1] == 0x00, __LINE__);
        constixel::image<constixel::format_4bit, 9, 4> nibbles;
        draw_pattern(nibbles);
        nibbles.plot(3, 2, constixel::color::BLUE);
        check(convert_rect<constixel::RGB565_8BIT_SERIAL>(nibbles, {3, 2, 1, 1}) == std::vector<uint8_t>{0x00, 0x1F},
              __LINE__);
        // packed formats: the bytes which hold the pixels of the rectangle
        const std::vector<uint8_t> raw = convert_all<constixel::STRAIGHT_THROUGH>(nibbles);
        check(convert_rect<constixel::STRAIGHT_THROUGH>(nibbles, {3, 1, 4, 2}) == cut(raw, 5, 1, 1, 3, 2), __LINE__);
        check_chunks<constixel::RGB666_8BIT_SERIAL_1>(nibbles, {1, 0, 7, 4});
    }

    // ---- 1-bit display layouts
    {
        constixel::image<constixel::format_1bit, 20, 19> mono_image;
        mono_image.plot(0, 0, 1);
        mono_image.plot(0, 7, 1);
        mono_image.plot(1, 1, 1);
        mono_image.plot(19, 8, 1);
        mono_image.plot(7, 0, 1);
        mono_image.plot(8, 1, 1);
        mono_image.plot(19, 18, 1);
        // SSD1306 / SH1106: 3 pages of 20 column bytes, LSB is the top row of a page
        const std::vector<uint8_t> pages = convert_all<constixel::MONO_1BIT_VERTICAL_PAGES>(mono_image);
        check(pages.size() == 3 * 20, __LINE__);
        check(pages[0] == 0x81 && pages[1] == 0x02 && pages[7] == 0x01 && pages[8] == 0x02 && pages[2] == 0x00,
              __LINE__);
        check(pages[20 + 19] == 0x01 && pages[20] == 0x00, __LINE__);
        check(pages[40 + 19] == 0x04 && pages[40] == 0x00, __LINE__);
        // a rectangle starts its own pages at its top row: rows 7 .. 14 from column 0
        const std::vector<uint8_t> part = convert_rect<constixel::MONO_1BIT_VERTICAL_PAGES>(mono_image, {0, 7, 20, 8});
        check(part.size() == 20 && part[0] == 0x01 && part[19] == 0x02 && part[1] == 0x00, __LINE__);
        // horizontal, LSB is the leftmost pixel: 3 bytes per row
        const std::vector<uint8_t> rows = convert_all<constixel::MONO_1BIT_HORIZONTAL_LSB>(mono_image);
        check(rows.size() == 19 * 3, __LINE__);
        check(rows[0] == 0x81 && rows[1] == 0x00 && rows[3] == 0x02 && rows[4] == 0x01, __LINE__);
        check(rows[8 * 3 + 2] == 0x08 && rows[18 * 3 + 2] == 0x08, __LINE__);
        // the existing MSB-first layout for comparison
        const std::vector<uint8_t> raw = convert_all<constixel::STRAIGHT_THROUGH>(mono_image);
        check(raw[0] == 0x81 && raw[3] == 0x40 && raw[4] == 0x80, __LINE__);
        check_chunks<constixel::MONO_1BIT_VERTICAL_PAGES>(mono_image, {0, 0, 20, 19});
        check_chunks<constixel::MONO_1BIT_VERTICAL_PAGES>(mono_image, {3, 2, 13, 11});
        check_chunks<constixel::MONO_1BIT_HORIZONTAL_LSB>(mono_image, {3, 2, 13, 11});
        check_chunks<constixel::STRAIGHT_THROUGH>(mono_image, {3, 2, 13, 11});
        // dynamic_image has the same functions
        constixel::dynamic_image<constixel::pixel_1bit> dynamic(20, 19);
        check(dynamic.copy(std::span<const uint8_t>{mono_image.data_ref()}), __LINE__);
        check(convert_all<constixel::MONO_1BIT_VERTICAL_PAGES>(dynamic) == pages, __LINE__);
    }
}

}  // namespace

int main() {
    test_image_basics();
    test_string_view();
    test_draw_string();
    test_convert();
    return failures == 0 ? 0 : 1;
}
