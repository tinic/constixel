#include "constixel.h"

#include <format>

//
// Font to show a specimen of
//
#include "fonts/cd_iconspc_mono.h"
using font = constixel::cd_iconspc_mono;

//--------------------------------------------------------------------

// For for label
#include "fonts/pxplus_hp_100lx_6x8_mono.h"
using label = constixel::pxplus_hp_100lx_6x8_mono;

// Function to convert UTF-32 to UTF-8
static void append_utf8(std::string& out, char32_t cp) {
    if (cp <= 0x7F) {
        out.push_back(static_cast<char>(cp));
    } else if (cp <= 0x7FF) {
        out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
        out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else if (cp <= 0xFFFF) {
        out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
        out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else if (cp <= 0x10FFFF) {
        out.push_back(static_cast<char>(0xF0 | (cp >> 18)));
        out.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    }
}

int main() {
    static constixel::image<constixel::format_1bit, 512, 600, 2> screen;

    int32_t y_pos = 8;
    int32_t x_pos = 8;

    for (auto c : font::glyph_table) {

        // Get a UTF-8 string with our character in it.
        std::string str{};
        append_utf8(str, c.first);

        // Get the width of the character
        int32_t cw = screen.string_width<font>(str.c_str()) + 2;
        cw = std::max(32, cw);

        // Draw character
        screen.draw_string_centered_mono<font>(x_pos + cw / 2, y_pos, str.c_str(), 1);

        // Draw label under it
        auto unicode = std::format("{:04x}", c.first);
        screen.draw_string_centered_mono<label>(x_pos + cw / 2, y_pos + font::total_height + 4, unicode.c_str(), 1);

        // Advance and break line
        x_pos += cw;
        if (x_pos > 460) {
            y_pos += font::total_height + 18;
            x_pos = 8;
        }
    }

    screen.sixel_to_cout();
}
