#include "constixel.hpp"

int main() {
    static constixel::image<constixel::format_8bit, 256, 256> image;

    for (int32_t y = 0; y < 16; y++) {
        for (int32_t x = 0; x < 16; x++) {
            image.fill_rect({.x = x * 16, .y = y * 16, .w = 16, .h = 16, .col = static_cast<uint8_t>(y * 16 + x)});
        }
    }

    image.sixel_to_cout();
}
