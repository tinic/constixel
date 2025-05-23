#include "constixel.hpp"

int main() {
    static constixel::image<constixel::format_8bit, 256, 256> image;

    for (int32_t y = 0; y < 16; y++) {
        for (int32_t x = 0; x < 16; x++) {
            image.rect(x * 16, y * 16, 16, 16).fill(static_cast<uint8_t>(y * 16 + x));
        }
    }

    image.sixel_to_cout();
}
