#define CONSTIXEL_ENABLE_COUT

#include "constixel.hpp"

#include "fonts/ibmplexsans_bold_48_aa.hpp"
using myfont = constixel::ibmplexsans_bold_48_aa;

static consteval auto create_static_panel_image() {
    using namespace constixel;

    image<format_8bit, 300, 256> image;

    const int32_t buttonHeight = 75;
    const int32_t textCenterY = ( buttonHeight + myfont::ascent) / 2 - myfont::line_height + 3;

    int32_t buttonY = 0;
    image.fill_round_rect_aa(0, buttonY, image.width(), 75, 20, image.get_nearest_color(0x2F,0x7F,0x2F));
    image.stroke_round_rect_aa(0, buttonY, image.width(), 75, 20, color::WHITE, 2);
    image.draw_string_centered_aa<myfont,true>(image.width()/2, buttonY + textCenterY, "Start", color::WHITE);

    buttonY += 86;
    image.fill_round_rect_aa(0, buttonY, image.width(), 75, 20, image.get_nearest_color(0x2F,0x2F,0x7F));
    image.stroke_round_rect_aa(0, buttonY, image.width(), 75, 20, color::WHITE, 2);
    image.draw_string_centered_aa<myfont,true>(image.width()/2, buttonY + textCenterY, "Pause", color::WHITE);

    buttonY += 86;
    image.fill_round_rect_aa(0, buttonY, image.width(), 75, 20, image.get_nearest_color(0x7F,0xF,0xF));
    image.stroke_round_rect_aa(0, buttonY, image.width(), 75, 20, color::WHITE, 2);
    image.draw_string_centered_aa<myfont,true>(image.width()/2, buttonY + textCenterY, "STOP", color::WHITE);

    return image;
}

int main() {
    static constexpr auto panel_image = create_static_panel_image();
    panel_image.sixel_to_cout();
}
