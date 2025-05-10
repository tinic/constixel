#include "constixel.h"
#include "fonts/interdisplay_bold_48_aa.h"

static consteval auto create_static_panel_image() {
    using namespace constixel;

    image<format_8bit, 300, 256, 1> image;

    using myfont = interdisplay_bold_48_aa;
    const int32_t buttonHeight = 75;
    const int32_t textCenterY = ( buttonHeight + myfont::ascent) / 2 - myfont::line_height + 3;

    int32_t buttonY = 0;
    const std::string buttonAStr("Start");
    image.fill_round_rect_aa(0, buttonY, image.width(), 75, 20, image.get_nearest_color(0x2F,0x7F,0x2F));
    int32_t widthA = image.string_width<myfont>(buttonAStr.c_str());
    image.draw_string_aa<myfont>((image.width() - widthA)/2, textCenterY, buttonAStr.c_str(), color::WHITE);

    buttonY += 86;
    const std::string buttonBStr("Pause");
    image.fill_round_rect_aa(0, buttonY, image.width(), 75, 20, image.get_nearest_color(0x2F,0x2F,0x7F));
    int32_t widthB = image.string_width<myfont>(buttonBStr.c_str());
    image.draw_string_aa<myfont>((image.width() - widthB)/2, buttonY + textCenterY, buttonBStr.c_str(), color::WHITE);

    buttonY += 86;
    const std::string buttonCStr("STOP");
    image.fill_round_rect_aa(0, buttonY, image.width(), 75, 20, image.get_nearest_color(0x7F,0xF,0xF));
    int32_t widthC = image.string_width<myfont>(buttonCStr.c_str());
    image.draw_string_aa<myfont>((image.width() - widthC)/2, buttonY + textCenterY, buttonCStr.c_str(), color::WHITE);

    return image;
}

int main() {
    static constexpr auto panel_image = create_static_panel_image();
    panel_image.sixel_to_cout();
}
