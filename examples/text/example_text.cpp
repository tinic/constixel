#define CONSTIXEL_ENABLE_COUT

#include <array>

#include "../example_output.hpp"
#include "../example_style.hpp"

using namespace constixel;
using namespace example;

// The finished panel, including antialiased type, is evaluated by the compiler.
static consteval auto create_static_panel_image() {
    image<format_8bit, 440, 344> img;
    const Theme theme(img);
    img.fill_rect(0, 0, img.width(), img.height(), theme.background);
    text<title>(img, 28, 49, "Session controls", theme.ink);
    text<small>(img, 28, 76, "Ready when you are.", theme.muted);

    constexpr std::array<const char*, 3> labels = {"Start session", "Pause session", "Stop session"};
    constexpr std::array<const char*, 3> shortcuts = {"01", "02", "03"};
    for (int i = 0; i < 3; ++i) {
        const int y = 100 + i * 68;
        const auto background = i == 0 ? theme.accent : theme.surface;
        const auto foreground = i == 0 ? theme.surface : (i == 2 ? theme.danger : theme.ink);
        img.fill_round_rect_aa(28, y, 384, 56, 12, background);
        text<heading>(img, 52, y + 35, labels[i], foreground);
        right_text<small>(img, 388, y + 34, shortcuts[i], foreground);
    }
    text<caption>(img, 28, 324, "CONSTEXPR / TYPE + SHAPES", theme.muted);
    return img;
}

int main(int argc, char** argv) {
    static constexpr auto panel_image = create_static_panel_image();
    return output(panel_image, argc, argv);
}
