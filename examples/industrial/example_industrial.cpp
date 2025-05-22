#include "constixel.hpp"
#include "fonts/ibmplexmono_bold_18_aa.hpp"
using myfont = constixel::ibmplexmono_bold_18_aa;
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <fstream>

using namespace constixel;
using UI = image<format_8bit, 640, 480>;

void drawStatusLED(UI &img, int cx, int cy, bool ok) {
    img.fill_circle_aa(cx, cy, 12, ok ? color::GREEN : color::RED);
    img.stroke_circle_aa(cx, cy, 12, color::WHITE, 2);
}

void drawGauge(UI &img, int x, int y, int w, int h, float pct, uint8_t bgCol, uint8_t fgCol) {
    // background track
    img.fill_round_rect_aa(x, y, w, h, h/2, bgCol);
    // filled portion
    int fw = int((w - 4) * pct);
    if (fw > 0)
        img.fill_round_rect_aa(x + 2, y + 2, fw, h - 4, (h - 4)/2, fgCol);
}

std::string nowString() {
    auto t = std::chrono::system_clock::now();
    std::time_t tt = std::chrono::system_clock::to_time_t(t);
    std::tm tm {}; 
#ifndef _MSC_VER
    localtime_r(&tt, &tm);
#endif  // #ifndef _MSC_VER
    std::ostringstream ss;
    ss << std::put_time(&tm, "%H:%M:%S");
    return ss.str();
}

int main() {
    UI img;
    img.clear();

    // 1) Top banner
    img.fill_round_rect_aa(0, 0, 640, 50, 0, color::DARK_BLUE);
    img.draw_string_aa<myfont>(20, 32-myfont::ascent, "Industrial Controller", color::WHITE);
    auto ts = nowString();
    img.draw_string_aa<myfont>(480, 32-myfont::ascent, ts.c_str(), color::GRAY_60);

    // 2) Motor panels
    const int xs[3] = {  15, 220, 425 };
    for (int i = 0; i < 3; ++i) {
        int x0 = xs[i], y0 = 70;
        bool ok = (i != 1);
        float curr = 0.3f + 0.2f * i;
        float volt = 0.4f + 0.1f * i;

        // panel background & border
        img.fill_round_rect_aa(x0+1, y0+1, 200-2, 200-2, 12, color::GRAY_60);
        img.stroke_round_rect_aa(x0, y0, 200, 200, 12, color::GRAY_40, 3);

        // label
        char lbl[3] = { char('1' + i), 0 };
        img.draw_string_aa<myfont>(x0 + 12, y0 + 12, lbl, color::WHITE);

        // status LED
        drawStatusLED(img, x0 + 180, y0 + 24, ok);

        // numeric readouts
        char buf[32];
        std::snprintf(buf, sizeof(buf), "I: %5.2f A", curr * 10.0f);
        img.draw_string_aa<myfont>(x0 + 12, y0 + 60, buf, color::WHITE);
        std::snprintf(buf, sizeof(buf), "V: %5.2f V", volt * 100.0f);
        img.draw_string_aa<myfont>(x0 + 12, y0 + 92, buf, color::WHITE);

        // gauges
        drawGauge(img, x0 + 12, y0 + 130, 176, 20, curr, color::GRAY_40, color::DARK_RED);
        drawGauge(img, x0 + 12, y0 + 165, 176, 20, volt, color::GRAY_40, color::BLUE);
    }

    // 3) Alarms strip
    img.fill_round_rect_aa(0, 430, 640, 50, 0, color::DARK_RED);
    img.draw_string_aa<myfont>(20, 440,
        "FAULT: Motor 2 Overcurrent", color::WHITE);

    // 4) Write out PNG
    std::ofstream fout("example_industrial.png", std::ios::binary);
    img.png([&](char c){ fout.put(c); });

    return 0;
}