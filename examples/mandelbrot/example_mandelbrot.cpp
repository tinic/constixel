#include "constixel.hpp"
#include "fonts/ibmplexsans_bold_32_aa.hpp"
#include <fstream>
#include <cmath>

using namespace constixel;
using UI = image<format_32bit, 640, 480>;
using myfont = ibmplexsans_bold_32_aa;

int main() {
    UI img;
    img.clear();

    // Shortcut to raw buffer
    auto *buf = img.data_ref().data();

    // Mandelbrot params
    const int W = 640, H = 480;
    const double xmin = -2.5, xmax = 1.0;
    const double ymin = -1.5, ymax = 1.5;
    const int maxIter = 500;

    // Draw fractal by direct BGRA writes
    for (int py = 0; py < H; ++py) {
        double y0 = ymin + (ymax - ymin) * py / (H - 1);
        for (int px = 0; px < W; ++px) {
            double x0 = xmin + (xmax - xmin) * px / (W - 1);
            double x = 0.0, y = 0.0;
            int iter = 0;
            while (x*x + y*y <= 4.0 && iter < maxIter) {
                double xt = x*x - y*y + x0;
                y = 2*x*y + y0;
                x = xt;
                ++iter;
            }

            // compute pixel index
            int idx = (py * W + px) * 4;

            if (iter == maxIter) {
                // black interior
                buf[idx+0] = 0;   // B
                buf[idx+1] = 0;   // G
                buf[idx+2] = 0;   // R
            } else {
                // smooth coloring: map iteration to hue
                double t = static_cast<double>(iter) / maxIter;
                uint8_t r = static_cast<uint8_t>(9*(1-t)*t*t*t*255);
                uint8_t g = static_cast<uint8_t>(15*(1-t)*(1-t)*t*t*255);
                uint8_t b = static_cast<uint8_t>(8.5*(1-t)*(1-t)*(1-t)*t*255);
                buf[idx+0] = b;
                buf[idx+1] = g;
                buf[idx+2] = r;
            }
            buf[idx+3] = 255;  // alpha
        }
    }

    // Draw a thin white border
    uint8_t cW = img.get_nearest_color(255,255,255);
    img.stroke_round_rect_aa(0, 0, W, H, 0, cW, 2);

    // Overlay title
    uint8_t cT = img.get_nearest_color(230,230,230);
    img.draw_string_centered_aa<myfont>(W/2, 32, "Mandelbrot Set", cT);

    // Write PNG
    std::ofstream fout("example_mandelbrot.png", std::ios::binary);
    img.png([&](char ch){ fout.put(ch); });

    return 0;
}
