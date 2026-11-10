#define CONSTIXEL_ENABLE_COUT

#include <algorithm>
#include <array>
#include <cmath>

#include "../example_output.hpp"
#include "../example_style.hpp"

using namespace constixel;
using namespace example;
using Canvas = image<format_32bit, 1024, 848>;

int main(int argc, char** argv) {
    static Canvas img;
    const Theme theme(img);
    header(img, "The Mandelbrot set", "An infinite boundary, drawn one pixel at a time", theme);
    card(img, 32, 108, 960, 676, theme);

    constexpr int left = 56, top = 132, width = 912, height = 592;
    constexpr int max_iterations = 500;
    constexpr double xmin = -2.6, xmax = 1.1;
    // Use the same scale on both axes so the set is not stretched.
    constexpr double scale = (xmax - xmin) / (width - 1);
    constexpr double ymax = scale * (height - 1) / 2;
    constexpr std::array<std::array<double, 3>, 6> colors = {
        {{18, 29, 48}, {31, 73, 113}, {46, 133, 151}, {172, 216, 202}, {249, 222, 169}, {222, 127, 89}}};
    auto& pixels = img.data_ref();
    for (int py = 0; py < height; ++py) {
        const double ci = ymax - py * scale;
        for (int px = 0; px < width; ++px) {
            const double cr = xmin + px * scale;
            double real = 0, imaginary = 0;
            int iteration = 0;
            while (real * real + imaginary * imaginary <= 4 && iteration < max_iterations) {
                const double next_real = real * real - imaginary * imaginary + cr;
                imaginary = 2 * real * imaginary + ci;
                real = next_real;
                ++iteration;
            }
            std::array<uint8_t, 3> rgb{13, 22, 34};
            if (iteration < max_iterations) {
                const double magnitude = std::sqrt(real * real + imaginary * imaginary);
                const double smooth = iteration + 1 - std::log2(std::log2(magnitude));
                const double t = std::clamp(std::log1p(std::max(0.0, smooth)) / std::log(120.0), 0.0, 1.0);
                const double position = t * (colors.size() - 1);
                const size_t first = std::min(static_cast<size_t>(position), colors.size() - 2);
                const double blend = position - static_cast<double>(first);
                for (size_t channel = 0; channel < 3; ++channel) {
                    rgb[channel] =
                        static_cast<uint8_t>(std::lerp(colors[first][channel], colors[first + 1][channel], blend));
                }
            }
            // format_32bit stores R, G, B, A bytes in this order.
            const size_t offset = (static_cast<size_t>(top + py) * img.width() + left + px) * 4;
            pixels[offset] = rgb[0];
            pixels[offset + 1] = rgb[1];
            pixels[offset + 2] = rgb[2];
            pixels[offset + 3] = 255;
        }
    }
    text<heading>(img, 56, 760, "z(n + 1) = z(n)^2 + c", theme.ink);
    right_text<small>(img, 968, 759, "500 iterations / smooth escape-time coloring", theme.muted);
    text<small>(img, 32, 821, "Real axis: -2.6 to 1.1", theme.muted);
    right_text<small>(img, 992, 821, "Equal axis scale / direct RGBA buffer writes", theme.muted);
    return output(img, argc, argv, "example_mandelbrot.png");
}
