#include "constixel.h"

#include "fonts/ibmplexsans_semibold_24_aa.h"
using large_font = constixel::ibmplexsans_semibold_24_aa;
#include "fonts/ibmplexsans_semibold_18_aa.h"
using small_font = constixel::ibmplexsans_semibold_18_aa;

#include <iomanip>
#include <sstream>

const char *x_label = "time (s)";
const char *y_label = "amplitude";

const float axis_xmin = 0.0f;
const float axis_xmax = 6.0f;
const float axis_ymin = -1.0f;
const float axis_ymax = +1.0f;

const uint8_t background_color = constixel::color::BLACK;
const uint8_t axis_color = constixel::color::WHITE;
const uint8_t graph_color = constixel::color::RED;

auto generate_damped_sine(float x_start, float x_end, float step, float decay_factor) {
    std::vector<std::pair<float, float>> points;
    for (float x = x_start; x <= x_end; x += step) {
        float y = std::exp(-decay_factor * x) * std::sin(5.0f * x);
        points.emplace_back(x, y);
    }
    return points;
}

auto generate_plot() {
    return generate_damped_sine(0.0f, 6.0f, 0.02f, 0.5f);
}

// ---------------------------

auto draw_graph(const char *x_label, const char *y_label, float x_min, float x_max, float y_min, float y_max,
                const std::vector<std::pair<float, float>> &points, uint8_t axis_col, uint8_t plot_col) {
    constixel::image<constixel::format_8bit, 1024, 1024, 1> img;

    const int32_t left_margin = 40 + img.string_width<large_font>(y_label);
    const int32_t right_margin = 40;
    const int32_t top_margin = 40;
    const int32_t bottom_margin = 40 + large_font::size + small_font::size;

    int32_t graph_x0 = left_margin;
    int32_t graph_y0 = top_margin;
    int32_t graph_w = img.width() - (left_margin + right_margin);
    int32_t graph_h = img.height() - (top_margin + bottom_margin);

    img.fill_rect(0, 0, img.width(), img.height(), background_color);

    img.draw_line_aa(graph_x0, graph_y0, graph_x0, graph_y0 + graph_h, axis_col);
    img.draw_line_aa(graph_x0, graph_y0 + graph_h, graph_x0 + graph_w, graph_y0 + graph_h, axis_col);

    img.draw_string_aa<large_font>(graph_x0 + graph_w / 2 - img.string_width<large_font>(x_label) / 2, graph_y0 + graph_h + 20 + small_font::size, x_label, axis_col);
    img.draw_string_aa<large_font>(graph_x0 - img.string_width<large_font>(y_label) - 10, graph_y0 + graph_h / 2, y_label, axis_col);

    auto scale_x = [&](float x) {
        return graph_x0 + (x - x_min) / (x_max - x_min) * graph_w;
    };

    auto scale_y = [&](float y) {
        return graph_y0 + graph_h - (y - y_min) / (y_max - y_min) * graph_h;
    };

    for (size_t i = 1; i < points.size(); ++i) {
        int32_t x0 = scale_x(points[i - 1].first);
        int32_t y0 = scale_y(points[i - 1].second);
        int32_t x1 = scale_x(points[i].first);
        int32_t y1 = scale_y(points[i].second);

        img.draw_line_aa(x0, y0, x1, y1, plot_col);
    }

    auto convert_to_string = [](float v) {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << v;
        return stream.str();
    };

    for (int32_t i = 0; i <= 5; ++i) {
        int32_t tx = graph_x0 + i * graph_w / 5;
        img.draw_line_aa(tx, graph_y0 + graph_h - 3, tx, graph_y0 + graph_h + 3, axis_col);        
        std::string x_val = convert_to_string(x_min + ( static_cast<float>(i) / 5.0f) * (x_max - x_min));
        int32_t wx = img.string_width<small_font>(x_val.c_str());
        img.draw_string_aa<small_font>(tx - wx/2, graph_y0 + graph_h + 8, x_val.c_str(), axis_col);
        int32_t ty = graph_y0 + graph_h - i * graph_h / 5;
        std::string y_val = convert_to_string(y_min + ( static_cast<float>(i) / 5.0f) * (y_max - y_min));
        int32_t wy = img.string_width<small_font>(y_val.c_str());
        img.draw_line_aa(graph_x0 - 3, ty, graph_x0 + 3, ty, axis_col);
        img.draw_string_aa<small_font>(graph_x0 - 8 - wy, ty - small_font::size / 2, y_val.c_str(), axis_col);
    }
    return img;
}

int main() {
    static auto plot_image = draw_graph(x_label, y_label, axis_xmin, axis_xmax, axis_ymin, axis_ymax, generate_plot(), axis_color, graph_color);
    plot_image.sixel_to_cout();
}
