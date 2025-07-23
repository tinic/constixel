#define CONSTIXEL_ENABLE_COUT

#include "constixel.hpp"
#include "fonts/ibmplexsans_semibold_24_aa.hpp"
#include "fonts/ibmplexsans_semibold_18_aa.hpp"
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace constixel;
using large_font = ibmplexsans_semibold_24_aa;
using small_font = ibmplexsans_semibold_18_aa;

constexpr char x_label[] = "Time (seconds)";
constexpr char y_label[] = "Amplitude";
constexpr float axis_xmin = 0.0f;
constexpr float axis_xmax = 6.0f;
constexpr float axis_ymin = -1.2f;
constexpr float axis_ymax = +1.2f;

const uint8_t col_bg_dark = color::BLACK;
const uint8_t col_bg_light = 8;
const uint8_t col_axis = color::WHITE;
const uint8_t col_grid = 7;
const uint8_t col_graph_primary = color::CYAN;
const uint8_t col_graph_secondary = color::YELLOW;
const uint8_t col_text = color::WHITE;


constexpr int bayer_matrix[4][4] = {
    { 0,  8,  2, 10},
    {12,  4, 14,  6},
    { 3, 11,  1,  9},
    {15,  7, 13,  5}
};

auto generate_damped_sine(float x_start, float x_end, float step, float decay_factor) {
    std::vector<std::pair<float, float>> points;
    for (float x = x_start; x <= x_end; x += step) {
        float y = std::exp(-decay_factor * x) * std::sin(5.0f * x);
        points.emplace_back(x, y);
    }
    return points;
}

auto generate_secondary_wave(float x_start, float x_end, float step) {
    std::vector<std::pair<float, float>> points;
    for (float x = x_start; x <= x_end; x += step) {
        float y = 0.6f * std::sin(2.0f * x + 1.0f) * std::cos(0.5f * x);
        points.emplace_back(x, y);
    }
    return points;
}

auto generate_primary_data() {
    return generate_damped_sine(axis_xmin, axis_xmax, 0.015f, 0.4f);
}

auto generate_secondary_data() {
    return generate_secondary_wave(axis_xmin, axis_xmax, 0.02f);
}

static image<format_8bit, 1024, 1024> img;

void draw_gradient_background() {
    for (int32_t y = 0; y < img.height(); ++y) {
        for (int32_t x = 0; x < img.width(); ++x) {
            float gradient = static_cast<float>(y) / img.height();
            int dither_threshold = bayer_matrix[y % 4][x % 4];
            float dither_value = gradient * 16.0f;
            
            uint8_t color = (dither_value > dither_threshold) ? col_bg_light : col_bg_dark;
            img.plot(x, y, color);
        }
    }
}

void draw_grid(int32_t x0, int32_t y0, int32_t w, int32_t h) {
    for (int i = 1; i < 6; ++i) {
        int32_t x = x0 + i * w / 6;
        for (int32_t y = y0; y < y0 + h; y += 3) {
            img.plot(x, y, col_grid);
        }
    }
    
    for (int i = 1; i < 6; ++i) {
        int32_t y = y0 + i * h / 6;
        for (int32_t x = x0; x < x0 + w; x += 3) {
            img.plot(x, y, col_grid);
        }
    }
}

auto draw_enhanced_graph() {
    const int32_t left_margin = 120;
    const int32_t right_margin = 60;
    const int32_t top_margin = 60;
    const int32_t bottom_margin = 80 + large_font::size + small_font::size;

    const int32_t graph_x0 = left_margin;
    const int32_t graph_y0 = top_margin;
    const int32_t graph_w = img.width() - (left_margin + right_margin);
    const int32_t graph_h = img.height() - (top_margin + bottom_margin);

    draw_gradient_background();
    
    draw_grid(graph_x0, graph_y0, graph_w, graph_h);

    for (int i = 0; i < 2; ++i) {
        img.draw_line_aa(graph_x0 - i, graph_y0, graph_x0 - i, graph_y0 + graph_h, col_axis);
        img.draw_line_aa(graph_x0, graph_y0 + graph_h + i, graph_x0 + graph_w, graph_y0 + graph_h + i, col_axis);
    }


    img.draw_string_aa<large_font>(
        graph_x0 + graph_w / 2 - img.string_width<large_font>(x_label) / 2, 
        graph_y0 + graph_h + 35 + small_font::size, 
        x_label, col_text);
    img.draw_string_aa<large_font, false, DEGREE_270>(
        graph_x0 - 80, 
        graph_y0 + graph_h / 2 + img.string_width<large_font>(y_label) / 2, 
        y_label, col_text);

    auto scale_x = [&](float x) -> int32_t {
        return graph_x0 + static_cast<int32_t>((x - axis_xmin) / (axis_xmax - axis_xmin) * graph_w);
    };

    auto scale_y = [&](float y) -> int32_t {
        return graph_y0 + graph_h - static_cast<int32_t>((y - axis_ymin) / (axis_ymax - axis_ymin) * graph_h);
    };

    auto primary_data = generate_primary_data();
    auto secondary_data = generate_secondary_data();

    for (size_t i = 1; i < primary_data.size(); ++i) {
        int32_t x0 = scale_x(primary_data[i - 1].first);
        int32_t y0 = scale_y(primary_data[i - 1].second);
        int32_t x1 = scale_x(primary_data[i].first);
        int32_t y1 = scale_y(primary_data[i].second);
        img.draw_line_aa(x0, y0, x1, y1, col_graph_primary);
    }

    for (size_t i = 1; i < secondary_data.size(); ++i) {
        int32_t x0 = scale_x(secondary_data[i - 1].first);
        int32_t y0 = scale_y(secondary_data[i - 1].second);
        int32_t x1 = scale_x(secondary_data[i].first);
        int32_t y1 = scale_y(secondary_data[i].second);
        img.draw_line_aa(x0, y0, x1, y1, col_graph_secondary);
    }


    auto format_value = [](float v) -> std::string {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << v;
        return stream.str();
    };

    for (int i = 0; i <= 6; ++i) {
        int32_t tx = graph_x0 + i * graph_w / 6;
        
        img.draw_line_aa(tx, graph_y0 + graph_h - 5, tx, graph_y0 + graph_h + 5, col_axis);
        
        float x_val = axis_xmin + (static_cast<float>(i) / 6.0f) * (axis_xmax - axis_xmin);
        std::string x_label_text = format_value(x_val);
        int32_t label_width = img.string_width<small_font>(x_label_text.c_str());
        img.draw_string_aa<small_font>(
            tx - label_width / 2, 
            graph_y0 + graph_h + 15, 
            x_label_text.c_str(), col_text);
    }

    for (int i = 0; i <= 6; ++i) {
        int32_t ty = graph_y0 + graph_h - i * graph_h / 6;
        
        img.draw_line_aa(graph_x0 - 5, ty, graph_x0 + 5, ty, col_axis);
        
        float y_val = axis_ymin + (static_cast<float>(i) / 6.0f) * (axis_ymax - axis_ymin);
        std::string y_label_text = format_value(y_val);
        int32_t label_width = img.string_width<small_font>(y_label_text.c_str());
        img.draw_string_aa<small_font>(
            graph_x0 - 15 - label_width, 
            ty - small_font::size / 2, 
            y_label_text.c_str(), col_text);
    }


    const int32_t legend_x = graph_x0 + graph_w - 200;
    const int32_t legend_y = graph_y0 + 20;
    
    img.draw_line_aa(legend_x, legend_y, legend_x + 30, legend_y, col_graph_primary);
    img.draw_string_aa<small_font>(legend_x + 35, legend_y - small_font::size / 2, "Damped Sine", col_text);
    
    img.draw_line_aa(legend_x, legend_y + 25, legend_x + 30, legend_y + 25, col_graph_secondary);
    img.draw_string_aa<small_font>(legend_x + 35, legend_y + 25 - small_font::size / 2, "Modulated Wave", col_text);

    return img;
}

int main() {
    static auto enhanced_plot = draw_enhanced_graph();
    enhanced_plot.sixel_to_cout();
    return 0;
}
