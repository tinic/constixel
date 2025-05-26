#include "constixel.hpp"

#include "fonts/ibmplexsans_bold_18_aa.hpp"
using titlefont = constixel::ibmplexsans_bold_18_aa;

#include "fonts/ibmplexsans_bold_14_aa.hpp"
using labelfont = constixel::ibmplexsans_bold_14_aa;

#include "fonts/ibmplexsans_regular_12_aa.hpp"
using textfont = constixel::ibmplexsans_regular_12_aa;

int main() {
    using namespace constixel;

    static image<format_8bit, 800, 800> image;

    // Title
    image.draw_string_aa<titlefont>(50, 30, "Dashed Lines Showcase", color::WHITE);

    int32_t y_pos = 80;
    int32_t x_margin = 60;
    int32_t spacing = 90;

    // 1. Basic Patterns
    {
        image.draw_string_aa<labelfont>(x_margin, y_pos, "Basic Patterns", color::CYAN);
        y_pos += 30;

        std::array<float, 16> simple = {20.0f, 10.0f};
        image.draw_dashed_line_aa(x_margin, y_pos, x_margin + 400, y_pos, color::WHITE, 3.0f, simple, 2);
        image.draw_string_aa<textfont>(x_margin + 420, y_pos - 5, "Simple: dash-gap", color::GRAY_60);

        y_pos += 25;
        std::array<float, 16> complex = {25.0f, 8.0f, 5.0f, 8.0f, 15.0f, 8.0f};
        image.draw_dashed_line_aa(x_margin, y_pos, x_margin + 400, y_pos, color::YELLOW, 3.0f, complex, 6);
        image.draw_string_aa<textfont>(x_margin + 420, y_pos - 5, "Complex: long-short-dot pattern", color::GRAY_60);

        y_pos += 25;
        std::array<float, 16> morse = {12.0f, 4.0f, 12.0f, 4.0f, 12.0f, 12.0f, 30.0f, 4.0f, 30.0f, 4.0f, 30.0f, 12.0f};
        image.draw_dashed_line_aa(x_margin, y_pos, x_margin + 400, y_pos, color::GREEN, 3.0f, morse, 12);
        image.draw_string_aa<textfont>(x_margin + 420, y_pos - 5, "Morse: SOS pattern", color::GRAY_60);

        y_pos += spacing;
    }

    // 2. Line Weights
    {
        image.draw_string_aa<labelfont>(x_margin, y_pos, "Variable Thickness", color::CYAN);
        y_pos += 30;

        std::array<float, 16> pattern = {18.0f, 10.0f};
        for (int i = 0; i < 5; i++) {
            float width = 1.0f + i * 1.2f;
            image.draw_dashed_line_aa(x_margin, y_pos + i * 20, x_margin + 350, y_pos + i * 20, 
                                     color::RED, width, pattern, 2);
            
            char label[10];
            snprintf(label, sizeof(label), "%.1fpx", width);
            image.draw_string_aa<textfont>(x_margin + 370, y_pos + i * 20 - 5, label, color::GRAY_60);
        }

        y_pos += spacing + 20;
    }

    // 3. Connected Paths
    {
        image.draw_string_aa<labelfont>(x_margin, y_pos, "Seamless Connections", color::CYAN);
        y_pos += 30;

        std::array<float, 16> pattern = {20.0f, 8.0f, 5.0f, 8.0f};
        float offset = 0.0f;

        // Create a path that shows pattern continuity
        int32_t path[][2] = {
            {x_margin, y_pos}, {x_margin + 120, y_pos},
            {x_margin + 120, y_pos + 25}, {x_margin + 250, y_pos + 25},
            {x_margin + 250, y_pos + 50}, {x_margin + 80, y_pos + 50},
            {x_margin + 80, y_pos + 75}, {x_margin + 200, y_pos + 75}
        };

        for (int i = 0; i < 7; i++) {
            offset = image.draw_dashed_line_aa(path[i][0], path[i][1], 
                                              path[i+1][0], path[i+1][1],
                                              color::MAGENTA, 4.0f, pattern, 4, offset);
        }

        image.draw_string_aa<textfont>(x_margin + 280, y_pos + 30, "Pattern flows smoothly", color::GRAY_60);
        image.draw_string_aa<textfont>(x_margin + 280, y_pos + 45, "across line segments", color::GRAY_60);

        y_pos += spacing + 10;
    }

    // 4. Creative Applications
    {
        image.draw_string_aa<labelfont>(x_margin, y_pos, "Creative Uses", color::CYAN);
        y_pos += 30;

        // Dashed border
        std::array<float, 16> border_pattern = {12.0f, 6.0f, 3.0f, 6.0f};
        int32_t box_x = x_margin;
        int32_t box_y = y_pos;
        int32_t box_w = 140;
        int32_t box_h = 80;
        
        float border_offset = 0.0f;
        border_offset = image.draw_dashed_line_aa(box_x, box_y, box_x + box_w, box_y, 
                                                 color::BLUE, 3.0f, border_pattern, 4, border_offset);
        border_offset = image.draw_dashed_line_aa(box_x + box_w, box_y, box_x + box_w, box_y + box_h, 
                                                 color::BLUE, 3.0f, border_pattern, 4, border_offset);
        border_offset = image.draw_dashed_line_aa(box_x + box_w, box_y + box_h, box_x, box_y + box_h, 
                                                 color::BLUE, 3.0f, border_pattern, 4, border_offset);
        border_offset = image.draw_dashed_line_aa(box_x, box_y + box_h, box_x, box_y, 
                                                 color::BLUE, 3.0f, border_pattern, 4, border_offset);

        image.fill_rect(box_x + 8, box_y + 8, box_w - 16, box_h - 16, color::GRAY_20);
        image.draw_string_centered_aa<textfont>(box_x + box_w/2, box_y + box_h/2 - 8, "Dashed", color::WHITE);
        image.draw_string_centered_aa<textfont>(box_x + box_w/2, box_y + box_h/2 + 8, "Border", color::WHITE);

        // Grid pattern
        std::array<float, 16> grid_pattern = {6.0f, 8.0f};
        int32_t grid_x = x_margin + 200;
        
        for (int i = 0; i <= 6; i++) {
            image.draw_dashed_line_aa(grid_x + i * 25, y_pos, grid_x + i * 25, y_pos + 80, 
                                     color::GRAY_40, 1.5f, grid_pattern, 2);
            if (i <= 3) {
                image.draw_dashed_line_aa(grid_x, y_pos + i * 20, grid_x + 150, y_pos + i * 20, 
                                         color::GRAY_40, 1.5f, grid_pattern, 2);
            }
        }

        image.draw_string_aa<textfont>(grid_x + 160, y_pos + 20, "Dashed Grid", color::GRAY_60);

        // Wave pattern overlaid on the grid
        std::array<float, 16> wave_pattern = {4.0f, 3.0f};
        for (int x = 0; x < 150; x += 3) {
            int32_t y1 = y_pos + 40 + static_cast<int32_t>(25 * sin(x * 0.04));
            int32_t y2 = y_pos + 40 + static_cast<int32_t>(25 * sin((x + 3) * 0.04));
            image.draw_dashed_line_aa(grid_x + x, y1, grid_x + x + 3, y2, 
                                     color::GREEN, 2.5f, wave_pattern, 2);
        }

        image.draw_string_aa<textfont>(grid_x + 160, y_pos + 40, "Dashed Grid", color::GRAY_60);
        image.draw_string_aa<textfont>(grid_x + 160, y_pos + 55, "with Sine Wave", color::GRAY_60);

        // Separate example area for other demonstrations
        int32_t demo_x = x_margin + 400;

        y_pos += spacing + 10;
    }

    // 5. API Styles
    {
        image.draw_string_aa<labelfont>(x_margin, y_pos, "Multiple API Styles", color::CYAN);
        y_pos += 30;

        std::array<float, 16> pattern = {15.0f, 7.0f};

        // Direct function call
        image.draw_dashed_line_aa(x_margin, y_pos, x_margin + 200, y_pos, color::WHITE, 2.5f, pattern, 2);
        image.draw_string_aa<textfont>(x_margin + 220, y_pos - 5, "Direct function call", color::GRAY_60);

        y_pos += 25;
        
        // Struct-based
        draw_dashed_line line_struct = {x_margin, y_pos, x_margin + 200, y_pos, color::CYAN, 2.5f, pattern, 2, 0.0f};
        image.draw_dashed_line_aa(line_struct);
        image.draw_string_aa<textfont>(x_margin + 220, y_pos - 5, "Struct-based call", color::GRAY_60);

        y_pos += 25;

        // Fluent API
        image.dashed_line_aa(x_margin, y_pos, x_margin + 200, y_pos, pattern, 2, 0.0f)
             .stroke(color::YELLOW, 2.5f);
        image.draw_string_aa<textfont>(x_margin + 220, y_pos - 5, "Fluent API", color::GRAY_60);
    }

    image.sixel_to_cout<
#ifdef __APPLE__
    2
#endif  // #ifdef __APPLE__
    >();
}