#include "constixel.hpp"
#include "fonts/ibmplexsans_bold_18_aa.hpp"
#include "fonts/ibmplexsans_bold_12_aa.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace constixel;
using myfont = constixel::ibmplexsans_bold_18_aa;
using smallfont = constixel::ibmplexsans_bold_12_aa;
using UI = image<format_8bit, 1024, 768>;

// Industrial HMI color scheme (inspired by classic industrial monitors)
struct Colors {
    static uint8_t background;      // Very dark blue-black #0a0f1a
    static uint8_t panel_bg;        // Dark blue panel #1a2332
    static uint8_t panel_border;    // Medium blue border #2d3f54
    static uint8_t text_primary;    // White text #ffffff
    static uint8_t text_secondary;  // Light gray #cccccc
    static uint8_t accent_blue;     // Subdued blue #4080c0
    static uint8_t accent_green;    // Subdued green #60a060
    static uint8_t accent_cyan;     // Subdued cyan #5099aa
    static uint8_t accent_orange;   // Subdued orange #cc7744
    static uint8_t danger_red;      // Subdued red #cc5555
    static uint8_t warning_yellow;  // Subdued yellow #ccaa44
    static uint8_t success_green;   // Subdued green #77cc77
    static uint8_t button_bg;       // Button dark blue #334455
    static uint8_t button_hover;    // Button highlight #4a6080
    
    static void init(const UI& img) {
        background = img.get_nearest_color(10, 15, 26);      // #0a0f1a - Very dark blue-black
        panel_bg = img.get_nearest_color(26, 35, 50);        // #1a2332 - Dark blue panel
        panel_border = img.get_nearest_color(45, 63, 84);    // #2d3f54 - Medium blue border
        text_primary = img.get_nearest_color(255, 255, 255); // #ffffff - White text
        text_secondary = img.get_nearest_color(204, 204, 204); // #cccccc - Light gray
        accent_blue = img.get_nearest_color(64, 128, 192);   // #4080c0 - Subdued blue
        accent_green = img.get_nearest_color(96, 160, 96);   // #60a060 - Subdued green
        accent_cyan = img.get_nearest_color(80, 153, 170);   // #5099aa - Subdued cyan
        accent_orange = img.get_nearest_color(204, 119, 68); // #cc7744 - Subdued orange
        danger_red = img.get_nearest_color(204, 85, 85);     // #cc5555 - Subdued red
        warning_yellow = img.get_nearest_color(204, 170, 68); // #ccaa44 - Subdued yellow
        success_green = img.get_nearest_color(119, 204, 119); // #77cc77 - Subdued green
        button_bg = img.get_nearest_color(51, 68, 85);       // #334455 - Button dark blue
        button_hover = img.get_nearest_color(74, 96, 128);   // #4a6080 - Button highlight
    }
};

uint8_t Colors::background;
uint8_t Colors::panel_bg;
uint8_t Colors::panel_border;
uint8_t Colors::text_primary;
uint8_t Colors::text_secondary;
uint8_t Colors::accent_blue;
uint8_t Colors::accent_green;
uint8_t Colors::accent_cyan;
uint8_t Colors::accent_orange;
uint8_t Colors::danger_red;
uint8_t Colors::warning_yellow;
uint8_t Colors::success_green;
uint8_t Colors::button_bg;
uint8_t Colors::button_hover;

void drawStatusLED(UI &img, int cx, int cy, int radius, bool status, uint8_t on_color = 0) {
    uint8_t led_color = status ? (on_color ? on_color : Colors::success_green) : Colors::danger_red;
    img.fill_circle_aa(cx, cy, radius, led_color);
    img.stroke_circle_aa(cx, cy, radius, Colors::text_primary, 1);
    if (status) {
        img.fill_circle_aa(cx - radius/3, cy - radius/3, radius/4, Colors::text_primary);
    }
}

void drawProgressBar(UI &img, int x, int y, int w, int h, float pct, uint8_t bgCol, uint8_t fgCol, const char* label = nullptr) {
    img.fill_round_rect_aa(x, y, w, h, h/2, bgCol);
    img.stroke_round_rect_aa(x, y, w, h, h/2, Colors::panel_border, 1);
    
    int fw = int((w - 4) * pct);
    if (fw > 0) {
        uint8_t bar_color = fgCol;
        if (pct > 0.8f) bar_color = Colors::danger_red;
        else if (pct > 0.6f) bar_color = Colors::warning_yellow;
        
        img.fill_round_rect_aa(x + 2, y + 2, fw, h - 4, (h - 4)/2, bar_color);
    }
    
    if (label) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%s: %3.0f%%", label, pct * 100.0f);
        img.draw_string_aa<smallfont, true>(x + 5, y + h/2 + 5 - smallfont::ascent, buf, Colors::text_primary);
    }
}

void drawGauge(UI &img, int cx, int cy, int radius, float value, float min_val, float max_val, 
               const char* label, const char* unit, uint8_t needle_color = 0) {
    const float start_angle = 3.14159f * 1.25f;
    const float end_angle = 3.14159f * 0.25f;
    const float angle_range = end_angle - start_angle;
    
    img.fill_circle_aa(cx, cy, radius, Colors::panel_bg);
    img.stroke_circle_aa(cx, cy, radius, Colors::panel_border, 2);
    
    for (int i = 0; i <= 10; i++) {
        float angle = start_angle + (angle_range * i / 10.0f);
        int x1 = cx + (radius - 15) * cos(angle);
        int y1 = cy + (radius - 15) * sin(angle);
        int x2 = cx + (radius - 5) * cos(angle);
        int y2 = cy + (radius - 5) * sin(angle);
        img.draw_line_aa(x1, y1, x2, y2, Colors::text_secondary, 2);
    }
    
    float normalized = (value - min_val) / (max_val - min_val);
    float needle_angle = start_angle + (angle_range * normalized);
    int nx = cx + (radius - 20) * cos(needle_angle);
    int ny = cy + (radius - 20) * sin(needle_angle);
    uint8_t final_needle_color = needle_color ? needle_color : Colors::danger_red;
    img.draw_line_aa(cx, cy, nx, ny, final_needle_color, 4);
    
    img.fill_circle_aa(cx, cy, 8, Colors::panel_border);
    img.stroke_circle_aa(cx, cy, 8, Colors::text_primary, 1);
    
    img.draw_string_aa<smallfont, true>(cx - 30, cy + radius + 20 - smallfont::ascent, label, Colors::text_primary);
    char val_str[32];
    snprintf(val_str, sizeof(val_str), "%.1f %s", value, unit);
    img.draw_string_aa<smallfont, true>(cx - 25, cy + radius + 35 - smallfont::ascent, val_str, Colors::accent_cyan);
}

void drawEquipmentPanel(UI &img, int x, int y, int w, int h, const char* name, 
                       bool online, float temp, float pressure, float flow) {
    img.fill_round_rect_aa(x, y, w, h, 12, Colors::panel_bg);
    img.stroke_round_rect_aa(x, y, w, h, 12, Colors::panel_border, 1);
    
    img.draw_string_aa<myfont, true>(x + 10, y + 25 - myfont::ascent, name, Colors::text_primary);
    
    drawStatusLED(img, x + w - 25, y + 20, 8, online);
    
    char status_text[16];
    strcpy(status_text, online ? "ONLINE" : "OFFLINE");
    img.draw_string_aa<smallfont, true>(x + w - 80, y + 40 - smallfont::ascent, status_text, online ? Colors::success_green : Colors::danger_red);
    
    drawProgressBar(img, x + 10, y + 50, w - 20, 15, temp / 100.0f, Colors::background, Colors::accent_blue, "TEMP");
    drawProgressBar(img, x + 10, y + 75, w - 20, 15, pressure / 10.0f, Colors::background, Colors::accent_green, "PRES");
    drawProgressBar(img, x + 10, y + 100, w - 20, 15, flow / 50.0f, Colors::background, Colors::accent_cyan, "FLOW");
}

void drawAlarmPanel(UI &img, int x, int y, int w, int h) {
    bool alarm_active = true;
    uint8_t bg_color = alarm_active ? Colors::danger_red : Colors::panel_bg;
    
    img.fill_round_rect_aa(x, y, w, h, 8, bg_color);
    img.stroke_round_rect_aa(x, y, w, h, 8, alarm_active ? Colors::danger_red : Colors::panel_border, 2);
    
    if (alarm_active) {
        drawStatusLED(img, x + 20, y + h/2, 8, true, Colors::danger_red);
        img.draw_string_aa<myfont, true>(x + 40, y + h/2 - myfont::ascent, "CRITICAL: Pump 2 Overpressure", Colors::text_primary);
    } else {
        img.draw_string_aa<myfont, true>(x + 20, y + h/2 - myfont::ascent, "All Systems Normal", Colors::success_green);
    }
}

std::string getCurrentTime() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%H:%M:%S");
    return ss.str();
}

int main() {
    static UI img;
    img.clear();
    
    // Initialize our modern color scheme
    Colors::init(img);
    
    
    img.fill_rect(0, 0, 1024, 768, Colors::background);
    
    img.fill_round_rect_aa(0, 0, 1024, 60, 0, Colors::panel_bg);
    img.draw_string_aa<myfont, true>(20, 35 - myfont::ascent, "ACME Industries - Control Center", Colors::text_primary);
    
    auto time_str = getCurrentTime();
    img.draw_string_aa<myfont, true>(800, 35 - myfont::ascent, time_str.c_str(), Colors::text_secondary);
    
    drawEquipmentPanel(img, 20, 80, 240, 130, "Pump Station 1", true, 72.5f, 8.2f, 35.7f);
    drawEquipmentPanel(img, 280, 80, 240, 130, "Pump Station 2", false, 45.3f, 12.8f, 0.0f);
    drawEquipmentPanel(img, 540, 80, 240, 130, "Compressor A", true, 88.1f, 6.4f, 42.3f);
    drawEquipmentPanel(img, 800, 80, 200, 130, "Cooling Tower", true, 34.7f, 2.1f, 28.9f);
    
    drawGauge(img, 120, 300, 80, 72.5f, 0.0f, 100.0f, "Temperature", "°C", Colors::accent_orange);
    drawGauge(img, 320, 300, 80, 8.2f, 0.0f, 15.0f, "Pressure", "PSI", Colors::accent_blue);
    drawGauge(img, 520, 300, 80, 35.7f, 0.0f, 50.0f, "Flow Rate", "L/s", Colors::accent_cyan);
    drawGauge(img, 720, 300, 80, 2.4f, 0.0f, 5.0f, "Power", "MW", Colors::warning_yellow);
    drawGauge(img, 920, 300, 80, 94.2f, 0.0f, 100.0f, "Efficiency", "%", Colors::success_green);
    
    int graph_x = 50, graph_y = 450, graph_w = 400, graph_h = 150;
    img.fill_round_rect_aa(graph_x, graph_y, graph_w, graph_h, 12, Colors::panel_bg);
    img.stroke_round_rect_aa(graph_x, graph_y, graph_w, graph_h, 12, Colors::panel_border, 1);
    img.draw_string_aa<smallfont, true>(graph_x + 10, graph_y + 20 - smallfont::ascent, "System Performance (24h)", Colors::text_primary);
    
    for (int i = 0; i < graph_w - 20; i += 2) {
        float t = i / float(graph_w - 20);
        float value = 0.5f + 0.3f * sin(t * 8.0f) + 0.1f * sin(t * 20.0f);
        int px = graph_x + 10 + i;
        int py = graph_y + graph_h - 20 - int(value * (graph_h - 40));
        if (i > 0) {
            float prev_value = 0.5f + 0.3f * sin((i-2) / float(graph_w - 20) * 8.0f) + 0.1f * sin((i-2) / float(graph_w - 20) * 20.0f);
            int prev_py = graph_y + graph_h - 20 - int(prev_value * (graph_h - 40));
            img.draw_line_aa(px - 2, prev_py, px, py, Colors::accent_cyan, 2);
        }
    }
    
    for (int i = 0; i <= 4; i++) {
        int y = graph_y + 30 + i * (graph_h - 50) / 4;
        img.draw_line_aa(graph_x + 10, y, graph_x + graph_w - 10, y, Colors::panel_border, 1);
    }
    
    int status_x = 470, status_y = 450, status_w = 300, status_h = 150;
    img.fill_round_rect_aa(status_x, status_y, status_w, status_h, 12, Colors::panel_bg);
    img.stroke_round_rect_aa(status_x, status_y, status_w, status_h, 12, Colors::panel_border, 1);
    img.draw_string_aa<smallfont, true>(status_x + 10, status_y + 20 - smallfont::ascent, "System Status", Colors::text_primary);
    
    const char* status_items[] = {
        "Primary Power: OK",
        "Backup Power: Standby", 
        "Network: Connected",
        "Safety Systems: Armed",
        "Environmental: Normal"
    };
    
    bool status_ok[] = {true, true, true, true, true};
    
    for (int i = 0; i < 5; i++) {
        int item_y = status_y + 40 + i * 20;
        drawStatusLED(img, status_x + 15, item_y, 5, status_ok[i]);
        img.draw_string_aa<smallfont, true>(status_x + 30, item_y + 5 - smallfont::ascent, status_items[i], Colors::text_primary);
    }
    
    int btn_y = 620;
    for (int i = 0; i < 4; i++) {
        int btn_x = 50 + i * 200;
        bool pressed = (i == 1);
        uint8_t btn_color = pressed ? Colors::button_hover : Colors::button_bg;
        
        img.fill_round_rect_aa(btn_x, btn_y, 180, 40, 12, btn_color);
        img.stroke_round_rect_aa(btn_x, btn_y, 180, 40, 12, pressed ? Colors::success_green : Colors::panel_border, 2);
        
        const char* btn_labels[] = {"START", "STOP", "EMERGENCY", "RESET"};
        uint8_t text_colors[] = {Colors::text_primary, Colors::text_primary, Colors::danger_red, Colors::text_primary};
        
        img.draw_string_aa<smallfont, true>(btn_x + 20, btn_y + 25 - smallfont::ascent, btn_labels[i], text_colors[i]);
    }
    
    drawAlarmPanel(img, 50, 680, 700, 40);
    
    int indicator_x = 800;
    for (int i = 0; i < 8; i++) {
        int led_x = indicator_x + (i % 4) * 40;
        int led_y = 680 + (i / 4) * 25;
        bool state = (i % 3 != 0);
        drawStatusLED(img, led_x, led_y, 6, state);
    }
    
    img.sixel_to_cout<1>();
    
    
    return 0;
}