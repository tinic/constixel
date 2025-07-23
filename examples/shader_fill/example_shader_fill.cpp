#define CONSTIXEL_ENABLE_COUT

#include "constixel.hpp"

#include "fonts/ibmplexsans_semibold_18_mono.hpp"
using font = constixel::ibmplexsans_semibold_18_mono;

#include <array>
#include <cmath>
#include <fstream>
#include <memory>

int main() {
    using namespace constixel;
    
    static auto img = std::make_unique<image<format_32bit, 1024, 768>>();
    
    // Clear background
    img->fill_rect(0, 0, static_cast<int32_t>(img->width()), static_cast<int32_t>(img->height()), color::BLACK);
    
    // Title
    img->draw_string_mono<font>(20, 20, "Constixel Shader Fill Demo - All Shapes", color::WHITE);
    
    // === RECTANGLES ===
    img->draw_string_mono<font>(20, 50, "Rectangles:", color::CYAN);
    
    // Linear gradient rectangle
    img->fill_rect(50, 70, 150, 80, [](float u, float /* v */, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        return {1.0f - u, 0.0f, u, 1.0f}; // Red to blue horizontally
    });
    img->draw_string_mono<font>(50, 155, "Linear Gradient", color::WHITE);
    
    // Radial gradient rectangle
    img->fill_rect(220, 70, 150, 80, [](float u, float v, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        float cx = 0.5f, cy = 0.5f;
        float dist = std::sqrt((u - cx) * (u - cx) + (v - cy) * (v - cy));
        float intensity = std::max(0.0f, 1.0f - dist * 2.0f);
        return {intensity, intensity * 0.8f, 1.0f, 1.0f};
    });
    img->draw_string_mono<font>(220, 155, "Radial Gradient", color::WHITE);
    
    // Checkerboard rectangle
    img->fill_rect(390, 70, 150, 80, [](float u, float v, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        int checker = (int(u * 8) + int(v * 8)) % 2;
        float gray = checker ? 1.0f : 0.2f;
        return {gray, gray, gray, 1.0f};
    });
    img->draw_string_mono<font>(390, 155, "Checkerboard", color::WHITE);
    
    // === CIRCLES ===
    img->draw_string_mono<font>(20, 185, "Circles:", color::CYAN);
    
    // Radial gradient circle
    img->fill_circle_aa(125, 240, 40, [](float u, float v, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        float cx = 0.5f, cy = 0.5f;
        float dist = std::sqrt((u - cx) * (u - cx) + (v - cy) * (v - cy));
        float intensity = std::max(0.0f, 1.0f - dist * 2.0f);
        return {intensity, intensity * 0.8f, 1.0f, 1.0f}; // Blue center fading out
    });
    img->draw_string_mono<font>(85, 290, "Radial Circle", color::WHITE);
    
    // Pattern circle
    img->fill_circle_aa(295, 240, 40, [](float u, float v, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        float rings = std::sin((u + v) * 10.0f) * 0.5f + 0.5f;
        return {rings, 1.0f - rings, rings * 0.5f, 1.0f};
    });
    img->draw_string_mono<font>(250, 290, "Pattern Circle", color::WHITE);
    
    // Spiral circle
    img->fill_circle_aa(465, 240, 40, [](float u, float v, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        float cx = 0.5f, cy = 0.5f;
        float dx = u - cx, dy = v - cy;
        float angle = std::atan2(dy, dx);
        float dist = std::sqrt(dx * dx + dy * dy);
        float spiral = std::sin(angle * 3.0f + dist * 15.0f) * 0.5f + 0.5f;
        return {spiral, spiral * 0.7f, 1.0f - spiral, 1.0f};
    });
    img->draw_string_mono<font>(425, 290, "Spiral Circle", color::WHITE);
    
    // === ROUNDED RECTANGLES ===
    img->draw_string_mono<font>(20, 320, "Rounded Rectangles:", color::CYAN);
    
    // Gradient rounded rectangle
    img->fill_round_rect_aa(50, 350, 150, 80, 20, [](float u, float v, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        return {u, v, 1.0f - u, 1.0f}; // Corner-to-corner gradient
    });
    img->draw_string_mono<font>(50, 440, "Gradient RoundRect", color::WHITE);
    
    // Wave pattern rounded rectangle
    img->fill_round_rect_aa(220, 350, 150, 80, 20, [](float u, float v, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        float wave = (std::sin(u * 8.0f) + std::sin(v * 8.0f)) * 0.5f + 0.5f;
        return {wave, wave * 0.5f, 1.0f - wave, 1.0f};
    });
    img->draw_string_mono<font>(220, 440, "Wave RoundRect", color::WHITE);
    
    // Complex pattern rounded rectangle
    img->fill_round_rect_aa(390, 350, 150, 80, 20, [](float u, float v, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        float x = (u - 0.5f) * 4.0f;
        float y = (v - 0.5f) * 4.0f;
        float r = std::sqrt(x * x + y * y);
        float intensity = std::sin(r * 2.0f) * std::exp(-r * 0.5f) * 0.5f + 0.5f;
        return {intensity, intensity * 0.7f, 1.0f - intensity, 1.0f};
    });
    img->draw_string_mono<font>(390, 440, "Math RoundRect", color::WHITE);
    
    // === LARGE SHOWCASE ===
    img->draw_string_mono<font>(20, 470, "Complex Showcase:", color::CYAN);
    
    // Large complex mathematical pattern
    img->fill_rect(50, 500, 350, 120, [](float u, float v, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        float x = (u - 0.5f) * 8.0f;
        float y = (v - 0.5f) * 8.0f;
        float r = std::sqrt(x * x + y * y);
        
        // Multiple overlapping effects
        float ripple = std::sin(r * 1.5f - 8.0f) * std::exp(-r * 0.3f);
        float spiral = std::sin(std::atan2(y, x) * 4.0f + r * 2.0f);
        float waves = std::sin(u * 12.0f) * std::sin(v * 12.0f);
        
        float red = (ripple + 1.0f) * 0.5f;
        float green = (spiral + 1.0f) * 0.3f + 0.2f;
        float blue = (waves + 1.0f) * 0.4f + 0.3f;
        
        return {red, green, blue, 1.0f};
    });
    img->draw_string_mono<font>(50, 630, "Combined: Ripple + Spiral + Waves", color::WHITE);
    
    // Large showcase circle
    img->fill_circle_aa(500, 560, 60, [](float u, float v, int32_t /* px */, int32_t /* py */) -> std::array<float, 4> {
        float cx = 0.5f, cy = 0.5f;
        float dx = u - cx, dy = v - cy;
        float dist = std::sqrt(dx * dx + dy * dy);
        float angle = std::atan2(dy, dx);
        
        // Mandala-like pattern
        float petals = std::sin(angle * 6.0f) * 0.5f + 0.5f;
        float rings = std::sin(dist * 12.0f) * 0.5f + 0.5f;
        float fade = std::max(0.0f, 1.0f - dist * 1.5f);
        
        return {petals * fade, rings * fade, fade, 1.0f};
    });
    img->draw_string_mono<font>(440, 630, "Mandala Circle", color::WHITE);
    
#ifdef _MSC_VER
    img->sixel_to_cout<1>();
#else  // #ifdef _MSC_VER
    img->png_to_iterm();
#endif  // #ifdef _MSC_VER
    
    return 0;
}