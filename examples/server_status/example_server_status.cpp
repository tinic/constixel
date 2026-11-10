#define CONSTIXEL_ENABLE_COUT

#include <algorithm>
#include <array>
#include <cstdio>
#include <numeric>

#include "../example_output.hpp"
#include "../example_style.hpp"

using namespace constixel;
using namespace example;
using Canvas = image<format_8bit, 1024, 768>;

namespace {
// A one-second histogram. Each bucket is 25 ms wide, starting at zero.
constexpr std::array<int, 8> latency = {310, 640, 240, 65, 19, 6, 3, 1};
constexpr int requests = std::accumulate(latency.begin(), latency.end(), 0);
constexpr int errors = 1;
constexpr int percentile95() {
    int count = 0;
    for (size_t i = 0; i < latency.size(); ++i) {
        count += latency[i];
        if (count * 100 >= requests * 95) {
            return (static_cast<int>(i) + 1) * 25;
        }
    }
    return static_cast<int>(latency.size()) * 25;
}

struct Service {
    const char* name;
    int latency_ms, available, desired;
};
constexpr std::array<Service, 4> services = {
    {{"api-gateway", 43, 4, 4}, {"authentication", 18, 3, 3}, {"catalog", 36, 2, 2}, {"media", 120, 1, 2}}};

struct Colors {
    uint8_t background, panel, line, ink, muted, mint, blue, amber;
    explicit Colors(const Canvas& img)
        : background(img.get_nearest_color(17, 17, 17)),
          panel(img.get_nearest_color(34, 34, 34)),
          line(img.get_nearest_color(68, 68, 68)),
          ink(img.get_nearest_color(238, 238, 238)),
          muted(img.get_nearest_color(170, 170, 170)),
          mint(img.get_nearest_color(130, 220, 190)),
          blue(img.get_nearest_color(100, 170, 230)),
          amber(img.get_nearest_color(230, 180, 80)) {
    }
};

void panel(Canvas& img, int x, int y, int w, int h, const Colors& col) {
    img.fill_round_rect_aa(x, y, w, h, 14, col.panel);
}

void metric(Canvas& img, const Colors& col, int x, const char* label, const char* value, const char* note,
            uint8_t accent) {
    panel(img, x, 108, 228, 132, col);
    text<small>(img, x + 22, 140, label, col.muted);
    text<display>(img, x + 22, 187, value, col.ink);
    img.fill_circle_aa(x + 26, 214, 3, accent);
    text<caption>(img, x + 38, 218, note, col.muted);
}
}  // namespace

int main(int argc, char** argv) {
    static Canvas img;
    const Colors col(img);
    int available = 0, desired = 0, degraded = 0;
    for (const auto& service : services) {
        available += service.available;
        desired += service.desired;
        degraded += service.available < service.desired ? 1 : 0;
    }
    char value[32], note[48];
    img.fill_rect(0, 0, 1024, 768, col.background);
    img.fill_round_rect_aa(32, 30, 6, 48, 3, col.mint);
    text<display>(img, 52, 54, "Edge / production", col.ink);
    text<small>(img, 54, 80, "Service health and request latency", col.muted);
    std::snprintf(note, sizeof(note), "%d service%s degraded", degraded, degraded == 1 ? "" : "s");
    const char* health = degraded == 0 ? "All services healthy" : note;
    const auto health_color = degraded == 0 ? col.mint : col.amber;
    img.fill_circle_aa(980 - img.string_width<small, true>(health), 45, 4, health_color);
    right_text<small>(img, 992, 50, health, health_color);
    right_text<caption>(img, 992, 77, "DEMO SNAPSHOT / 14:32:08", col.muted);

    std::snprintf(value, sizeof(value), "%d /s", requests);
    metric(img, col, 32, "REQUEST RATE", value, "Last one-second window", col.mint);
    std::snprintf(value, sizeof(value), "%d ms", percentile95());
    metric(img, col, 276, "P95 UPPER BOUND", value, "From 25 ms histogram buckets", col.blue);
    std::snprintf(value, sizeof(value), "%.2f%%", 100.0 * errors / requests);
    std::snprintf(note, sizeof(note), "%d failed of %d requests", errors, requests);
    metric(img, col, 520, "ERROR RATE", value, note, col.mint);
    std::snprintf(value, sizeof(value), "%d / %d", available, desired);
    std::snprintf(note, sizeof(note), "%d instance%s unavailable", desired - available,
                  desired - available == 1 ? "" : "s");
    metric(img, col, 764, "HEALTHY INSTANCES", value, note, health_color);

    panel(img, 32, 260, 620, 244, col);
    text<heading>(img, 56, 294, "Request latency", col.ink);
    text<small>(img, 56, 319, "Distribution / last second", col.muted);
    right_text<caption>(img, 626, 294, "REQUESTS", col.muted);
    constexpr int left = 92, top = 340, bottom = 452, step = 66;
    for (int i = 0; i <= 2; ++i) {
        const int y = bottom - i * 56;
        img.draw_line(left, y, 620, y, col.line);
        const std::array<const char*, 3> labels = {"0", "350", "700"};
        right_text<caption>(img, 78, y + 4, labels[i], col.muted);
    }
    for (size_t i = 0; i < latency.size(); ++i) {
        const int x = left + static_cast<int>(i) * step;
        const int h = std::max(2, latency[i] * (bottom - top) / 700);
        img.fill_round_rect_aa(x + 8, bottom - h, 46, h, std::min(4, h / 2), i < 4 ? col.mint : col.amber);
        std::snprintf(value, sizeof(value), "%d", latency[i]);
        centered_text<caption>(img, x + 31, bottom - h - 8, value, col.muted);
        std::snprintf(value, sizeof(value), "%d", (static_cast<int>(i) + 1) * 25);
        centered_text<caption>(img, x + 31, 476, value, col.muted);
    }
    right_text<caption>(img, 626, 495, "Upper bucket bound (ms)", col.muted);

    panel(img, 676, 260, 316, 244, col);
    text<heading>(img, 700, 294, "Resource use", col.ink);
    constexpr std::array<const char*, 3> names = {"CPU", "Memory", "Connection pool"};
    constexpr std::array<int, 3> percentages = {58, 72, 41};
    for (int i = 0; i < 3; ++i) {
        const int y = 335 + i * 62;
        text<small>(img, 700, y, names[i], col.muted);
        std::snprintf(value, sizeof(value), "%d%%", percentages[i]);
        right_text<small>(img, 968, y, value, col.ink);
        img.fill_round_rect_aa(700, y + 12, 268, 6, 3, col.line);
        const int width = 268 * std::clamp(percentages[i], 0, 100) / 100;
        img.fill_round_rect_aa(700, y + 12, width, 6, 3, col.blue);
    }

    panel(img, 32, 528, 960, 208, col);
    text<caption>(img, 56, 556, "SERVICE", col.muted);
    text<caption>(img, 390, 556, "STATUS", col.muted);
    right_text<caption>(img, 740, 556, "MEDIAN LATENCY", col.muted);
    right_text<caption>(img, 968, 556, "INSTANCES", col.muted);
    for (size_t i = 0; i < services.size(); ++i) {
        const auto& service = services[i];
        const int y = 594 + static_cast<int>(i) * 39;
        const bool healthy = service.available == service.desired;
        text<heading>(img, 56, y, service.name, col.ink);
        img.fill_circle_aa(396, y - 5, 4, healthy ? col.mint : col.amber);
        text<small>(img, 411, y, healthy ? "Healthy" : "Degraded", healthy ? col.mint : col.amber);
        std::snprintf(value, sizeof(value), "%d ms", service.latency_ms);
        right_text<small>(img, 740, y, value, col.ink);
        std::snprintf(value, sizeof(value), "%d / %d", service.available, service.desired);
        right_text<small>(img, 968, y, value, col.ink);
        if (i < services.size() - 1) {
            img.draw_line(56, y + 12, 968, y + 12, col.line);
        }
    }
    text<caption>(img, 32, 758, "CLUSTER: EU-WEST / 4 SERVICES", col.muted);
    right_text<caption>(img, 992, 758, "STATIC DEMO / NO NETWORK CONNECTION", col.muted);
    return output(img, argc, argv);
}
