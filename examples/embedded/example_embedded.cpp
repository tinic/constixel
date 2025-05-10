#include "constixel.h"

#include <string>
#include <stdio.h>

#include "fonts/interdisplay_bold_18_mono.h"
using font = constixel::interdisplay_bold_18_mono;

//
// This an example which targets a screen like an ERM19264DNS-4
// We show the status of a nitrogen generator on this screen.
// In production this runs on a STM32G030F6 (32K Flash, 8K RAM)
//

static constixel::image<constixel::format_1bit, 192, 64, 2> screen;

bool in_fault_state = false;
float refill_ellapsed_time = 500.0;
float psi_air = 95;
float psi_nitrogen = 75;
bool solenoid_air = true;
bool solenoid_nitrogen = true;
float duty_cycle_average_air = 0.125;
float duty_cycle_average_nitrogen = 0.115;
int32_t system_time = 12356;

void draw_status() {

    screen.clear();

    static char output[64] = {};
    snprintf(output, sizeof(output), "Air:");
    screen.draw_string_mono<font>(0, -2, output, 1);
    snprintf(output, sizeof(output), "%dpsi", static_cast<int>(psi_air));
    screen.draw_string_mono<font>(192 / 2 - screen.string_width<font>(output) - 4, -2, output, 1);
    snprintf(output, sizeof(output), "N2:");
    screen.draw_string_mono<font>(192 / 2 + 4, -2, output, 1);
    snprintf(output, sizeof(output), "%dpsi", static_cast<int>(psi_nitrogen));
    screen.draw_string_mono<font>(191 - screen.string_width<font>(output), -2, output, 1);

    snprintf(output, sizeof(output), "%s", solenoid_air ? "Open" : "Clsd");
    screen.draw_string_mono<font>(0, 18, output, 1);
    snprintf(output, sizeof(output), "%s", solenoid_nitrogen ? "Open" : "Clsd");
    screen.draw_string_mono<font>(192 / 2 + 4, 18, output, 1);

    screen.draw_line(192 / 2 + 1, 0, 192 / 2 + 1, 42, 1);
    screen.draw_line(0, 42, 192, 42, 1);

    snprintf(output, sizeof(output), "%02d%%", static_cast<int>(duty_cycle_average_air * 100.0f));
    screen.draw_string_mono<font>(192 / 2 - screen.string_width<font>(output) - 4, 18, output, 1);
    snprintf(output, sizeof(output), "%02d%%", static_cast<int>(duty_cycle_average_nitrogen * 100.0f));
    screen.draw_string_mono<font>(191 - screen.string_width<font>(output), 18, output, 1);

    const int32_t h = (static_cast<int32_t>(system_time) / 3600);
    const int32_t m = (static_cast<int32_t>(system_time) / 60) % 60;
    const int32_t s = (static_cast<int32_t>(system_time)) % 60;
    snprintf(output, sizeof(output), "T%03d:%02d:%02d", static_cast<int>(h), static_cast<int>(m), static_cast<int>(s));
    screen.draw_string_mono<font>(0, 43, output, 1);

    if (in_fault_state) {
        snprintf(output, sizeof(output), "⚠Fault!");
        screen.draw_string_mono<font>(191 - screen.string_width<font>(output), 43, output, 1);
    } else {
        const int32_t em = (static_cast<int32_t>(refill_ellapsed_time) / 60) % 60;
        const int32_t es = (static_cast<int32_t>(refill_ellapsed_time)) % 60;
        snprintf(output, sizeof(output), "R%03d:%02d", static_cast<int>(em), static_cast<int>(es));
        screen.draw_string_mono<font>(190 - screen.string_width<font>(output), 43, output, 1);
    }
}

int main() {
    draw_status();

    size_t flash_memory_used = sizeof(font::glyph_bitmap) + sizeof(font::char_table) + sizeof(font::glyph_tree);
    size_t ram_memory_used = sizeof(screen);

    screen.sixel_to_cout();

    printf("This example occupies %d bytes of Flash ROM and %d bytes of RAM on a device\n", int(flash_memory_used), int(ram_memory_used));
}
