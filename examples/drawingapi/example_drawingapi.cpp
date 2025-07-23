
#define CONSTIXEL_ENABLE_COUT

#include "constixel.hpp"

#include "fonts/pxplus_toshibasat_9x16_mono.hpp"
using myfont = constixel::pxplus_toshibasat_9x16_mono;

#include "fonts/ibmplexsans_bold_32_aa.hpp"
using aafont = constixel::ibmplexsans_bold_32_aa;

#include "fonts/ibmplexsans_bold_32_mono.hpp"
using monofont = constixel::ibmplexsans_bold_32_mono;

#include "fonts/ibmplexsans_bold_18_aa.hpp"
using aafontsmall = constixel::ibmplexsans_bold_18_aa;

#include "fonts/ibmplexsans_bold_18_mono.hpp"
using monofontsmall = constixel::ibmplexsans_bold_18_mono;

int main() {
    using namespace constixel;

    static image<format_8bit, 1024, 600> image;

    int32_t x_pos = 20;
    int32_t y_pos = 20;
    image.fill_rect(x_pos,y_pos, 128, 64, color::RED);
    image.draw_string_centered_mono<myfont>(x_pos+64,70+y_pos,"fill_rect",color::WHITE);

    x_pos += 200;
    image.fill_round_rect(x_pos, y_pos, 128, 64, 16, color::BLUE);
    image.draw_string_centered_mono<myfont>(x_pos+64,70+y_pos,"fill_round_rect",color::WHITE);

    x_pos += 200;
    image.fill_round_rect_aa(x_pos, y_pos, 128, 64, 16, color::GREEN);
    image.draw_string_centered_mono<myfont>(x_pos+64,70+y_pos,"fill_round_rect_aa",color::WHITE);

    x_pos += 200;
    image.fill_circle(x_pos+16, y_pos+32, 32, color::CYAN);
    image.draw_string_centered_mono<myfont>(x_pos+16,70+y_pos,"fill_circle",color::WHITE);

    x_pos += 200 - 64;
    image.fill_circle_aa(x_pos+16, y_pos+32, 32, color::MAGENTA);
    image.draw_string_centered_mono<myfont>(x_pos+16,70+y_pos,"fill_circle_aa",color::WHITE);

    x_pos = 20;
    y_pos += 100;
    image.draw_line(x_pos, y_pos, x_pos+128, y_pos+64, color::GRAY_20, 3);
    image.draw_string_centered_mono<myfont>(x_pos+64,70+y_pos,"draw_line",color::WHITE);

    x_pos += 200;
    image.draw_line_aa(x_pos, y_pos, x_pos+128, y_pos+64, color::GRAY_20);
    image.draw_string_centered_mono<myfont>(x_pos+64,70+y_pos,"draw_line_aa",color::WHITE);

    x_pos += 200;
    image.stroke_rect(x_pos, y_pos, 128, 64, color::GRAY_20, 4);
    image.draw_string_centered_mono<myfont>(x_pos+64,70+y_pos,"stroke_rect",color::WHITE);

    x_pos = 20;
    y_pos += 100;
    image.stroke_circle(x_pos + 64, y_pos + 32, 32, color::WHITE, 8);
    image.draw_string_centered_mono<myfont>(x_pos+64,70+y_pos,"stroke_circle",color::WHITE);

    x_pos += 200;
    image.stroke_circle_aa(x_pos + 64, y_pos + 32, 32, color::WHITE, 8);
    image.draw_string_centered_mono<myfont>(x_pos+64,70+y_pos,"stroke_circle_aa",color::WHITE);

    x_pos += 200;
    image.stroke_round_rect(x_pos, y_pos, 128, 64, 24, color::YELLOW, 6);
    image.draw_string_centered_mono<myfont>(x_pos+64,70+y_pos,"stroke_round_rect",color::WHITE);

    x_pos += 200;
    image.stroke_round_rect_aa(x_pos, y_pos, 128, 64, 24, color::YELLOW, 6);
    image.draw_string_centered_mono<myfont>(x_pos+64,70+y_pos,"stroke_round_rect_aa",color::WHITE);

    x_pos = 20;
    y_pos += 200;
    image.draw_string_mono<monofont, true>(x_pos, y_pos + 20, "ABCabcAVA", color::GRAY_20);
    image.draw_string_centered_mono<myfont>(x_pos+80,70+y_pos,"draw_string_mono",color::WHITE);

    x_pos += 200;
    image.draw_string_aa<aafont, true>(x_pos, y_pos + 20, "ABCabcAVA", color::GRAY_20);
    image.draw_string_centered_mono<myfont>(x_pos+80,70+y_pos,"draw_string_aa",color::WHITE);

    image.draw_string_mono<monofontsmall, true, DEGREE_90>(700, y_pos + 20, "ABCabcAVA", color::GRAY_20);
    image.draw_string_mono<monofontsmall, true, DEGREE_270>(700, y_pos + 20, "ABCabcAVA", color::GRAY_20);
    image.draw_string_mono<monofontsmall, true, DEGREE_180>(700, y_pos + 20, "ABCabcAVA", color::GRAY_20);
    image.draw_string_centered_mono<myfont>(700,130+y_pos,"draw_string_mono",color::WHITE);

    image.draw_string_aa<aafontsmall, true, DEGREE_90>(900, y_pos + 20, "ABCabcAVA", color::GRAY_20);
    image.draw_string_aa<aafontsmall, true, DEGREE_270>(900, y_pos + 20, "ABCabcAVA", color::GRAY_20);
    image.draw_string_aa<aafontsmall, true, DEGREE_180>(900, y_pos + 20, "ABCabcAVA", color::GRAY_20);
    image.draw_string_centered_mono<myfont>(900,130+y_pos,"draw_string_aa",color::WHITE);
    
    image.sixel_to_cout<
#ifdef __APPLE__
    2
#endif  // #ifdef __APPLE__
    >();
}
