#ifndef _SIXEL_TOOLS_H_
#define _SIXEL_TOOLS_H_
#include <cstdint>
#include <array>

namespace sixel {

    template<size_t W, size_t H> class format_1bit {
    public:
        static constexpr size_t bits_per_pixel = 1;
        static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
        static constexpr size_t image_size = H * bytes_per_line;
        static constexpr std::array<uint8_t, 3 * (1UL << bits_per_pixel)> palette = { 
            0x00, 0x00, 0x00,
            0xff, 0xff, 0xff};
        static constexpr void plot(std::array<uint8_t, image_size> &data, size_t x0, size_t y, uint32_t col) {
            col %= 1UL<<bits_per_pixel;
            size_t x8 = x0 / 8; x0 %= 8;
            uint8_t *yptr = &data.data()[y * bytes_per_line];
            yptr[x8] &= ~static_cast<uint8_t>(1UL << (7-x0));
            yptr[x8] |=  static_cast<uint8_t>(col << (7-x0));
        }
        static constexpr void span(std::array<uint8_t, image_size> &data, size_t xl0, size_t xr0, size_t y, uint32_t col) {
            col %= 1UL<<bits_per_pixel;
            size_t xl8 = xl0 / 8; xl0 %= 8;
            size_t xr8 = xr0 / 8; xr0 %= 8;
            size_t xs8 = xr8 - xl8;
            uint8_t c8 = col << 7 | col << 6 | col << 5 | col << 4 | col << 3 | col << 2 | col << 1 | col << 0;
            constexpr uint8_t ml[] = { 0b11111111, 0b01111111, 0b00111111, 0b00011111, 0b00001111, 0b00000111, 0b00000011, 0b00000001 };
            constexpr uint8_t mr[] = { 0b00000000, 0b10000000, 0b11000000, 0b11100000, 0b11110000, 0b11111000, 0b11111100, 0b11111110 };
            uint8_t *yptr = &data.data()[y * bytes_per_line];
            if (xs8 > 0) {
                yptr[xl8] &= ~ml[xl0]; 
                yptr[xl8] |=  ml[xl0] & c8; 
                for (size_t x = xl8+1; x < xr8; x++) {
                    yptr[+ x] = c8;
                }
                yptr[xr8] &= ~mr[xr0]; 
                yptr[xr8] |=  mr[xr0] & c8; 
            } else {
                yptr[xl8] &= ~(ml[xl0] & mr[xr0]); 
                yptr[xl8] |=  (ml[xl0] & mr[xr0] & c8); 
            }
        }
    };

    template<size_t W, size_t H> class format_2bit {
    public:
        static constexpr size_t bits_per_pixel = 2;
        static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
        static constexpr size_t image_size = H * bytes_per_line;
        static constexpr std::array<uint8_t, 3 * (1UL << bits_per_pixel)> palette = { 
            0x00, 0x00, 0x00,
            0xff, 0xff, 0xff,
            0xff, 0x00, 0x00,
            0x00, 0xff, 0x0};
        static constexpr void plot(std::array<uint8_t, image_size> &data, size_t x0, size_t y, uint32_t col) {
            col %= 1UL<<bits_per_pixel;
            size_t x4 = x0 / 4; x0 %= 4;
            uint8_t *yptr = &data.data()[y * bytes_per_line];
            yptr[x4] &= ~static_cast<uint8_t>(3UL << (6-x0*2));
            yptr[x4] |=  static_cast<uint8_t>(col << (6-x0*2));
        }
        static constexpr void span(std::array<uint8_t, image_size> &data, size_t xl0, size_t xr0, size_t y, uint32_t col) {
            col %= 1UL<<bits_per_pixel;
            size_t xl4 = xl0 / 4; xl0 %= 4;
            size_t xr4 = xr0 / 4; xr0 %= 4;
            size_t xs4 = xr4 - xl4;
            uint8_t c4 = col << 6 | col << 4 | col << 2 | col << 0;
            constexpr uint8_t ml[] = { 0b11111111, 0b00111111, 0b00001111, 0b00000011 };
            constexpr uint8_t mr[] = { 0b00000000, 0b11000000, 0b11110000, 0b11111100 };
            uint8_t *yptr = &data.data()[y * bytes_per_line];
            if (xs4 > 0) {
                yptr[xl4] &= ~ml[xl0]; 
                yptr[xl4] |=  ml[xl0] & c4; 
                for (size_t x = xl4+1; x < xr4; x++) {
                    yptr[+ x] = c4;
                }
                yptr[xr4] &= ~mr[xr0]; 
                yptr[xr4] |=  mr[xr0] & c4; 
            } else {
                yptr[xl4] &= ~(ml[xl0] & mr[xr0]); 
                yptr[xl4] |=  (ml[xl0] & mr[xr0] & c4); 
            }
        }
    };

    template<size_t W, size_t H> class format_4bit {
    public:
        static constexpr size_t bits_per_pixel = 3;
        static constexpr size_t bytes_per_line = (W * bits_per_pixel + 7) / 8;
        static constexpr size_t image_size = H * bytes_per_line;
        static constexpr std::array<uint8_t, 3 * (1UL << bits_per_pixel)> palette = { 
            0x00, 0x00, 0x00,
            0xff, 0xff, 0xff,
            0xff, 0x00, 0x00,
            0x00, 0xff, 0x00,
            0x00, 0x00, 0xff,
            0xff, 0xff, 0x00,
            0x00, 0xff, 0xff,
            0xff, 0x00, 0xff};
        static constexpr void plot(std::array<uint8_t, image_size> &data, size_t x0, size_t y, uint32_t col) {
            col %= 1UL<<bits_per_pixel;
            size_t x2 = x0 / 2; x0 %= 2;
            uint8_t *yptr = &data.data()[y * bytes_per_line];
            yptr[x2] &= ~static_cast<uint8_t>(7UL << (2-x0*2));
            yptr[x2] |=  static_cast<uint8_t>(col << (2-x0*2));
        }
        static constexpr void span(std::array<uint8_t, image_size> &data, size_t xl0, size_t xr0, size_t y, uint32_t col) {
            col %= 1UL<<bits_per_pixel;
            size_t xl2 = xl0 / 2; xl0 %= 2;
            size_t xr2 = xr0 / 2; xr0 %= 2;
            size_t xs2 = xr2 - xl2;
            uint8_t c2 = col << 4 | col << 0;
            constexpr uint8_t ml[] = { 0b11111111, 0b00001111 };
            constexpr uint8_t mr[] = { 0b00000000, 0b11110000 };
            uint8_t *yptr = &data.data()[y * bytes_per_line];
            if (xs2 > 0) {
                yptr[xl2] &= ~ml[xl0]; 
                yptr[xl2] |=  ml[xl0] & c2; 
                for (size_t x = xl2+1; x < xr2; x++) {
                    yptr[+ x] = c2;
                }
                yptr[xr2] &= ~mr[xr0]; 
                yptr[xr2] |=  mr[xr0] & c2; 
            } else {
                yptr[xl2] &= ~(ml[xl0] & mr[xr0]); 
                yptr[xl2] |=  (ml[xl0] & mr[xr0] & c2); 
            }
        }
    };

    template<template<size_t, size_t> class T, size_t W, size_t H> class image {
        static_assert(sizeof(W)>=sizeof(uint32_t));
        static_assert(sizeof(H)>=sizeof(uint32_t));
        static_assert(W<=16384 && H<=16384);
    public:
        void plot(int32_t x, int32_t y, uint32_t col) {
            size_t _x = static_cast<size_t>(x); _x %= W;
            size_t _y = static_cast<size_t>(y); _y %= H;
            T<W, H>::plot(data, _x, _y, col);
        }
        void span(int32_t xl,  int32_t xr, int32_t y, uint32_t col) {
            size_t _xl = static_cast<size_t>(xl); _xl %= W;
            size_t _xr = static_cast<size_t>(xr); _xr %= W;
            size_t _y  = static_cast<size_t>(y ); _y  %= H;
            T<W, H>::span(data, _xl, _xr, _y, col);
        }
    private:
        std::array<uint8_t, T<W, H>::image_size> data;
        T<W, H> format;
    };

}
#endif  // #ifndef _SIXEL_TOOLS_H_
