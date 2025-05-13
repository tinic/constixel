namespace constixel {

struct retro_pixel_petty_5h_mono {

    static constexpr const char *name = "Retro Pixel Petty 5H";
    static constexpr const char *style = "Regular";
    static constexpr int32_t size = 5;
    static constexpr int32_t ascent = 7;
    static constexpr int32_t descent = -2;
    static constexpr int32_t line_height = 9;
    static constexpr int32_t total_height = 5;

    using lookup_type = uint16_t;
    static constexpr std::array<std::pair<uint16_t, uint16_t>, 69> glyph_table{{
        { uint16_t{0x0020}, uint16_t{0x0000} },
        { uint16_t{0x0021}, uint16_t{0x0001} },
        { uint16_t{0x0022}, uint16_t{0x0002} },
        { uint16_t{0x0023}, uint16_t{0x0003} },
        { uint16_t{0x0024}, uint16_t{0x0004} },
        { uint16_t{0x0025}, uint16_t{0x0005} },
        { uint16_t{0x0026}, uint16_t{0x0006} },
        { uint16_t{0x0027}, uint16_t{0x0007} },
        { uint16_t{0x0028}, uint16_t{0x0008} },
        { uint16_t{0x0029}, uint16_t{0x0009} },
        { uint16_t{0x002a}, uint16_t{0x000a} },
        { uint16_t{0x002b}, uint16_t{0x000b} },
        { uint16_t{0x002c}, uint16_t{0x000c} },
        { uint16_t{0x002d}, uint16_t{0x000d} },
        { uint16_t{0x002e}, uint16_t{0x000e} },
        { uint16_t{0x002f}, uint16_t{0x000f} },
        { uint16_t{0x0030}, uint16_t{0x0010} },
        { uint16_t{0x0031}, uint16_t{0x0011} },
        { uint16_t{0x0032}, uint16_t{0x0012} },
        { uint16_t{0x0033}, uint16_t{0x0013} },
        { uint16_t{0x0034}, uint16_t{0x0014} },
        { uint16_t{0x0035}, uint16_t{0x0015} },
        { uint16_t{0x0036}, uint16_t{0x0016} },
        { uint16_t{0x0037}, uint16_t{0x0017} },
        { uint16_t{0x0038}, uint16_t{0x0018} },
        { uint16_t{0x0039}, uint16_t{0x0019} },
        { uint16_t{0x003a}, uint16_t{0x001a} },
        { uint16_t{0x003b}, uint16_t{0x001b} },
        { uint16_t{0x003c}, uint16_t{0x001c} },
        { uint16_t{0x003d}, uint16_t{0x001d} },
        { uint16_t{0x003e}, uint16_t{0x001e} },
        { uint16_t{0x003f}, uint16_t{0x001f} },
        { uint16_t{0x0040}, uint16_t{0x0020} },
        { uint16_t{0x005b}, uint16_t{0x0021} },
        { uint16_t{0x005c}, uint16_t{0x0022} },
        { uint16_t{0x005d}, uint16_t{0x0023} },
        { uint16_t{0x005e}, uint16_t{0x0024} },
        { uint16_t{0x005f}, uint16_t{0x0025} },
        { uint16_t{0x0060}, uint16_t{0x0026} },
        { uint16_t{0x0061}, uint16_t{0x0027} },
        { uint16_t{0x0062}, uint16_t{0x0028} },
        { uint16_t{0x0063}, uint16_t{0x0029} },
        { uint16_t{0x0064}, uint16_t{0x002a} },
        { uint16_t{0x0065}, uint16_t{0x002b} },
        { uint16_t{0x0066}, uint16_t{0x002c} },
        { uint16_t{0x0067}, uint16_t{0x002d} },
        { uint16_t{0x0068}, uint16_t{0x002e} },
        { uint16_t{0x0069}, uint16_t{0x002f} },
        { uint16_t{0x006a}, uint16_t{0x0030} },
        { uint16_t{0x006b}, uint16_t{0x0031} },
        { uint16_t{0x006c}, uint16_t{0x0032} },
        { uint16_t{0x006d}, uint16_t{0x0033} },
        { uint16_t{0x006e}, uint16_t{0x0034} },
        { uint16_t{0x006f}, uint16_t{0x0035} },
        { uint16_t{0x0070}, uint16_t{0x0036} },
        { uint16_t{0x0071}, uint16_t{0x0037} },
        { uint16_t{0x0072}, uint16_t{0x0038} },
        { uint16_t{0x0073}, uint16_t{0x0039} },
        { uint16_t{0x0074}, uint16_t{0x003a} },
        { uint16_t{0x0075}, uint16_t{0x003b} },
        { uint16_t{0x0076}, uint16_t{0x003c} },
        { uint16_t{0x0077}, uint16_t{0x003d} },
        { uint16_t{0x0078}, uint16_t{0x003e} },
        { uint16_t{0x0079}, uint16_t{0x003f} },
        { uint16_t{0x007a}, uint16_t{0x0040} },
        { uint16_t{0x007b}, uint16_t{0x0041} },
        { uint16_t{0x007c}, uint16_t{0x0042} },
        { uint16_t{0x007d}, uint16_t{0x0043} },
        { uint16_t{0x007e}, uint16_t{0x0044} }
    }};

    static constexpr hextree<hextree<0, uint16_t>::size(glyph_table), uint16_t> glyph_tree{glyph_table};

    using kerning_lookup_type = uint32_t;
    using kerning_amount_type = int32_t;
    static constexpr size_t kerning_code_shift = 16;
    static constexpr size_t kerning_amount_offset = 0x40000000;
    static constexpr hextree<0, uint32_t> kerning_tree{};

    static constexpr std::array<char_info, 69> char_table{{
        { int16_t{  17}, int16_t{  31}, int16_t{   1}, int16_t{   1}, int16_t{   5}, int16_t{   0}, int16_t{   7} },
        { int16_t{  34}, int16_t{  12}, int16_t{   1}, int16_t{   5}, int16_t{   2}, int16_t{   0}, int16_t{   2} },
        { int16_t{   4}, int16_t{  30}, int16_t{   3}, int16_t{   2}, int16_t{   4}, int16_t{   0}, int16_t{   2} },
        { int16_t{   0}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  26}, int16_t{  15}, int16_t{   3}, int16_t{   5}, int16_t{   4}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{  20}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  18}, int16_t{  30}, int16_t{   1}, int16_t{   2}, int16_t{   2}, int16_t{   0}, int16_t{   2} },
        { int16_t{  32}, int16_t{  15}, int16_t{   2}, int16_t{   5}, int16_t{   3}, int16_t{   0}, int16_t{   2} },
        { int16_t{  32}, int16_t{  20}, int16_t{   2}, int16_t{   5}, int16_t{   3}, int16_t{   0}, int16_t{   2} },
        { int16_t{   0}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  33}, int16_t{   5}, int16_t{   3}, int16_t{   3}, int16_t{   4}, int16_t{   0}, int16_t{   3} },
        { int16_t{  10}, int16_t{  30}, int16_t{   2}, int16_t{   2}, int16_t{   3}, int16_t{   0}, int16_t{   5} },
        { int16_t{  14}, int16_t{  31}, int16_t{   3}, int16_t{   1}, int16_t{   4}, int16_t{   0}, int16_t{   4} },
        { int16_t{  19}, int16_t{  30}, int16_t{   1}, int16_t{   1}, int16_t{   2}, int16_t{   0}, int16_t{   6} },
        { int16_t{  29}, int16_t{  15}, int16_t{   3}, int16_t{   5}, int16_t{   4}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{  25}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{   9}, int16_t{  20}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{   9}, int16_t{  25}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  10}, int16_t{   0}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  10}, int16_t{   5}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  10}, int16_t{  10}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  10}, int16_t{  15}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  13}, int16_t{  20}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  13}, int16_t{  25}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  14}, int16_t{   0}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  34}, int16_t{  27}, int16_t{   1}, int16_t{   4}, int16_t{   2}, int16_t{   0}, int16_t{   3} },
        { int16_t{  33}, int16_t{   8}, int16_t{   2}, int16_t{   4}, int16_t{   3}, int16_t{   0}, int16_t{   3} },
        { int16_t{  29}, int16_t{  20}, int16_t{   3}, int16_t{   5}, int16_t{   4}, int16_t{   0}, int16_t{   2} },
        { int16_t{  35}, int16_t{   0}, int16_t{   3}, int16_t{   3}, int16_t{   4}, int16_t{   0}, int16_t{   3} },
        { int16_t{  29}, int16_t{  25}, int16_t{   3}, int16_t{   5}, int16_t{   4}, int16_t{   0}, int16_t{   2} },
        { int16_t{  14}, int16_t{   5}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{   0}, int16_t{  15}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  32}, int16_t{  25}, int16_t{   2}, int16_t{   5}, int16_t{   3}, int16_t{   0}, int16_t{   2} },
        { int16_t{  30}, int16_t{   0}, int16_t{   3}, int16_t{   5}, int16_t{   4}, int16_t{   0}, int16_t{   2} },
        { int16_t{  33}, int16_t{   0}, int16_t{   2}, int16_t{   5}, int16_t{   3}, int16_t{   0}, int16_t{   2} },
        { int16_t{   7}, int16_t{  30}, int16_t{   3}, int16_t{   2}, int16_t{   4}, int16_t{   0}, int16_t{   2} },
        { int16_t{  14}, int16_t{  30}, int16_t{   4}, int16_t{   1}, int16_t{   5}, int16_t{   0}, int16_t{   6} },
        { int16_t{  12}, int16_t{  30}, int16_t{   2}, int16_t{   2}, int16_t{   3}, int16_t{   0}, int16_t{   2} },
        { int16_t{  14}, int16_t{  10}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  14}, int16_t{  15}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  17}, int16_t{  20}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  17}, int16_t{  25}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  18}, int16_t{   0}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  18}, int16_t{   5}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  18}, int16_t{  10}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  18}, int16_t{  15}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  34}, int16_t{  17}, int16_t{   1}, int16_t{   5}, int16_t{   2}, int16_t{   0}, int16_t{   2} },
        { int16_t{  21}, int16_t{  20}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  21}, int16_t{  25}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  22}, int16_t{   0}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{   0}, int16_t{  20}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  22}, int16_t{   5}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  22}, int16_t{  10}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  22}, int16_t{  15}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  25}, int16_t{  20}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  25}, int16_t{  25}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  26}, int16_t{   0}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{   0}, int16_t{  25}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  26}, int16_t{   5}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{  15}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  26}, int16_t{  10}, int16_t{   4}, int16_t{   5}, int16_t{   5}, int16_t{   0}, int16_t{   2} },
        { int16_t{  30}, int16_t{   5}, int16_t{   3}, int16_t{   5}, int16_t{   4}, int16_t{   0}, int16_t{   2} },
        { int16_t{  34}, int16_t{  22}, int16_t{   1}, int16_t{   5}, int16_t{   2}, int16_t{   0}, int16_t{   2} },
        { int16_t{  30}, int16_t{  10}, int16_t{   3}, int16_t{   5}, int16_t{   4}, int16_t{   0}, int16_t{   2} },
        { int16_t{   0}, int16_t{  30}, int16_t{   4}, int16_t{   2}, int16_t{   5}, int16_t{   0}, int16_t{   3} }
    }};

    static constexpr bool mono = true;
    static constexpr size_t glyph_bitmap_width = 38;
    static constexpr size_t glyph_bitmap_height = 32;
    static constexpr size_t glyph_bitmap_stride = 5;

    static constexpr std::array<uint8_t, 160> glyph_bitmap{{
        0x54,0x79,0xbe,0x1e,0x7c,0xfc,0x46,0x62,0x22,0x20,0x52,0x99,0xfa,0x19,0x3c,0xfa,
        0x84,0x62,0x04,0xa0,0x51,0x39,0xbf,0xf8,0xe0,0x75,0x4d,0xbe,0x65,0xa0,0xa5,0x56,
        0x63,0x65,0x70,0x75,0x64,0xba,0xe6,0x20,0x2a,0xbc,0x22,0x65,0x20,0x72,0x84,0xa2,
        0x59,0x80,0x24,0x7d,0x9d,0xbf,0x20,0xaa,0xa2,0x62,0x49,0x40,0x71,0x3b,0xee,0x50,
        0xa0,0xaa,0x86,0x66,0x61,0x20,0x24,0x7a,0x5d,0xbf,0x20,0x7c,0x5b,0xa7,0xa9,0x40,
        0xaa,0xa2,0x66,0x49,0xa0,0xb9,0x3b,0xbf,0x92,0xa0,0x81,0x26,0x66,0x24,0xa0,0x71,
        0x1b,0xa6,0x2c,0x60,0x8a,0x17,0xb8,0xb1,0xa0,0xdd,0x70,0xc0,0xca,0x60,0xaa,0x91,
        0x40,0xcc,0x60,0x8d,0x12,0x44,0xd2,0x60,0x8a,0xfa,0x3b,0x29,0xa0,0xfb,0x73,0x74,
        0xf4,0xe0,0x25,0x8c,0xcd,0x4a,0xa0,0x24,0xb3,0x4e,0x71,0xa0,0x26,0xc4,0xcd,0x4a,
        0x80,0x23,0x7b,0x74,0xcc,0xc0,0x5a,0x9b,0xf0,0x00,0x20,0xab,0x67,0xa0,0x00,0x00
        }};

};

}
