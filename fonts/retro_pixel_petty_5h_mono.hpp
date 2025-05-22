namespace constixel {

struct retro_pixel_petty_5h_mono {

    static constexpr const char *name = "Retro Pixel Petty 5H";
    static constexpr const char *style = "Regular";
    static constexpr int32_t size = 5;
    static constexpr int32_t ascent = 7;
    static constexpr int32_t descent = -2;
    static constexpr int32_t line_height = 9;
    static constexpr int32_t total_height = 5;

    using lookup_type = uint8_t;
    static constexpr std::array<std::pair<uint8_t, uint8_t>, 69> glyph_table{{
        { uint8_t{0x20}, uint8_t{0x00} },
        { uint8_t{0x21}, uint8_t{0x01} },
        { uint8_t{0x22}, uint8_t{0x02} },
        { uint8_t{0x23}, uint8_t{0x03} },
        { uint8_t{0x24}, uint8_t{0x04} },
        { uint8_t{0x25}, uint8_t{0x05} },
        { uint8_t{0x26}, uint8_t{0x06} },
        { uint8_t{0x27}, uint8_t{0x07} },
        { uint8_t{0x28}, uint8_t{0x08} },
        { uint8_t{0x29}, uint8_t{0x09} },
        { uint8_t{0x2a}, uint8_t{0x0a} },
        { uint8_t{0x2b}, uint8_t{0x0b} },
        { uint8_t{0x2c}, uint8_t{0x0c} },
        { uint8_t{0x2d}, uint8_t{0x0d} },
        { uint8_t{0x2e}, uint8_t{0x0e} },
        { uint8_t{0x2f}, uint8_t{0x0f} },
        { uint8_t{0x30}, uint8_t{0x10} },
        { uint8_t{0x31}, uint8_t{0x11} },
        { uint8_t{0x32}, uint8_t{0x12} },
        { uint8_t{0x33}, uint8_t{0x13} },
        { uint8_t{0x34}, uint8_t{0x14} },
        { uint8_t{0x35}, uint8_t{0x15} },
        { uint8_t{0x36}, uint8_t{0x16} },
        { uint8_t{0x37}, uint8_t{0x17} },
        { uint8_t{0x38}, uint8_t{0x18} },
        { uint8_t{0x39}, uint8_t{0x19} },
        { uint8_t{0x3a}, uint8_t{0x1a} },
        { uint8_t{0x3b}, uint8_t{0x1b} },
        { uint8_t{0x3c}, uint8_t{0x1c} },
        { uint8_t{0x3d}, uint8_t{0x1d} },
        { uint8_t{0x3e}, uint8_t{0x1e} },
        { uint8_t{0x3f}, uint8_t{0x1f} },
        { uint8_t{0x40}, uint8_t{0x20} },
        { uint8_t{0x5b}, uint8_t{0x21} },
        { uint8_t{0x5c}, uint8_t{0x22} },
        { uint8_t{0x5d}, uint8_t{0x23} },
        { uint8_t{0x5e}, uint8_t{0x24} },
        { uint8_t{0x5f}, uint8_t{0x25} },
        { uint8_t{0x60}, uint8_t{0x26} },
        { uint8_t{0x61}, uint8_t{0x27} },
        { uint8_t{0x62}, uint8_t{0x28} },
        { uint8_t{0x63}, uint8_t{0x29} },
        { uint8_t{0x64}, uint8_t{0x2a} },
        { uint8_t{0x65}, uint8_t{0x2b} },
        { uint8_t{0x66}, uint8_t{0x2c} },
        { uint8_t{0x67}, uint8_t{0x2d} },
        { uint8_t{0x68}, uint8_t{0x2e} },
        { uint8_t{0x69}, uint8_t{0x2f} },
        { uint8_t{0x6a}, uint8_t{0x30} },
        { uint8_t{0x6b}, uint8_t{0x31} },
        { uint8_t{0x6c}, uint8_t{0x32} },
        { uint8_t{0x6d}, uint8_t{0x33} },
        { uint8_t{0x6e}, uint8_t{0x34} },
        { uint8_t{0x6f}, uint8_t{0x35} },
        { uint8_t{0x70}, uint8_t{0x36} },
        { uint8_t{0x71}, uint8_t{0x37} },
        { uint8_t{0x72}, uint8_t{0x38} },
        { uint8_t{0x73}, uint8_t{0x39} },
        { uint8_t{0x74}, uint8_t{0x3a} },
        { uint8_t{0x75}, uint8_t{0x3b} },
        { uint8_t{0x76}, uint8_t{0x3c} },
        { uint8_t{0x77}, uint8_t{0x3d} },
        { uint8_t{0x78}, uint8_t{0x3e} },
        { uint8_t{0x79}, uint8_t{0x3f} },
        { uint8_t{0x7a}, uint8_t{0x40} },
        { uint8_t{0x7b}, uint8_t{0x41} },
        { uint8_t{0x7c}, uint8_t{0x42} },
        { uint8_t{0x7d}, uint8_t{0x43} },
        { uint8_t{0x7e}, uint8_t{0x44} }
    }};

    static constexpr hextree<hextree<0, uint8_t>::size(glyph_table), uint8_t> glyph_tree{glyph_table};

    using kerning_lookup_type = uint8_t;
    using kerning_amount_type = int8_t;
    static constexpr size_t kerning_code_shift = 4;
    static constexpr size_t kerning_amount_offset = 0x40;
    static constexpr hextree<0, uint8_t> kerning_tree{};

    using char_info_type = int8_t;
    static constexpr std::array<char_info<int8_t>, 69> char_table{{
        { int8_t{  17}, int8_t{  31}, int8_t{   1}, int8_t{   1}, int8_t{   5}, int8_t{   0}, int8_t{   7} },
        { int8_t{  34}, int8_t{  12}, int8_t{   1}, int8_t{   5}, int8_t{   2}, int8_t{   0}, int8_t{   2} },
        { int8_t{   4}, int8_t{  30}, int8_t{   3}, int8_t{   2}, int8_t{   4}, int8_t{   0}, int8_t{   2} },
        { int8_t{   0}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  26}, int8_t{  15}, int8_t{   3}, int8_t{   5}, int8_t{   4}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{  20}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  18}, int8_t{  30}, int8_t{   1}, int8_t{   2}, int8_t{   2}, int8_t{   0}, int8_t{   2} },
        { int8_t{  32}, int8_t{  15}, int8_t{   2}, int8_t{   5}, int8_t{   3}, int8_t{   0}, int8_t{   2} },
        { int8_t{  32}, int8_t{  20}, int8_t{   2}, int8_t{   5}, int8_t{   3}, int8_t{   0}, int8_t{   2} },
        { int8_t{   0}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  33}, int8_t{   5}, int8_t{   3}, int8_t{   3}, int8_t{   4}, int8_t{   0}, int8_t{   3} },
        { int8_t{  10}, int8_t{  30}, int8_t{   2}, int8_t{   2}, int8_t{   3}, int8_t{   0}, int8_t{   5} },
        { int8_t{  14}, int8_t{  31}, int8_t{   3}, int8_t{   1}, int8_t{   4}, int8_t{   0}, int8_t{   4} },
        { int8_t{  19}, int8_t{  30}, int8_t{   1}, int8_t{   1}, int8_t{   2}, int8_t{   0}, int8_t{   6} },
        { int8_t{  29}, int8_t{  15}, int8_t{   3}, int8_t{   5}, int8_t{   4}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{  25}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{   9}, int8_t{  20}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{   9}, int8_t{  25}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  10}, int8_t{   0}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  10}, int8_t{   5}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  10}, int8_t{  10}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  10}, int8_t{  15}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  13}, int8_t{  20}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  13}, int8_t{  25}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  14}, int8_t{   0}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  34}, int8_t{  27}, int8_t{   1}, int8_t{   4}, int8_t{   2}, int8_t{   0}, int8_t{   3} },
        { int8_t{  33}, int8_t{   8}, int8_t{   2}, int8_t{   4}, int8_t{   3}, int8_t{   0}, int8_t{   3} },
        { int8_t{  29}, int8_t{  20}, int8_t{   3}, int8_t{   5}, int8_t{   4}, int8_t{   0}, int8_t{   2} },
        { int8_t{  35}, int8_t{   0}, int8_t{   3}, int8_t{   3}, int8_t{   4}, int8_t{   0}, int8_t{   3} },
        { int8_t{  29}, int8_t{  25}, int8_t{   3}, int8_t{   5}, int8_t{   4}, int8_t{   0}, int8_t{   2} },
        { int8_t{  14}, int8_t{   5}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{   0}, int8_t{  15}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  32}, int8_t{  25}, int8_t{   2}, int8_t{   5}, int8_t{   3}, int8_t{   0}, int8_t{   2} },
        { int8_t{  30}, int8_t{   0}, int8_t{   3}, int8_t{   5}, int8_t{   4}, int8_t{   0}, int8_t{   2} },
        { int8_t{  33}, int8_t{   0}, int8_t{   2}, int8_t{   5}, int8_t{   3}, int8_t{   0}, int8_t{   2} },
        { int8_t{   7}, int8_t{  30}, int8_t{   3}, int8_t{   2}, int8_t{   4}, int8_t{   0}, int8_t{   2} },
        { int8_t{  14}, int8_t{  30}, int8_t{   4}, int8_t{   1}, int8_t{   5}, int8_t{   0}, int8_t{   6} },
        { int8_t{  12}, int8_t{  30}, int8_t{   2}, int8_t{   2}, int8_t{   3}, int8_t{   0}, int8_t{   2} },
        { int8_t{  14}, int8_t{  10}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  14}, int8_t{  15}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  17}, int8_t{  20}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  17}, int8_t{  25}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  18}, int8_t{   0}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  18}, int8_t{   5}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  18}, int8_t{  10}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  18}, int8_t{  15}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  34}, int8_t{  17}, int8_t{   1}, int8_t{   5}, int8_t{   2}, int8_t{   0}, int8_t{   2} },
        { int8_t{  21}, int8_t{  20}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  21}, int8_t{  25}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  22}, int8_t{   0}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{   0}, int8_t{  20}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  22}, int8_t{   5}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  22}, int8_t{  10}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  22}, int8_t{  15}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  25}, int8_t{  20}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  25}, int8_t{  25}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  26}, int8_t{   0}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{   0}, int8_t{  25}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  26}, int8_t{   5}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{  15}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  26}, int8_t{  10}, int8_t{   4}, int8_t{   5}, int8_t{   5}, int8_t{   0}, int8_t{   2} },
        { int8_t{  30}, int8_t{   5}, int8_t{   3}, int8_t{   5}, int8_t{   4}, int8_t{   0}, int8_t{   2} },
        { int8_t{  34}, int8_t{  22}, int8_t{   1}, int8_t{   5}, int8_t{   2}, int8_t{   0}, int8_t{   2} },
        { int8_t{  30}, int8_t{  10}, int8_t{   3}, int8_t{   5}, int8_t{   4}, int8_t{   0}, int8_t{   2} },
        { int8_t{   0}, int8_t{  30}, int8_t{   4}, int8_t{   2}, int8_t{   5}, int8_t{   0}, int8_t{   3} }
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
