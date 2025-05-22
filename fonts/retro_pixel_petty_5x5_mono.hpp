namespace constixel {

struct retro_pixel_petty_5x5_mono {

    static constexpr const char *name = "Retro Pixel Petty 5x5";
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
        { int8_t{   9}, int8_t{  31}, int8_t{   1}, int8_t{   1}, int8_t{   6}, int8_t{   0}, int8_t{   7} },
        { int8_t{  37}, int8_t{  21}, int8_t{   1}, int8_t{   5}, int8_t{   6}, int8_t{   2}, int8_t{   2} },
        { int8_t{   0}, int8_t{  30}, int8_t{   3}, int8_t{   2}, int8_t{   6}, int8_t{   1}, int8_t{   2} },
        { int8_t{   0}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   0}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   0}, int8_t{  15}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  15}, int8_t{  30}, int8_t{   1}, int8_t{   2}, int8_t{   6}, int8_t{   2}, int8_t{   2} },
        { int8_t{  38}, int8_t{  15}, int8_t{   2}, int8_t{   5}, int8_t{   6}, int8_t{   2}, int8_t{   2} },
        { int8_t{  33}, int8_t{  18}, int8_t{   2}, int8_t{   5}, int8_t{   6}, int8_t{   1}, int8_t{   2} },
        { int8_t{   0}, int8_t{  20}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  35}, int8_t{  18}, int8_t{   3}, int8_t{   3}, int8_t{   6}, int8_t{   1}, int8_t{   3} },
        { int8_t{  11}, int8_t{  30}, int8_t{   2}, int8_t{   2}, int8_t{   6}, int8_t{   1}, int8_t{   5} },
        { int8_t{   6}, int8_t{  31}, int8_t{   3}, int8_t{   1}, int8_t{   6}, int8_t{   1}, int8_t{   4} },
        { int8_t{  10}, int8_t{  31}, int8_t{   1}, int8_t{   1}, int8_t{   6}, int8_t{   2}, int8_t{   6} },
        { int8_t{   0}, int8_t{  25}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{  15}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{  20}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{   5}, int8_t{  25}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  10}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  10}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  10}, int8_t{  15}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  34}, int8_t{  23}, int8_t{   1}, int8_t{   4}, int8_t{   6}, int8_t{   2}, int8_t{   3} },
        { int8_t{  35}, int8_t{  21}, int8_t{   2}, int8_t{   4}, int8_t{   6}, int8_t{   1}, int8_t{   3} },
        { int8_t{  45}, int8_t{   0}, int8_t{   3}, int8_t{   5}, int8_t{   6}, int8_t{   1}, int8_t{   2} },
        { int8_t{  38}, int8_t{  20}, int8_t{   3}, int8_t{   3}, int8_t{   6}, int8_t{   1}, int8_t{   3} },
        { int8_t{  45}, int8_t{   5}, int8_t{   3}, int8_t{   5}, int8_t{   6}, int8_t{   1}, int8_t{   2} },
        { int8_t{  10}, int8_t{  20}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  10}, int8_t{  25}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  45}, int8_t{  10}, int8_t{   3}, int8_t{   5}, int8_t{   6}, int8_t{   1}, int8_t{   2} },
        { int8_t{  40}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  30}, int8_t{  15}, int8_t{   3}, int8_t{   5}, int8_t{   6}, int8_t{   1}, int8_t{   2} },
        { int8_t{   3}, int8_t{  30}, int8_t{   3}, int8_t{   2}, int8_t{   6}, int8_t{   1}, int8_t{   2} },
        { int8_t{   6}, int8_t{  30}, int8_t{   5}, int8_t{   1}, int8_t{   6}, int8_t{   0}, int8_t{   6} },
        { int8_t{  13}, int8_t{  30}, int8_t{   2}, int8_t{   2}, int8_t{   6}, int8_t{   2}, int8_t{   2} },
        { int8_t{  15}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  15}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  15}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  15}, int8_t{  15}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  15}, int8_t{  20}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  15}, int8_t{  25}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  20}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  25}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  30}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  35}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  40}, int8_t{   0}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  20}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  20}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  20}, int8_t{  15}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  20}, int8_t{  20}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  20}, int8_t{  25}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  25}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  30}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  35}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  40}, int8_t{   5}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  25}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  25}, int8_t{  15}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  25}, int8_t{  20}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  25}, int8_t{  25}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  30}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  35}, int8_t{  10}, int8_t{   5}, int8_t{   5}, int8_t{   6}, int8_t{   0}, int8_t{   2} },
        { int8_t{  30}, int8_t{  20}, int8_t{   3}, int8_t{   5}, int8_t{   6}, int8_t{   1}, int8_t{   2} },
        { int8_t{  33}, int8_t{  23}, int8_t{   1}, int8_t{   5}, int8_t{   6}, int8_t{   2}, int8_t{   2} },
        { int8_t{  30}, int8_t{  25}, int8_t{   3}, int8_t{   5}, int8_t{   6}, int8_t{   1}, int8_t{   2} },
        { int8_t{  33}, int8_t{  15}, int8_t{   5}, int8_t{   3}, int8_t{   6}, int8_t{   0}, int8_t{   3} }
    }};

    static constexpr bool mono = true;
    static constexpr size_t glyph_bitmap_width = 48;
    static constexpr size_t glyph_bitmap_height = 32;
    static constexpr size_t glyph_bitmap_stride = 6;

    static constexpr std::array<uint8_t, 192> glyph_bitmap{{
        0x53,0x9c,0xe7,0xc7,0xe1,0x89,0xfc,0xe1,0x18,0x44,0x81,0x92,0x55,0x7d,0xf9,0xfc,
        0x81,0xe4,0xfe,0x63,0x18,0xc4,0x91,0x92,0x53,0x9d,0x17,0xc7,0xee,0x89,0x71,0x3f,
        0xe8,0x3b,0xcf,0xfc,0xa7,0x03,0x18,0x46,0x30,0x22,0x71,0x05,0xe8,0x57,0xce,0x21,
        0x29,0x09,0x18,0x4a,0x21,0x22,0x77,0xc9,0xef,0xb6,0x3e,0x24,0x4f,0x9c,0xf8,0xc6,
        0x3f,0x87,0x10,0x63,0x0d,0xc5,0x42,0x44,0x23,0x9d,0x0a,0xc4,0x84,0x24,0x44,0x23,
        0x08,0xc4,0x88,0x14,0x97,0xdc,0xf8,0xb8,0x9f,0x0f,0x43,0x9d,0xe8,0xc7,0xa1,0x00,
        0xa4,0x63,0x1c,0xc4,0xd6,0x00,0x69,0x9f,0x1a,0xa8,0x8a,0x00,0x94,0x43,0x19,0xa8,
        0xca,0x00,0x6b,0x9d,0xe8,0x93,0xbd,0x00,0x21,0x9d,0xf7,0x55,0xab,0x80,0xaa,0xa3,
        0x08,0xd5,0x2c,0x00,0x74,0x8d,0xe8,0xd6,0x47,0x80,0xaf,0xc1,0x08,0xa9,0x6c,0x00,
        0x20,0x89,0xf7,0x29,0xd0,0x00,0x0f,0xdf,0xff,0x47,0x44,0x00,0x14,0x2b,0x08,0xa9,
        0x60,0x00,0x27,0xaf,0xef,0x10,0xc0,0x00,0x40,0x61,0x08,0x29,0x00,0x00,0x87,0x9d,
        0x08,0x47,0x00,0x00,0xab,0xed,0x00,0x00,0x00,0x00,0xb7,0xb3,0x00,0x00,0x00,0x00
    }};

};

}
