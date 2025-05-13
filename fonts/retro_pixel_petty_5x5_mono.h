namespace constixel {

struct retro_pixel_petty_5x5_mono {

    static constexpr const char *name = "Retro Pixel Petty 5x5";
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
        { int16_t{   9}, int16_t{  31}, int16_t{   1}, int16_t{   1}, int16_t{   6}, int16_t{   0}, int16_t{   7} },
        { int16_t{  37}, int16_t{  21}, int16_t{   1}, int16_t{   5}, int16_t{   6}, int16_t{   2}, int16_t{   2} },
        { int16_t{   0}, int16_t{  30}, int16_t{   3}, int16_t{   2}, int16_t{   6}, int16_t{   1}, int16_t{   2} },
        { int16_t{   0}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   0}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   0}, int16_t{  15}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  15}, int16_t{  30}, int16_t{   1}, int16_t{   2}, int16_t{   6}, int16_t{   2}, int16_t{   2} },
        { int16_t{  38}, int16_t{  15}, int16_t{   2}, int16_t{   5}, int16_t{   6}, int16_t{   2}, int16_t{   2} },
        { int16_t{  33}, int16_t{  18}, int16_t{   2}, int16_t{   5}, int16_t{   6}, int16_t{   1}, int16_t{   2} },
        { int16_t{   0}, int16_t{  20}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  35}, int16_t{  18}, int16_t{   3}, int16_t{   3}, int16_t{   6}, int16_t{   1}, int16_t{   3} },
        { int16_t{  11}, int16_t{  30}, int16_t{   2}, int16_t{   2}, int16_t{   6}, int16_t{   1}, int16_t{   5} },
        { int16_t{   6}, int16_t{  31}, int16_t{   3}, int16_t{   1}, int16_t{   6}, int16_t{   1}, int16_t{   4} },
        { int16_t{  10}, int16_t{  31}, int16_t{   1}, int16_t{   1}, int16_t{   6}, int16_t{   2}, int16_t{   6} },
        { int16_t{   0}, int16_t{  25}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{  15}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{  20}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{   5}, int16_t{  25}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  10}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  10}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  10}, int16_t{  15}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  34}, int16_t{  23}, int16_t{   1}, int16_t{   4}, int16_t{   6}, int16_t{   2}, int16_t{   3} },
        { int16_t{  35}, int16_t{  21}, int16_t{   2}, int16_t{   4}, int16_t{   6}, int16_t{   1}, int16_t{   3} },
        { int16_t{  45}, int16_t{   0}, int16_t{   3}, int16_t{   5}, int16_t{   6}, int16_t{   1}, int16_t{   2} },
        { int16_t{  38}, int16_t{  20}, int16_t{   3}, int16_t{   3}, int16_t{   6}, int16_t{   1}, int16_t{   3} },
        { int16_t{  45}, int16_t{   5}, int16_t{   3}, int16_t{   5}, int16_t{   6}, int16_t{   1}, int16_t{   2} },
        { int16_t{  10}, int16_t{  20}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  10}, int16_t{  25}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  45}, int16_t{  10}, int16_t{   3}, int16_t{   5}, int16_t{   6}, int16_t{   1}, int16_t{   2} },
        { int16_t{  40}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  30}, int16_t{  15}, int16_t{   3}, int16_t{   5}, int16_t{   6}, int16_t{   1}, int16_t{   2} },
        { int16_t{   3}, int16_t{  30}, int16_t{   3}, int16_t{   2}, int16_t{   6}, int16_t{   1}, int16_t{   2} },
        { int16_t{   6}, int16_t{  30}, int16_t{   5}, int16_t{   1}, int16_t{   6}, int16_t{   0}, int16_t{   6} },
        { int16_t{  13}, int16_t{  30}, int16_t{   2}, int16_t{   2}, int16_t{   6}, int16_t{   2}, int16_t{   2} },
        { int16_t{  15}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  15}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  15}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  15}, int16_t{  15}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  15}, int16_t{  20}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  15}, int16_t{  25}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  20}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  25}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  30}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  35}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  40}, int16_t{   0}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  20}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  20}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  20}, int16_t{  15}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  20}, int16_t{  20}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  20}, int16_t{  25}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  25}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  30}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  35}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  40}, int16_t{   5}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  25}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  25}, int16_t{  15}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  25}, int16_t{  20}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  25}, int16_t{  25}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  30}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  35}, int16_t{  10}, int16_t{   5}, int16_t{   5}, int16_t{   6}, int16_t{   0}, int16_t{   2} },
        { int16_t{  30}, int16_t{  20}, int16_t{   3}, int16_t{   5}, int16_t{   6}, int16_t{   1}, int16_t{   2} },
        { int16_t{  33}, int16_t{  23}, int16_t{   1}, int16_t{   5}, int16_t{   6}, int16_t{   2}, int16_t{   2} },
        { int16_t{  30}, int16_t{  25}, int16_t{   3}, int16_t{   5}, int16_t{   6}, int16_t{   1}, int16_t{   2} },
        { int16_t{  33}, int16_t{  15}, int16_t{   5}, int16_t{   3}, int16_t{   6}, int16_t{   0}, int16_t{   3} }
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
