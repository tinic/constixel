namespace constixel {

struct px437_ibm_vga_9x8_mono {

    static constexpr const char *name = "Px437 IBM VGA 9x8";
    static constexpr const char *style = "Regular";
    static constexpr int32_t size = 8;
    static constexpr int32_t ascent = 7;
    static constexpr int32_t descent = -1;
    static constexpr int32_t line_height = 8;
    static constexpr int32_t total_height = 8;

    using lookup_type = uint8_t;
    static constexpr std::array<std::pair<uint8_t, uint8_t>, 95> glyph_table{{
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
        { uint8_t{0x41}, uint8_t{0x21} },
        { uint8_t{0x42}, uint8_t{0x22} },
        { uint8_t{0x43}, uint8_t{0x23} },
        { uint8_t{0x44}, uint8_t{0x24} },
        { uint8_t{0x45}, uint8_t{0x25} },
        { uint8_t{0x46}, uint8_t{0x26} },
        { uint8_t{0x47}, uint8_t{0x27} },
        { uint8_t{0x48}, uint8_t{0x28} },
        { uint8_t{0x49}, uint8_t{0x29} },
        { uint8_t{0x4a}, uint8_t{0x2a} },
        { uint8_t{0x4b}, uint8_t{0x2b} },
        { uint8_t{0x4c}, uint8_t{0x2c} },
        { uint8_t{0x4d}, uint8_t{0x2d} },
        { uint8_t{0x4e}, uint8_t{0x2e} },
        { uint8_t{0x4f}, uint8_t{0x2f} },
        { uint8_t{0x50}, uint8_t{0x30} },
        { uint8_t{0x51}, uint8_t{0x31} },
        { uint8_t{0x52}, uint8_t{0x32} },
        { uint8_t{0x53}, uint8_t{0x33} },
        { uint8_t{0x54}, uint8_t{0x34} },
        { uint8_t{0x55}, uint8_t{0x35} },
        { uint8_t{0x56}, uint8_t{0x36} },
        { uint8_t{0x57}, uint8_t{0x37} },
        { uint8_t{0x58}, uint8_t{0x38} },
        { uint8_t{0x59}, uint8_t{0x39} },
        { uint8_t{0x5a}, uint8_t{0x3a} },
        { uint8_t{0x5b}, uint8_t{0x3b} },
        { uint8_t{0x5c}, uint8_t{0x3c} },
        { uint8_t{0x5d}, uint8_t{0x3d} },
        { uint8_t{0x5e}, uint8_t{0x3e} },
        { uint8_t{0x5f}, uint8_t{0x3f} },
        { uint8_t{0x60}, uint8_t{0x40} },
        { uint8_t{0x61}, uint8_t{0x41} },
        { uint8_t{0x62}, uint8_t{0x42} },
        { uint8_t{0x63}, uint8_t{0x43} },
        { uint8_t{0x64}, uint8_t{0x44} },
        { uint8_t{0x65}, uint8_t{0x45} },
        { uint8_t{0x66}, uint8_t{0x46} },
        { uint8_t{0x67}, uint8_t{0x47} },
        { uint8_t{0x68}, uint8_t{0x48} },
        { uint8_t{0x69}, uint8_t{0x49} },
        { uint8_t{0x6a}, uint8_t{0x4a} },
        { uint8_t{0x6b}, uint8_t{0x4b} },
        { uint8_t{0x6c}, uint8_t{0x4c} },
        { uint8_t{0x6d}, uint8_t{0x4d} },
        { uint8_t{0x6e}, uint8_t{0x4e} },
        { uint8_t{0x6f}, uint8_t{0x4f} },
        { uint8_t{0x70}, uint8_t{0x50} },
        { uint8_t{0x71}, uint8_t{0x51} },
        { uint8_t{0x72}, uint8_t{0x52} },
        { uint8_t{0x73}, uint8_t{0x53} },
        { uint8_t{0x74}, uint8_t{0x54} },
        { uint8_t{0x75}, uint8_t{0x55} },
        { uint8_t{0x76}, uint8_t{0x56} },
        { uint8_t{0x77}, uint8_t{0x57} },
        { uint8_t{0x78}, uint8_t{0x58} },
        { uint8_t{0x79}, uint8_t{0x59} },
        { uint8_t{0x7a}, uint8_t{0x5a} },
        { uint8_t{0x7b}, uint8_t{0x5b} },
        { uint8_t{0x7c}, uint8_t{0x5c} },
        { uint8_t{0x7d}, uint8_t{0x5d} },
        { uint8_t{0x7e}, uint8_t{0x5e} }
    }};

    static constexpr hextree<hextree<0, uint8_t>::size(glyph_table), uint8_t> glyph_tree{glyph_table};

    using kerning_lookup_type = uint8_t;
    using kerning_amount_type = int8_t;
    static constexpr size_t kerning_code_shift = 4;
    static constexpr size_t kerning_amount_offset = 0x40;
    static constexpr hextree<0, uint8_t> kerning_tree{};

    using char_info_type = int8_t;
    static constexpr std::array<char_info<int8_t>, 95> char_table{{
        { int8_t{   7}, int8_t{   0}, int8_t{   1}, int8_t{   1}, int8_t{   9}, int8_t{   0}, int8_t{   7} },
        { int8_t{  60}, int8_t{   7}, int8_t{   4}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{  14}, int8_t{  61}, int8_t{   5}, int8_t{   3}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{   0}, int8_t{   0}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  48}, int8_t{   7}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  14}, int8_t{  14}, int8_t{   7}, int8_t{   6}, int8_t{   9}, int8_t{   0}, int8_t{   1} },
        { int8_t{   0}, int8_t{   7}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  30}, int8_t{  48}, int8_t{   3}, int8_t{   3}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  60}, int8_t{  35}, int8_t{   4}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{  60}, int8_t{  42}, int8_t{   4}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{  46}, int8_t{  35}, int8_t{   8}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   1} },
        { int8_t{  28}, int8_t{  34}, int8_t{   6}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   1} },
        { int8_t{  30}, int8_t{  51}, int8_t{   3}, int8_t{   3}, int8_t{   9}, int8_t{   1}, int8_t{   5} },
        { int8_t{   8}, int8_t{  63}, int8_t{   6}, int8_t{   1}, int8_t{   9}, int8_t{   0}, int8_t{   3} },
        { int8_t{  24}, int8_t{  44}, int8_t{   2}, int8_t{   2}, int8_t{   9}, int8_t{   2}, int8_t{   5} },
        { int8_t{   0}, int8_t{  14}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   0}, int8_t{  21}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  41}, int8_t{  13}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  35}, int8_t{  15}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  54}, int8_t{   7}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   0}, int8_t{  28}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  47}, int8_t{  14}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  41}, int8_t{  20}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  53}, int8_t{  14}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  47}, int8_t{  21}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  53}, int8_t{  21}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  58}, int8_t{  46}, int8_t{   2}, int8_t{   6}, int8_t{   9}, int8_t{   2}, int8_t{   1} },
        { int8_t{  32}, int8_t{  20}, int8_t{   3}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   1} },
        { int8_t{  59}, int8_t{  14}, int8_t{   5}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  14}, int8_t{  52}, int8_t{   6}, int8_t{   4}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  59}, int8_t{  21}, int8_t{   5}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{  14}, int8_t{  20}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   0}, int8_t{  35}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  20}, int8_t{  20}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   0}, int8_t{  42}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   0}, int8_t{  49}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   0}, int8_t{  56}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   8}, int8_t{   0}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  15}, int8_t{   0}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  22}, int8_t{   0}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  26}, int8_t{  20}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  60}, int8_t{  49}, int8_t{   4}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{  29}, int8_t{   0}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  36}, int8_t{   0}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  43}, int8_t{   0}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  50}, int8_t{   0}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  57}, int8_t{   0}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   7}, int8_t{   7}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   7}, int8_t{  14}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  35}, int8_t{  22}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   7}, int8_t{  21}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  41}, int8_t{  27}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  47}, int8_t{  28}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  53}, int8_t{  28}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  14}, int8_t{  27}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   7}, int8_t{  28}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   7}, int8_t{  35}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  20}, int8_t{  27}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   7}, int8_t{  42}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  20}, int8_t{  44}, int8_t{   4}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{   7}, int8_t{  49}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  14}, int8_t{  44}, int8_t{   4}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{  26}, int8_t{  44}, int8_t{   7}, int8_t{   4}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{   0}, int8_t{  63}, int8_t{   8}, int8_t{   1}, int8_t{   9}, int8_t{   0}, int8_t{   7} },
        { int8_t{  20}, int8_t{  51}, int8_t{   3}, int8_t{   3}, int8_t{   9}, int8_t{   2}, int8_t{   0} },
        { int8_t{  46}, int8_t{  40}, int8_t{   7}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{   7}, int8_t{  56}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  34}, int8_t{  36}, int8_t{   6}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  14}, int8_t{   7}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  28}, int8_t{  39}, int8_t{   6}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  26}, int8_t{  27}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  21}, int8_t{  14}, int8_t{   7}, int8_t{   6}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  21}, int8_t{   7}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  33}, int8_t{  46}, int8_t{   4}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{  35}, int8_t{   7}, int8_t{   6}, int8_t{   8}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  28}, int8_t{   7}, int8_t{   7}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  26}, int8_t{  48}, int8_t{   4}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{  53}, int8_t{  41}, int8_t{   7}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  34}, int8_t{  41}, int8_t{   6}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  40}, int8_t{  41}, int8_t{   6}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  28}, int8_t{  14}, int8_t{   7}, int8_t{   6}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  41}, int8_t{   7}, int8_t{   7}, int8_t{   6}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  14}, int8_t{  34}, int8_t{   7}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  46}, int8_t{  45}, int8_t{   6}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  59}, int8_t{  28}, int8_t{   5}, int8_t{   7}, int8_t{   9}, int8_t{   1}, int8_t{   0} },
        { int8_t{  21}, int8_t{  34}, int8_t{   7}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  52}, int8_t{  46}, int8_t{   6}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  14}, int8_t{  39}, int8_t{   7}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  21}, int8_t{  39}, int8_t{   7}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  54}, int8_t{  35}, int8_t{   6}, int8_t{   6}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  14}, int8_t{  56}, int8_t{   6}, int8_t{   5}, int8_t{   9}, int8_t{   0}, int8_t{   2} },
        { int8_t{  34}, int8_t{  29}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  32}, int8_t{  27}, int8_t{   2}, int8_t{   7}, int8_t{   9}, int8_t{   3}, int8_t{   0} },
        { int8_t{  40}, int8_t{  34}, int8_t{   6}, int8_t{   7}, int8_t{   9}, int8_t{   0}, int8_t{   0} },
        { int8_t{  19}, int8_t{  61}, int8_t{   7}, int8_t{   2}, int8_t{   9}, int8_t{   0}, int8_t{   0} }
    }};

    static constexpr bool mono = true;
    static constexpr size_t glyph_bitmap_width = 64;
    static constexpr size_t glyph_bitmap_height = 64;
    static constexpr size_t glyph_bitmap_stride = 8;

    static constexpr std::array<uint8_t, 512> glyph_bitmap{{
        0x6c,0xff,0xfc,0xf0,0xfe,0x7e,0x31,0xe3,0x6c,0x62,0xc5,0x98,0x66,0x6c,0x3b,0xf3,
        0xfe,0x68,0xd3,0x00,0x66,0xcc,0x3f,0xfb,0x6c,0x78,0xf3,0x00,0x67,0x8c,0x3f,0xef,
        0xfe,0x68,0xd3,0x3e,0x66,0xcc,0x75,0xe7,0x6c,0x62,0xc1,0x9e,0x66,0x6c,0xf1,0xe3,
        0x6c,0xff,0xe0,0xfb,0xce,0x7f,0xf1,0xe3,0x38,0x70,0x77,0x0e,0x01,0xbb,0x31,0xe6,
        0x6c,0xd8,0x33,0x06,0x00,0x66,0x7f,0x3f,0x39,0x8c,0x33,0x66,0x61,0xe6,0xc0,0x3f,
        0x77,0x8d,0xf3,0xb6,0xc1,0xbe,0x78,0xe6,0xdd,0x8f,0x33,0x37,0x81,0x86,0x0c,0x36,
        0xcc,0xdb,0x33,0x36,0xd9,0x8f,0xfb,0x30,0x76,0x71,0xdf,0x3e,0x79,0x98,0x31,0xe6,
        0x07,0xfb,0x1b,0xbd,0xcf,0x39,0xff,0xe3,0x0c,0xcf,0x36,0x66,0x6f,0x19,0x86,0x66,
        0x18,0xcc,0x66,0x66,0x79,0x99,0xf0,0x6c,0x30,0xf8,0xc3,0xe7,0xc1,0x98,0x18,0xd8,
        0x60,0xc1,0x98,0x66,0x07,0x18,0x19,0x8c,0xc0,0xc3,0x1f,0xcf,0x0c,0x7f,0x99,0x86,
        0x81,0xe1,0xe3,0x33,0x79,0x9c,0xf1,0x83,0x7d,0xfb,0x37,0xb3,0x7f,0xb0,0xf3,0xd8,
        0xc6,0xcc,0x3c,0xf3,0x0f,0x61,0x9e,0x6c,0xce,0xcc,0x6c,0xff,0x19,0xfd,0x9e,0x66,
        0xde,0xf8,0xcf,0xf3,0x79,0xe6,0xf3,0xe3,0xf6,0xd8,0x0c,0xf3,0x79,0xe7,0x98,0x66,
        0xe6,0xcc,0xcc,0xf3,0xdb,0xbd,0x98,0xcc,0x7d,0xcf,0x3c,0xce,0xcf,0x3c,0xf3,0x98,
        0x1d,0x8f,0x3c,0xdb,0xc3,0xe7,0xfe,0x64,0x3d,0x8f,0x3c,0xd8,0xc7,0x71,0x6e,0x6c,
        0x6d,0x8f,0x37,0xbc,0x0c,0x38,0x66,0x7f,0xcd,0xaf,0x33,0x18,0xcc,0x0e,0x66,0x6c,
        0xff,0xfd,0xe3,0x18,0xf8,0x66,0x66,0x6c,0x0d,0xdc,0xc7,0xbc,0xcc,0x3c,0x66,0x6d,
        0x1f,0x8f,0x76,0x63,0x0c,0xe0,0xf7,0xe6,0x7d,0x8d,0xde,0x63,0x07,0x31,0x9b,0x33,
        0xc7,0x8d,0x9e,0x6f,0xde,0x30,0xf3,0x36,0xde,0xd9,0x86,0x63,0x33,0x1f,0xff,0x3c,
        0xde,0x73,0xc3,0xb3,0x30,0x30,0xf1,0xfc,0xde,0x73,0x1e,0x37,0xb3,0x31,0x98,0x3c,
        0xc0,0xdb,0x5b,0x6c,0xde,0xe1,0xe3,0xe6,0x79,0x8f,0xf9,0xcf,0xfe,0x78,0x36,0x63,
        0xfd,0xff,0xfb,0x6c,0x33,0xcd,0xf7,0xfc,0x67,0x8d,0xb6,0x37,0xb3,0xcf,0x37,0xf6,
        0x67,0x1b,0xcf,0xc4,0x33,0xcd,0xde,0xb3,0x7c,0x30,0xcc,0xce,0x33,0x79,0xf6,0x33,
        0x66,0x64,0xcc,0x1b,0x30,0x03,0x0c,0xf3,0x66,0xcc,0xcc,0x31,0x80,0x01,0xec,0xf6,
        0xfd,0xfc,0xcc,0x39,0xf0,0x00,0x3c,0xcc,0x3d,0x80,0xcc,0x19,0xb0,0x03,0xe7,0x8f,
        0x66,0xc3,0xcf,0x1b,0x30,0x00,0x03,0x36,0xc0,0x60,0x0c,0x19,0xb0,0x00,0x00,0x36,
        0xc0,0x33,0xfc,0x19,0xf8,0x00,0x00,0x06,0xc0,0x18,0x06,0x1b,0x00,0x00,0x00,0x06,
        0x66,0x0c,0x00,0x3c,0x00,0x00,0x00,0x06,0x3c,0x07,0xf0,0x00,0x00,0x00,0x00,0x0f,
        0xf9,0xc3,0xf0,0x00,0x00,0x00,0x00,0x00,0x6c,0xc2,0x60,0x00,0x00,0x00,0x00,0x00,
        0x66,0xc0,0xc0,0x00,0x00,0x00,0x00,0x00,0x66,0xf9,0x90,0x00,0x00,0x00,0x00,0x00,
        0x66,0xcf,0xf0,0x00,0x00,0x00,0x00,0x00,0x6c,0xcf,0x6e,0xc0,0x00,0x00,0x00,0x00,
        0xf9,0xbb,0x7b,0x80,0x00,0x00,0x00,0x00,0xff,0xff,0x60,0x00,0x00,0x00,0x00,0x00
    }};

};

}
