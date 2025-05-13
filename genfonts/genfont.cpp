
#include <charconv>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "../constixel.h"
#include "fontbm/src/external/cxxopts.hpp"
#include "fontbm/src/external/lodepng/lodepng.h"

struct Char {
    uint32_t id = 0;
    uint16_t x = 0, y = 0, width = 0, height = 0;
    int16_t xoffset = 0, yoffset = 0, xadvance = 0;
    uint16_t page = 0;
    uint8_t chnl = 0;
};

struct Kerning {
    uint32_t first = 0;
    uint32_t second = 0;
    int32_t amount = 0;
};

struct Font {
    std::string face{};
    std::string style{};
    int size = 0;
    int lineHeight = 0;
    int base = 0;
    int descent = 0;
    int scaleW = 0, scaleH = 0;
    int pages = 0;
    std::string page_file{};
    std::vector<Char> chars{};
    std::vector<Kerning> kernings{};
    int smooth = 0;
    int totalHeight = 0;
};

static std::string_view trim(std::string_view sv) {
    const auto begin = sv.find_first_not_of(' ');
    const auto end = sv.find_last_not_of(' ');
    return begin == sv.npos ? "" : sv.substr(begin, end - begin + 1);
}

static std::unordered_map<std::string_view, std::string_view> split_kv(std::string_view line) {
    std::unordered_map<std::string_view, std::string_view> kv;
    for (std::size_t i = 0; i < line.size();) {
        while (i < line.size() && line[i] == ' ') ++i;  // skip blanks

        const std::size_t key_beg = i;
        while (i < line.size() && line[i] != '=') ++i;
        if (i == line.size())
            break;
        const std::size_t key_len = i++ - key_beg;  // past '='

        std::string_view val;
        if (i < line.size() && line[i] == '"') {  // quoted value
            const std::size_t val_beg = ++i;
            while (i < line.size() && line[i] != '"') ++i;
            val = line.substr(val_beg, i - val_beg);
            if (i < line.size())
                ++i;  // skip closing quote
        } else {      // unquoted value
            const std::size_t val_beg = i;
            while (i < line.size() && line[i] != ' ') ++i;
            val = line.substr(val_beg, i - val_beg);
        }

        kv.emplace(line.substr(key_beg, key_len), val);
    }
    return kv;
}

static int to_int(std::string_view sv) {
    int v = 0;
    std::from_chars(sv.data(), sv.data() + sv.size(), v);
    return v;
}

static Font parse_fnt(const std::filesystem::path& path) {
    std::ifstream in(path);
    Font font{};
    std::string line;
    while (std::getline(in, line)) {
        const auto sv = trim(line);
        const auto sp = sv.find(' ');
        const auto cmd = sv.substr(0, sp);
        const auto args = sp == sv.npos ? "" : sv.substr(sp + 1);
        const auto kv = split_kv(args);

        if (cmd == "info") {
            font.face = std::string(kv.at("face"));
            font.size = to_int(kv.at("size"));
            font.smooth = to_int(kv.at("smooth"));
            font.style = std::string(kv.at("style"));
        } else if (cmd == "common") {
            font.lineHeight = to_int(kv.at("lineHeight"));
            font.base = to_int(kv.at("base"));
            font.scaleW = to_int(kv.at("scaleW"));
            font.scaleH = to_int(kv.at("scaleH"));
            font.pages = to_int(kv.at("pages"));
            font.totalHeight = to_int(kv.at("totalHeight"));
            font.descent = to_int(kv.at("descent"));
        } else if (cmd == "page") {
            font.page_file = std::string(kv.at("file"));
        } else if (cmd == "char") {
            Char c{};
            c.id = to_int(kv.at("id"));
            c.x = to_int(kv.at("x"));
            c.y = to_int(kv.at("y"));
            c.width = to_int(kv.at("width"));
            c.height = to_int(kv.at("height"));
            c.xoffset = to_int(kv.at("xoffset"));
            c.yoffset = to_int(kv.at("yoffset"));
            c.xadvance = to_int(kv.at("xadvance"));
            c.page = to_int(kv.at("page"));
            c.chnl = to_int(kv.at("chnl"));
            font.chars.push_back(c);
        } else if (cmd == "kerning") {
            Kerning k{};
            k.first = to_int(kv.at("first"));
            k.second = to_int(kv.at("second"));
            k.amount = to_int(kv.at("amount"));
            font.kernings.push_back(k);
        }
    }
    return font;
}

struct Config {
    std::string font_file;
    std::string out_dir;
} config;

static void parseCommandLine(int argc, char* argv[]) {
    try {
        cxxopts::Options options("genfont", "Command line font generator for constixel, compatible with bmfont");
        options.add_options()("font-file", "path to ttf/otf file, required", cxxopts::value<std::string>(config.font_file))(
            "out-dir", "path to output path, required", cxxopts::value<std::string>(config.out_dir));
        auto result = options.parse(argc, argv);

        if (!result.count("font-file")) {
            std::cout << options.help() << std::endl;
            throw std::runtime_error("--font-file required");
        }

        if (!result.count("out-dir")) {
            std::cout << options.help() << std::endl;
            throw std::runtime_error("--out-dir required");
        }
    } catch (const cxxopts::OptionException& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        throw std::exception();
    }
}

std::string slugify(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return c == '-' ? '_' : (c == ' ' ? '_' : std::tolower(c));
    });
    return s;
}

int main(int argc, char* argv[]) {
    parseCommandLine(argc, argv);
    Font font = parse_fnt(config.font_file);
    std::print("Converting font: {}\n", config.font_file);
    if (font.page_file.size() > 0) {
        if (font.face.size() == 0) {
            throw std::runtime_error("Font has no name.");
        }
        std::vector<uint8_t> rgbaimage;
        uint32_t w = 0;
        uint32_t h = 0;
        auto path = std::filesystem::path{config.font_file}.parent_path() / font.page_file;
        if (lodepng::decode(rgbaimage, w, h, path) == 0) {
            bool is_mono = font.smooth == 0;

            std::stringstream ss;

            std::string name = std::format("{}_{}", slugify(std::filesystem::path(config.font_file).stem().string()), is_mono ? "mono" : "aa");

            ss << std::format("namespace constixel {{\n\n", slugify(font.face), abs(font.size));
            ss << std::format("struct {} {{\n\n", name);

            ss << std::format("    static constexpr const char *name = \"{}\";\n", font.face);
            ss << std::format("    static constexpr const char *style = \"{}\";\n", font.style);
            ss << std::format("    static constexpr int32_t size = {};\n", std::abs(font.size));
            ss << std::format("    static constexpr int32_t ascent = {};\n", font.base);
            ss << std::format("    static constexpr int32_t descent = {};\n", font.descent);
            ss << std::format("    static constexpr int32_t line_height = {};\n", font.lineHeight);
            ss << std::format("    static constexpr int32_t total_height = {};\n\n", font.totalHeight);

            uint32_t max_id = 0;
            for (auto c : font.chars) {
                max_id = std::max(c.id, max_id);
            }
            if (font.chars.size() < 0xFF && max_id < 0xFF) {
                ss << std::format("    using lookup_type = uint8_t;\n", name);
                ss << std::format("    static constexpr std::array<std::pair<uint8_t, uint8_t>, {}> glyph_table{{{{\n", font.chars.size());
                for (size_t c = 0; c < font.chars.size(); c++) {
                    ss << std::format("        {{ uint8_t{{0x{:02x}}}, uint8_t{{0x{:02x}}} }}{}\n", font.chars[c].id, c,
                                      (c < font.chars.size() - 1) ? "," : "");
                }
            } else if (font.chars.size() < 0xFFFF && max_id < 0xFFFF) {
                ss << std::format("    using lookup_type = uint16_t;\n", name);
                ss << std::format("    static constexpr std::array<std::pair<uint16_t, uint16_t>, {}> glyph_table{{{{\n", font.chars.size());
                for (size_t c = 0; c < font.chars.size(); c++) {
                    ss << std::format("        {{ uint16_t{{0x{:04x}}}, uint16_t{{0x{:04x}}} }}{}\n", font.chars[c].id, c,
                                      (c < font.chars.size() - 1) ? "," : "");
                }
            } else {
                ss << std::format("    using lookup_type = uint32_t;\n", name);
                ss << std::format("    static constexpr std::array<std::pair<uint32_t, uint32_t>, {}> glyph_table{{{{\n", font.chars.size());
                for (size_t c = 0; c < font.chars.size(); c++) {
                    ss << std::format("        {{ uint32_t{{0x{:08x}}}, uint32_t{{0x{:08x}}} }}{}\n", font.chars[c].id, c,
                                      (c < font.chars.size() - 1) ? "," : "");
                }
            }
            ss << std::format("    }}}};\n\n");
            if (font.chars.size() < 0xFF && max_id < 0xFF) {
                ss << std::format("    static constexpr hextree<hextree<0, uint8_t>::size(glyph_table), uint8_t> glyph_tree{{glyph_table}};\n\n",
                                  font.chars.size());
            } else if (font.chars.size() < 0xFFFF && max_id < 0xFFFF) {
                ss << std::format("    static constexpr hextree<hextree<0, uint16_t>::size(glyph_table), uint16_t> glyph_tree{{glyph_table}};\n\n",
                                  font.chars.size());
            } else {
                ss << std::format("    static constexpr hextree<hextree<0, uint32_t>::size(glyph_table), uint32_t> glyph_tree{{glyph_table}};\n\n",
                                  font.chars.size());
            }

            uint32_t max_utf32 = 0;
            for (auto c : font.kernings) {
                max_utf32 = std::max(c.first, max_utf32);
                max_utf32 = std::max(c.second, max_utf32);
            }

            if (font.kernings.size() > 0 && font.kernings.size() < 0xFFFF && max_utf32 < 0xFF) {
                ss << std::format("    using kerning_lookup_type = uint16_t;\n", name);
                ss << std::format("    using kerning_amount_type = int16_t;\n", name);
                ss << std::format("    static constexpr size_t kerning_code_shift = 8;\n", name);
                ss << std::format("    static constexpr int16_t kerning_amount_offset = 0x4000;\n", name);
                ss << std::format("    static constexpr std::array<std::pair<uint16_t, uint16_t>, {}> kerning_table{{{{\n", font.kernings.size());
                for (size_t c = 0; c < font.kernings.size(); c++) {
                    ss << std::format("        {{ uint16_t{{0x{:04x}}}, uint16_t{{0x{:04x}}} }}{}\n", font.kernings[c].first << 8 | font.kernings[c].second,
                                      static_cast<uint16_t>(font.kernings[c].amount + 0x4000), (c < font.kernings.size() - 1) ? "," : "");
                }
                ss << std::format("    }}}};\n\n");

                ss << std::format("    static constexpr hextree<hextree<0, uint16_t>::size(kerning_table), uint16_t> kerning_tree{{kerning_table}};\n\n");

            } else if (font.kernings.size() > 0 && max_utf32 < 65535) {
                ss << std::format("    using kerning_lookup_type = uint32_t;\n", name);
                ss << std::format("    using kerning_amount_type = int32_t;\n", name);
                ss << std::format("    static constexpr size_t kerning_code_shift = 16;\n", name);
                ss << std::format("    static constexpr int32_t kerning_amount_offset = 0x40000000;\n", name);
                ss << std::format("    static constexpr std::array<std::pair<uint32_t, uint32_t>, {}> kerning_table{{{{\n", font.kernings.size());
                for (size_t c = 0; c < font.kernings.size(); c++) {
                    ss << std::format("        {{ uint32_t{{0x{:08x}}}, uint32_t{{0x{:08x}}} }}{}\n", font.kernings[c].first << 16 | font.kernings[c].second,
                                      static_cast<uint16_t>(font.kernings[c].amount + 0x4000), (c < font.kernings.size() - 1) ? "," : "");
                }
                ss << std::format("    }}}};\n\n");

                ss << std::format("    static constexpr hextree<hextree<0, uint32_t>::size(kerning_table), uint32_t> kerning_tree{{kerning_table}};\n\n");
            } else {
                ss << std::format("    using kerning_lookup_type = uint8_t;\n", name);
                ss << std::format("    using kerning_amount_type = int8_t;\n", name);
                ss << std::format("    static constexpr size_t kerning_code_shift = 4;\n", name);
                ss << std::format("    static constexpr size_t kerning_amount_offset = 0x40;\n", name);

                ss << std::format("    static constexpr hextree<0, uint8_t> kerning_tree{{}};\n\n");
            }

            int32_t max_value = 0;
            for (auto c : font.chars) {
                max_value = std::max(std::abs(static_cast<int32_t>(c.x)), max_value);
                max_value = std::max(std::abs(static_cast<int32_t>(c.y)), max_value);
                max_value = std::max(std::abs(static_cast<int32_t>(c.width)), max_value);
                max_value = std::max(std::abs(static_cast<int32_t>(c.height)), max_value);
                max_value = std::max(std::abs(static_cast<int32_t>(c.xadvance)), max_value);
                max_value = std::max(std::abs(static_cast<int32_t>(c.xoffset)), max_value);
                max_value = std::max(std::abs(static_cast<int32_t>(c.yoffset)), max_value);
            }
            const char* chars_unit = max_value > 127 ? (max_value > 32767 ? "int32_t" : "int16_t") : "int8_t";

            ss << std::format("    using char_info_type = {};\n", chars_unit);
            ss << std::format("    static constexpr std::array<char_info<{}>, {}> char_table{{{{\n", chars_unit, font.chars.size());
            for (size_t c = 0; c < font.chars.size(); c++) {
                ss << std::format("        {{ {}{{{:4}}}, {}{{{:4}}}, {}{{{:4}}}, {}{{{:4}}}, {}{{{:4}}}, {}{{{:4}}}, {}{{{:4}}} }}{}\n", chars_unit,
                                  font.chars[c].x, chars_unit, font.chars[c].y, chars_unit, font.chars[c].width, chars_unit, font.chars[c].height, chars_unit,
                                  font.chars[c].xadvance, chars_unit, font.chars[c].xoffset, chars_unit, font.chars[c].yoffset,
                                  (c < font.chars.size() - 1) ? "," : "");
            }

            ss << std::format("    }}}};\n\n");

            if (is_mono) {
                size_t bpr = ((w + 7) / 8);
                std::vector<uint8_t> bitmap;
                bitmap.assign(h * bpr, 0);

                ss << std::format("    static constexpr bool mono = true;\n");
                ss << std::format("    static constexpr size_t glyph_bitmap_width = {};\n", w);
                ss << std::format("    static constexpr size_t glyph_bitmap_height = {};\n", h);
                ss << std::format("    static constexpr size_t glyph_bitmap_stride = {};\n\n", bpr);

                ss << std::format("    static constexpr std::array<uint8_t, {}> glyph_bitmap{{{{\n", h * bpr);

                auto plot = [&bitmap, bpr](size_t x, size_t y) mutable {
                    uint8_t* ptr = &bitmap.data()[y * bpr + x / 8];
                    *ptr |= 1UL << (7 - x % 8);
                };

                int32_t breakline = 1;
                for (size_t y = 0; y < h; y++) {
                    for (size_t x = 0; x < w; x++) {
                        if (rgbaimage[y * w * 4 + x * 4 + 3] > 0) {
                            plot(x, y);
                        }
                    }
                }

                ss << std::format("        ");
                for (size_t c = 0; c < bitmap.size(); c++) {
                    ss << std::format("0x{:02x}{}", bitmap.data()[c], (c < bitmap.size() - 1) ? "," : "");
                    if ((breakline++ & 0xf) == 0) {
                        ss << std::format("\n        ");
                    }
                }

                ss << std::format("}}}};\n\n");

            } else {
                size_t bpr = ((w + 1) / 2);
                std::vector<uint8_t> bitmap;
                bitmap.assign(h * bpr, 0);

                ss << std::format("    static constexpr bool mono = false;\n");
                ss << std::format("    static constexpr size_t glyph_bitmap_width = {};\n", w);
                ss << std::format("    static constexpr size_t glyph_bitmap_height = {};\n", h);
                ss << std::format("    static constexpr size_t glyph_bitmap_stride = {};\n\n", bpr);

                ss << std::format("    static constexpr std::array<uint8_t, {}> glyph_bitmap{{{{\n", h * bpr);

                auto plot = [&bitmap, bpr](size_t x, size_t y, uint8_t a) mutable {
                    uint8_t* ptr = &bitmap.data()[y * bpr + x / 2];
                    *ptr |= (a >> 4) << ((1 - x % 2) * 4);
                };

                int32_t breakline = 1;
                for (size_t y = 0; y < h; y++) {
                    for (size_t x = 0; x < w; x++) {
                        plot(x, y, rgbaimage[y * w * 4 + x * 4 + 3]);
                    }
                }

                ss << std::format("        ");
                for (size_t c = 0; c < bitmap.size(); c++) {
                    ss << std::format("0x{:02x}{}", bitmap.data()[c], (c < bitmap.size() - 1) ? "," : "");
                    if ((breakline++ & 0xf) == 0) {
                        if (c == (bitmap.size() - 1)) {
                            ss << std::format("\n");
                        } else {
                            ss << std::format("\n        ");
                        }
                    }
                }

                ss << std::format("    }}}};\n\n");
            }

            ss << std::format("}};\n\n");
            ss << std::format("}}\n");

            // std::cout << ss.str();

            std::filesystem::path dir = config.out_dir;
            std::filesystem::path file = std::format("{}.h", name);

            std::ofstream out(dir / file);
            out << ss.str();

        } else {
            throw std::runtime_error("Can not open/decode font texture.");
        }
    } else {
        throw std::runtime_error("Font file does not contain texture file name.");
    }
    return 0;
}
