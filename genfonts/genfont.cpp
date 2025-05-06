
#include "../constixel.h"


#include <string>
#include <string_view>
#include <vector>
#include <fstream>
#include <sstream>
#include <charconv>
#include <unordered_map>
#include <filesystem>
#include <cstdint>
#include <format>

#include "fontbm/src/external/cxxopts.hpp"
#include "fontbm/src/external/lodepng/lodepng.h"

struct Char {
    uint32_t id;
    uint16_t x, y, width, height;
    int16_t xoffset, yoffset, xadvance;
    uint16_t page;
    uint8_t chnl;
};

struct Font {
    std::string face;
    int size;
    int lineHeight;
    int base;
    int scaleW, scaleH;
    int pages;
    std::string page_file;
    std::vector<Char> chars;
    int smooth;
};

static std::string_view trim(std::string_view sv) {
    const auto begin = sv.find_first_not_of(' ');
    const auto   end = sv.find_last_not_of(' ');
    return begin == sv.npos ? "" : sv.substr(begin, end - begin + 1);
}

static std::unordered_map<std::string_view, std::string_view>
split_kv(std::string_view line) {
    std::unordered_map<std::string_view, std::string_view> kv;
    for (std::size_t i = 0; i < line.size();) {
        while (i < line.size() && line[i] == ' ') ++i;               // skip blanks

        const std::size_t key_beg = i;
        while (i < line.size() && line[i] != '=') ++i;
        if (i == line.size()) break;
        const std::size_t key_len = i++ - key_beg;                   // past '='

        std::string_view val;
        if (i < line.size() && line[i] == '"') {                     // quoted value
            const std::size_t val_beg = ++i;
            while (i < line.size() && line[i] != '"') ++i;
            val = line.substr(val_beg, i - val_beg);
            if (i < line.size()) ++i;                                // skip closing quote
        } else {                                                     // unquoted value
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
    std::string   line;
    while (std::getline(in, line)) {
        const auto sv   = trim(line);
        const auto sp   = sv.find(' ');
        const auto cmd  = sv.substr(0, sp);
        const auto args = sp == sv.npos ? "" : sv.substr(sp + 1);
        const auto kv   = split_kv(args);

        if (cmd == "info") {
            font.face = std::string(kv.at("face"));
            font.size = to_int(kv.at("size"));
            font.smooth = to_int(kv.at("smooth"));
        } else if (cmd == "common") {
            font.lineHeight = to_int(kv.at("lineHeight"));
            font.base       = to_int(kv.at("base"));
            font.scaleW     = to_int(kv.at("scaleW"));
            font.scaleH     = to_int(kv.at("scaleH"));
            font.pages      = to_int(kv.at("pages"));
        } else if (cmd == "page") {
            font.page_file = std::string(kv.at("file"));
        } else if (cmd == "char") {
            Char c{};
            c.id       = to_int(kv.at("id"));
            c.x        = to_int(kv.at("x"));
            c.y        = to_int(kv.at("y"));
            c.width    = to_int(kv.at("width"));
            c.height   = to_int(kv.at("height"));
            c.xoffset  = to_int(kv.at("xoffset"));
            c.yoffset  = to_int(kv.at("yoffset"));
            c.xadvance = to_int(kv.at("xadvance"));
            c.page     = to_int(kv.at("page"));
            c.chnl     = to_int(kv.at("chnl"));
            font.chars.push_back(c);
        }
    }
    return font;
}

struct Config {
    std::string font_file;
    std::string out_dir;
} config;

static void parseCommandLine(int argc, char* argv[])
{
    try
    {
        cxxopts::Options options("genfont", "Command line font generator for constixel, compatible with bmfont");
        options.add_options()
            ("font-file", "path to ttf/otf file, required", cxxopts::value<std::string>(config.font_file))
            ("out-dir", "path to output path, required", cxxopts::value<std::string>(config.out_dir));
        auto result = options.parse(argc, argv);

        if (!result.count("font-file"))
        {
            std::cout << options.help() << std::endl;
            throw std::runtime_error("--font-file required");
        }

        if (!result.count("out-dir"))
        {
            std::cout << options.help() << std::endl;
            throw std::runtime_error("--out-dir required");
        }
    }
    catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        throw std::exception();
    }
}

std::string slugify(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return c == '-' ? '_' : (c == ' ' ? '_' : std::tolower(c)); });
    return s;
}

int main(int argc, char* argv[]) {
    parseCommandLine(argc, argv);
    Font font = parse_fnt(config.font_file);
    std::print("Converting font: {}\n",config.font_file);
    if (font.page_file.size() > 0) {
        if (font.face.size() == 0) {
            throw std::runtime_error("Font has no name.");
        }
        std::vector<uint8_t> rgbaimage;
        uint32_t w = 0;
        uint32_t h = 0;
        auto path = std::filesystem::path{config.font_file}.parent_path() / font.page_file;
        if (lodepng::decode(rgbaimage, w, h, path) == 0) {
            bool is_mono = font.smooth == 1;

            std::stringstream ss;

            std::string name = std::format("{}_{}", slugify(std::filesystem::path(config.font_file).stem().string()), is_mono ? "mono" : "grey");

            ss << std::format("namespace constixel {{\n\n", slugify(font.face), abs(font.size));
            ss << std::format("struct {} {{\n\n", name);

            uint32_t max_id = 0;
            for (size_t c = 0; c <font.chars.size(); c++) {
                max_id = std::max(font.chars[c].id, max_id);
            }
            if (font.chars.size() < 65536 && max_id < 65536) {
                ss << std::format("    static constexpr std::array<std::pair<uint16_t, uint16_t>> glyph_table{{{{\n");
            } else {
                ss << std::format("    static constexpr std::array<std::pair<uint32_t, uint32_t>> glyph_table{{{{\n");
            }
            for (size_t c = 0; c <font.chars.size(); c++) {
                ss << std::format("        {{ 0x{:06x}, 0x{:06} }},\n", font.chars[c].id, c);
            }
            ss << std::format("    }}}};\n\n");
            if (font.chars.size() < 65536 && max_id < 65536) {
                ss << std::format("    static constexpr hextree<hextree<0, uint16_t>, uint16_t> glyph_tree(glyph_table);\n\n");
            } else {
                ss << std::format("    static constexpr hextree<hextree<0, uint32_t>, uint32_t> glyph_tree(glyph_table);\n\n");
            }

            ss << std::format("    static constexpr std::array<char_info, {}> char_table({{{{\n", font.chars.size());
                 //                  { 0x000000,    0,    0,    0,    0,    6,    0,    0 },
            ss << std::format("      // | unicode|    x|    y|    w|    h| xadv| xoff| yoff|\n");
            for (size_t c = 0; c <font.chars.size(); c++) {
                ss << std::format("        {{ 0x{:06x}, {:4}, {:4}, {:4}, {:4}, {:4}, {:4}, {:4} }},\n", 
                    font.chars[c].id, 
                    font.chars[c].x, 
                    font.chars[c].y, 
                    font.chars[c].width, 
                    font.chars[c].height,
                    font.chars[c].xadvance,
                    font.chars[c].xoffset,
                    font.chars[c].yoffset
                 );
            }

            ss << std::format("    }}}};\n\n");

            if (font.smooth == 1) {
                // mono
                
                size_t bpr = ((w + 7) / 8);
                std::vector<uint8_t> bitmap;
                bitmap.assign(h * bpr, 0);

                ss << std::format("    static constexpr mono = true;\n\n");

                ss << std::format("    static constexpr std::array<uint8_t, {}> glyph_bitmap{{{{\n", h * bpr);

                auto plot = [&bitmap, bpr] (size_t x, size_t y) mutable {
                    uint8_t *ptr = &bitmap.data()[y * bpr + x / 8];
                    *ptr |= 1UL << (7 - x%8);
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
                for (uint8_t p : bitmap) {
                    ss << std::format("0x{:02x},", p);
                    if ((breakline++&0xf)==0) {
                        ss << std::format("\n        ");
                    }
                }

                ss << std::format("}}}};\n\n");

            } else {

                size_t bpr = ((w + 1) / 2);
                std::vector<uint8_t> bitmap;
                bitmap.assign(h * bpr, 0);

                ss << std::format("    static constexpr mono = false;\n\n");

                ss << std::format("    static constexpr std::array<uint8_t, {}> glyph_bitmap{{{{\n", h * bpr);

                auto plot = [&bitmap, bpr] (size_t x, size_t y, uint8_t a) mutable {
                    uint8_t *ptr = &bitmap.data()[y * bpr + x / 2];
                    *ptr |= (a >> 4) << ((1 - x%2) * 4);
                };

                int32_t breakline = 1;
                for (size_t y = 0; y < h; y++) {
                    for (size_t x = 0; x < w; x++) {
                        plot(x, y, rgbaimage[y * w * 4 + x * 4 + 0]);
                    }
                }

                ss << std::format("        ");
                for (uint8_t p : bitmap) {
                    ss << std::format("0x{:02x},", p);
                    if ((breakline++&0xf)==0) {
                        ss << std::format("\n        ");
                    }
                }

                ss << std::format("    }}}};\n\n");
            }

            ss << std::format("}};\n\n");
            ss << std::format("}};\n");

            //std::cout << ss.str();

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

/*
cd ../otf ; ../build/genfont --font-file SF-Compact-Display-Regular.otf.fnt --out-dir ../../fonts ; cd ../build
*/
