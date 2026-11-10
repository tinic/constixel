#pragma once

#include <fstream>
#include <iostream>

namespace example {
// Pass a filename to render a reproducible PNG without terminal escape codes.
// With no argument, preserve the example's terminal or default-file output.
template <typename Image>
int output(const Image& img, int argc, char** argv, const char* default_path = nullptr) {
    if (argc > 2) {
        std::cerr << "Usage: " << argv[0] << " [output.png]\n";
        return 1;
    }
    const char* path = argc == 2 ? argv[1] : default_path;
    if (path == nullptr) {
        img.template sixel_to_cout<1>();
        return 0;
    }
    std::ofstream file(path, std::ios::binary);
    if (file) {
        img.png([&](char ch) {
            file.put(ch);
        });
        file.close();
    }
    if (!file) {
        std::cerr << "Could not write PNG: " << path << '\n';
        return 1;
    }
    return 0;
}
}  // namespace example
