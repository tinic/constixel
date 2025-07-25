#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

#if __cplusplus >= 202302L
#include <format>
#include <print>
#endif  // #if __cplusplus >= 202302L

#define CONSTIXEL_ENABLE_COUT

#include "constixel.hpp"
#include "fontslist.hpp"

//#define PNG_OUTPUT

// -----------------------------------------------------------------------------
// for_each_type: invokes F::operator()<Idx, T>() for every Ts...
// -----------------------------------------------------------------------------
template <typename F, typename... Ts, std::size_t... Is>
constexpr void for_each_type_impl(F&& f, std::index_sequence<Is...>) {
    (f.template operator()<Is, Ts>(), ...);
}

template <typename... Ts, typename F>
constexpr void for_each_type(F&& f) {
    for_each_type_impl<F, Ts...>(std::forward<F>(f), std::index_sequence_for<Ts...>{});
}

// ------------------------------------------------------------
int main() {

#ifdef PNG_OUTPUT
    auto fonts = std::make_unique<constixel::image<constixel::format_8bit, 1024, 14000>> ();
#else // #ifdef PNG_OUTPUT
    auto fonts = std::make_unique<constixel::image<constixel::format_8bit, 1024, 2048>> ();
#endif // #ifdef PNG_OUTPUT

    using namespace constixel;
    // compile‑time index
    //    for_each_type<ALL_FONTS>(print_with_index{});

    char size[10] = {};
    // runtime counter you can update
    int32_t y_pos = 0;
    for_each_type<AA_HEADERS>([&]<std::size_t, typename T>() {
        std::string text;
        text.append(T::name);
        text.append(" ");
        text.append(T::style);
        text.append(" ");
        snprintf(size, 10, "%d", T::size);
        text.append(size);
        text.append(" ");
        text.append("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz01234567890!@#$%^&*()");
        fonts->draw_string_aa<T, true>(0, y_pos, text.c_str(), 1);
        y_pos += T::total_height;
#ifndef PNG_OUTPUT
        if (y_pos + 100 > fonts->height()) {
            fonts->sixel_to_cout();
            fonts->clear();
            y_pos = 0;
        }
#endif  // #ifndef PNG_OUTPUT
    });

    for_each_type<MONO_HEADERS>([&]<std::size_t, typename T>() {
        std::string text;
        text.append(T::name);
        text.append(" ");
        text.append(T::style);
        text.append(" ");
        snprintf(size, 10, "%d", T::size);
        text.append(size);
        text.append(" ");
        text.append("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz01234567890!@#$%^&*()░▒▓███▚▚╚═╝");
        fonts->draw_string_mono<T, true>(0, y_pos, text.c_str(), 1);
        y_pos += T::total_height;
#ifndef PNG_OUTPUT
        if (y_pos + 100 > fonts->height()) {
            fonts->sixel_to_cout();
            fonts->clear();
            y_pos = 0;
        }
#endif  // #ifndef PNG_OUTPUT
    });

#ifdef PNG_OUTPUT
    std::vector<char> out{};
    fonts->png([&out](char ch) mutable {
        out.push_back(ch);
    });

    std::ofstream file("constixel_fonts.png", std::ios::binary);
    file.write(reinterpret_cast<const char*>(out.data()), static_cast<std::streamsize>(out.size()));
#else  // #ifdef PNG_OUTPUT
    fonts->sixel_to_cout();
#endif  // #ifdef PNG_OUTPUT

#if __cplusplus >= 202302L
    for_each_type<AA_HEADERS>([&]<std::size_t, typename T>() {
        std::print("kerning_tree: {:8d} glyph_bitmap: {:8d} char_table: {:8d} glyph_tree: {:8d} total: {:8d} {} {} {} {} \n", sizeof(T::kerning_tree),
                   sizeof(T::glyph_bitmap), sizeof(T::char_table), sizeof(T::glyph_tree),
                   sizeof(T::kerning_tree) + sizeof(T::glyph_bitmap) + sizeof(T::char_table) + sizeof(T::glyph_tree), T::name, T::style, T::size,
                   T::mono ? "mono" : "aa");
    });

    for_each_type<MONO_HEADERS>([&]<std::size_t, typename T>() {
        std::print("kerning_tree: {:8d} glyph_bitmap: {:8d} char_table: {:8d} glyph_tree: {:8d} total: {:8d} {} {} {} {} \n", sizeof(T::kerning_tree),
                   sizeof(T::glyph_bitmap), sizeof(T::char_table), sizeof(T::glyph_tree),
                   sizeof(T::kerning_tree) + sizeof(T::glyph_bitmap) + sizeof(T::char_table) + sizeof(T::glyph_tree), T::name, T::style, T::size,
                   T::mono ? "mono" : "aa");
    });
#endif  // #if __cplusplus >= 202302L
}
