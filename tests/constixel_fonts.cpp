#include <filesystem>
#include <fstream>
#include <iostream>

#include "constixel.h"
#include "fontslist.h"

static constixel::image<constixel::format_8bit, 1024, 2048> fonts;

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
        fonts.draw_string_aa<T>(0, y_pos, text.c_str(), 1);
        y_pos += T::total_height;
        if (y_pos + 100 > fonts.height()) {
            fonts.sixel_to_cout();
            fonts.clear();
            y_pos = 0;
        }
    });

    fonts.sixel_to_cout();
    fonts.clear();
    y_pos = 0;

    for_each_type<MONO_HEADERS>([&]<std::size_t, typename T>() {
        std::string text;
        text.append(T::name);
        text.append(" ");
        text.append(T::style);
        text.append(" ");
        snprintf(size, 10, "%d", T::size);
        text.append(size);
        text.append(" ");
        text.append("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz01234567890!@#$%^&*()");
        fonts.draw_string_mono<T>(0, y_pos, text.c_str(), 1);
        y_pos += T::total_height;
        if (y_pos + 100 > fonts.height()) {
            fonts.sixel_to_cout();
            fonts.clear();
            y_pos = 0;
        }
    });

    fonts.sixel_to_cout();
}
