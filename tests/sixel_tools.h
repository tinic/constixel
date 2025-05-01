#ifndef _SIXEL_TOOLS_H_
#define _SIXEL_TOOLS_H_

#include <sys/ioctl.h>
#include <termios.h>
#include <unistd.h>

#include "sixel.h"

namespace sixel {

template <template <size_t, size_t> class T, size_t W, size_t H>
class progressbar {
   public:
    progressbar(auto &_image) : image(_image) {
    }

    void start(int32_t _w, int32_t _h) {
        charw = _w;
        charh = _h;
        cellw = 0;
        cellh = 0;
        scroll_fix(32, cellw, cellh);
        printf("\0337");
    }

    void update(float v) {
        printf("\0338");

        int32_t x = 0;
        int32_t y = 0;
        int32_t w = cellw * charw;
        int32_t h = cellh * charh;
        image.fillrect(x, y, w, h, 1);
        image.fillrect(x + 2, y + 2, w - 4, h - 4, 0);
        image.fillrect(x + 4, y + 4, static_cast<int32_t>(ceilf(v * static_cast<float>(w - 8))), h - 8, 3);

        image.sixel(
            [](uint8_t ch) {
                putc(ch, stdout);
            },
            sixel::rect<int32_t>{0, 0, w, h});

        printf("\0338");
        printf("\033[%dC", (w / cellw) + 1);
        printf("%4.1f%%", (v * 100.0f));
    }

    void end() {
        for (size_t c = 0; c < charh; c++) {
            printf("\n");
        }
    }

   private:
    int32_t charw = 0;
    int32_t charh = 0;
    int32_t cellw = 0;
    int32_t cellh = 0;

    sixel::image<T, W, H> image;

    void scroll_fix(int32_t pixelheight, int32_t &cel_w, int32_t &cel_h) {
        struct winsize size = {0, 0, 0, 0};
        struct termios old_termios;
        struct termios new_termios;
        int row = 0;
        int col = 0;
        int32_t cellheight;
        int32_t scroll;

        ioctl(STDOUT_FILENO, TIOCGWINSZ, &size);
        if (size.ws_ypixel <= 0) {
            printf("\033[H\0337");
            return;
        }

        cel_w = size.ws_xpixel / size.ws_col;
        cel_h = size.ws_ypixel / size.ws_row;

        /* set the terminal to cbreak mode */
        tcgetattr(STDIN_FILENO, &old_termios);
        memcpy(&new_termios, &old_termios, sizeof(old_termios));
        new_termios.c_lflag &= ~(ECHO | ICANON);
        new_termios.c_cc[VMIN] = 1;
        new_termios.c_cc[VTIME] = 0;
        tcsetattr(STDIN_FILENO, TCSAFLUSH, &new_termios);

        /* request cursor position report */
        printf("\033[6n");
        if (usleep(1000) != (-1)) { /* wait 1 sec */
            if (scanf("\033[%d;%dR", &row, &col) == 2) {
                cellheight = pixelheight * size.ws_row / size.ws_ypixel + 1;
                scroll = cellheight + row - size.ws_row + 1;
                printf("\033[%dS\033[%dA", scroll, scroll);
                printf("\0337");
            } else {
                printf("\033[H\0337");
            }
        }

        tcsetattr(STDIN_FILENO, TCSAFLUSH, &old_termios);
    }
};
};  // namespace sixel

#endif  // #ifndef _SIXEL_TOOLS_H_
