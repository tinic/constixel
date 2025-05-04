# constixel

constixel is a constexpr minimal graphics library with the ability to output to a terminal using the sixel protocol. 

Primary features and goals:

- Completely constexpr. All operations, including the sixel output stream can be generated during compilation.
- No dynamic allocations. The backbuffer and interal data structures can live as global static variables.
- Minimalistic interface and single header implementation so it can be used without fuzz in any modern C++ project.
- Ability to render text from fontbm generated font files. All constexpr safe.
- 1, 2, 4 and 8bit palette based back buffers for minimal memory usage.

Applications:

- Interface rendering on embedded devices.
- Send graphical data to your terminal through a serial connection from any embedded device.
- Help debug dynamic memory handling issues in any complex C++ projects.
- Add graphical output to unit tests.
- Programmatically render static assets.
