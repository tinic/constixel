#define CONSTIXEL_ENABLE_COUT

#include "constixel.hpp"
#include "fonts/ibmplexsans_bold_24_aa.hpp"
#include "fonts/ibmplexsans_semibold_18_aa.hpp"
#include "fonts/ibmplexsans_regular_12_aa.hpp"
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

using namespace constixel;
using UI = image<format_8bit, 1024, 768>;
using title_font = constixel::ibmplexsans_bold_24_aa;
using label_font = constixel::ibmplexsans_semibold_18_aa;
using small_font = constixel::ibmplexsans_regular_12_aa;

struct Point3D {
    float x, y, z;
    Point3D(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}
};

struct Point2D {
    int x, y;
    Point2D(int x = 0, int y = 0) : x(x), y(y) {}
};

class Matrix3D {
public:
    float m[4][4];
    
    Matrix3D() {
        identity();
    }
    
    void identity() {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                m[i][j] = (i == j) ? 1.0f : 0.0f;
            }
        }
    }
    
    static Matrix3D rotationY(float angle) {
        Matrix3D mat;
        float c = cos(angle);
        float s = sin(angle);
        mat.m[0][0] = c;  mat.m[0][2] = s;
        mat.m[2][0] = -s; mat.m[2][2] = c;
        return mat;
    }
    
    static Matrix3D rotationX(float angle) {
        Matrix3D mat;
        float c = cos(angle);
        float s = sin(angle);
        mat.m[1][1] = c;  mat.m[1][2] = -s;
        mat.m[2][1] = s;  mat.m[2][2] = c;
        return mat;
    }
    
    Point3D transform(const Point3D& p) const {
        return Point3D(
            m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z,
            m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z,
            m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z
        );
    }
    
    Matrix3D operator*(const Matrix3D& other) const {
        Matrix3D result;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                result.m[i][j] = 0;
                for (int k = 0; k < 4; k++) {
                    result.m[i][j] += m[i][k] * other.m[k][j];
                }
            }
        }
        return result;
    }
};

class Renderer3D {
private:
    UI& img;
    int centerX, centerY;
    float scale;
    float viewDistance;
    
public:
    Renderer3D(UI& image) : img(image) {
        centerX = img.width() / 2;
        centerY = img.height() / 2 - 50; // Offset up for labels
        scale = 180.0f;
        viewDistance = 8.0f;
    }
    
    Point2D project(const Point3D& p3d) const {
        // Simple perspective projection
        float perspective = viewDistance / (viewDistance + p3d.z);
        return Point2D(
            centerX + (int)(p3d.x * scale * perspective),
            centerY - (int)(p3d.y * scale * perspective)
        );
    }
    
    bool isVisible(const Point3D& p) const {
        Point2D p2d = project(p);
        return p2d.x >= 0 && p2d.x < img.width() && 
               p2d.y >= 0 && p2d.y < img.height();
    }
};

struct Colors {
    static uint8_t background;
    static uint8_t grid_light;
    static uint8_t grid_dark;
    static uint8_t surface_high;
    static uint8_t surface_mid;
    static uint8_t surface_low;
    static uint8_t text_primary;
    static uint8_t text_secondary;
    static uint8_t axes;
    
    static void init(const UI& img) {
        background = img.get_nearest_color(8, 12, 20);        // Dark blue-black
        grid_light = img.get_nearest_color(40, 50, 70);       // Light grid
        grid_dark = img.get_nearest_color(25, 30, 45);        // Dark grid
        surface_high = img.get_nearest_color(255, 100, 100);  // Red peaks
        surface_mid = img.get_nearest_color(100, 200, 255);   // Blue mid
        surface_low = img.get_nearest_color(50, 255, 150);    // Green valleys
        text_primary = img.get_nearest_color(255, 255, 255);  // White
        text_secondary = img.get_nearest_color(180, 190, 200); // Light gray
        axes = img.get_nearest_color(200, 200, 100);          // Yellow axes
    }
};

uint8_t Colors::background;
uint8_t Colors::grid_light;
uint8_t Colors::grid_dark;
uint8_t Colors::surface_high;
uint8_t Colors::surface_mid;
uint8_t Colors::surface_low;
uint8_t Colors::text_primary;
uint8_t Colors::text_secondary;
uint8_t Colors::axes;

float mathFunction(float x, float z) {
    // Impressive mathematical surface: combination of sine waves
    return 0.6f * sin(sqrt(x*x + z*z) * 2.0f) * exp(-0.3f * sqrt(x*x + z*z)) + 
           0.3f * sin(x * 3.0f) * cos(z * 2.0f);
}

uint8_t getHeightColor(float height, float minH, float maxH) {
    float normalized = (height - minH) / (maxH - minH);
    if (normalized < 0.33f) return Colors::surface_low;
    else if (normalized < 0.66f) return Colors::surface_mid;
    else return Colors::surface_high;
}

void drawSurface(UI& img, Renderer3D& renderer, const Matrix3D& transform) {
    const int resolution = 40;
    const float range = 3.0f;
    const float step = (2.0f * range) / resolution;
    
    float minHeight = 1000.0f, maxHeight = -1000.0f;
    
    // First pass: find height range
    for (int i = 0; i < resolution; i++) {
        for (int j = 0; j < resolution; j++) {
            float x = -range + i * step;
            float z = -range + j * step;
            float y = mathFunction(x, z);
            minHeight = std::min(minHeight, y);
            maxHeight = std::max(maxHeight, y);
        }
    }
    
    // Draw wireframe surface
    for (int i = 0; i < resolution - 1; i++) {
        for (int j = 0; j < resolution - 1; j++) {
            float x1 = -range + i * step;
            float z1 = -range + j * step;
            float x2 = -range + (i + 1) * step;
            float z2 = -range + (j + 1) * step;
            
            float y1 = mathFunction(x1, z1);
            float y2 = mathFunction(x2, z1);
            float y3 = mathFunction(x1, z2);
            
            Point3D p1 = transform.transform(Point3D(x1, y1, z1));
            Point3D p2 = transform.transform(Point3D(x2, y2, z1));
            Point3D p3 = transform.transform(Point3D(x1, y3, z2));
            
            if (renderer.isVisible(p1) || renderer.isVisible(p2)) {
                Point2D proj1 = renderer.project(p1);
                Point2D proj2 = renderer.project(p2);
                uint8_t color = getHeightColor((y1 + y2) / 2.0f, minHeight, maxHeight);
                img.draw_line_aa(proj1.x, proj1.y, proj2.x, proj2.y, color, 1);
            }
            
            if (renderer.isVisible(p1) || renderer.isVisible(p3)) {
                Point2D proj1 = renderer.project(p1);
                Point2D proj3 = renderer.project(p3);
                uint8_t color = getHeightColor((y1 + y3) / 2.0f, minHeight, maxHeight);
                img.draw_line_aa(proj1.x, proj1.y, proj3.x, proj3.y, color, 1);
            }
        }
    }
}

void drawAxes(UI& img, Renderer3D& renderer, const Matrix3D& transform) {
    // Draw 3D coordinate axes
    Point3D origin = transform.transform(Point3D(0, 0, 0));
    Point3D xAxis = transform.transform(Point3D(2.5f, 0, 0));
    Point3D yAxis = transform.transform(Point3D(0, 2.5f, 0));
    Point3D zAxis = transform.transform(Point3D(0, 0, 2.5f));
    
    Point2D oProj = renderer.project(origin);
    Point2D xProj = renderer.project(xAxis);
    Point2D yProj = renderer.project(yAxis);
    Point2D zProj = renderer.project(zAxis);
    
    // X axis (red)
    uint8_t xColor = img.get_nearest_color(255, 100, 100);
    img.draw_line_aa(oProj.x, oProj.y, xProj.x, xProj.y, xColor, 2);
    img.draw_string_aa<small_font>(xProj.x + 5, xProj.y - small_font::ascent, "X", Colors::text_secondary);
    
    // Y axis (green)
    uint8_t yColor = img.get_nearest_color(100, 255, 100);
    img.draw_line_aa(oProj.x, oProj.y, yProj.x, yProj.y, yColor, 2);
    img.draw_string_aa<small_font>(yProj.x + 5, yProj.y - small_font::ascent, "Y", Colors::text_secondary);
    
    // Z axis (blue)
    uint8_t zColor = img.get_nearest_color(100, 100, 255);
    img.draw_line_aa(oProj.x, oProj.y, zProj.x, zProj.y, zColor, 2);
    img.draw_string_aa<small_font>(zProj.x + 5, zProj.y - small_font::ascent, "Z", Colors::text_secondary);
}

void drawLabels(UI& img) {
    // Title
    const char* title = "3D Mathematical Surface Visualization";
    int titleWidth = img.string_width<title_font>(title);
    img.draw_string_aa<title_font>((img.width() - titleWidth) / 2, 40 - title_font::ascent, title, Colors::text_primary);
    
    // Function equation
    const char* equation = "f(x,z) = 0.6×sin(√(x²+z²)×2)×e^(-0.3√(x²+z²)) + 0.3×sin(3x)×cos(2z)";
    int eqWidth = img.string_width<small_font>(equation);
    img.draw_string_aa<small_font>((img.width() - eqWidth) / 2, 70 - small_font::ascent, equation, Colors::text_secondary);
    
    // Legend
    int legendY = img.height() - 120;
    img.draw_string_aa<label_font>(50, legendY - label_font::ascent, "Height Legend:", Colors::text_primary);
    
    // Color legend boxes
    img.fill_rect(50, legendY + 10, 20, 15, Colors::surface_high);
    img.draw_string_aa<small_font>(80, legendY + 15 - small_font::ascent / 2, "High", Colors::text_secondary);
    
    img.fill_rect(50, legendY + 35, 20, 15, Colors::surface_mid);
    img.draw_string_aa<small_font>(80, legendY + 40 - small_font::ascent / 2, "Medium", Colors::text_secondary);
    
    img.fill_rect(50, legendY + 60, 20, 15, Colors::surface_low);
    img.draw_string_aa<small_font>(80, legendY + 65 - small_font::ascent / 2, "Low", Colors::text_secondary);
    
    // Rotation info
    const char* rotInfo = "Rotation: Y-axis animation with X-tilt";
    img.draw_string_aa<small_font>(img.width() - img.string_width<small_font>(rotInfo) - 50, 
                                   img.height() - 30 - small_font::ascent, rotInfo, Colors::text_secondary);
}

int main() {
    static UI img;
    img.clear();
    
    Colors::init(img);
    
    // Clear background
    img.fill_rect(0, 0, img.width(), img.height(), Colors::background);
    
    // Create animated rotation (simulating frame from animation)
    float rotationY = 0.8f; // Approximately 45 degrees
    float rotationX = -0.3f; // Slight downward tilt
    
    Matrix3D transform = Matrix3D::rotationY(rotationY) * Matrix3D::rotationX(rotationX);
    
    Renderer3D renderer(img);
    
    // Draw the surface
    drawSurface(img, renderer, transform);
    
    // Draw coordinate axes
    drawAxes(img, renderer, transform);
    
    // Add labels and title
    drawLabels(img);
    
    // Output
    img.sixel_to_cout<1>();
    
    return 0;
}