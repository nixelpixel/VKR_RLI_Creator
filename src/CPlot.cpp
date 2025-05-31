#include "CPlot.hpp"
#include "helpers.hpp"
#include <assert.h>
#include <stdlib.h>

CPlot::CPlot(void) {
    Width = 0;
    Height = 0;
    pFrameBuffer = NULL;
}

CPlot::~CPlot(void) {
    if (pFrameBuffer) {
        free(pFrameBuffer);
        pFrameBuffer = NULL;
    }
}

//==================================
// plot2d()
//==================================
int CPlot::plot2d(float *data, float ymin, float ymax, int len) {
    // Clear(240);

    float k = 1.0f - frameWidth;

    ymax -= ymin;

    for (int i = 0; i < len - 1; i++) {
        float y1 = 2 * (data[i] - ymin) / ymax - 1.0f;
        float y2 = 2 * (data[i + 1] - ymin) / ymax - 1.0f;
        float x1 = 2 * (float)i / len - 1.0f;
        float x2 = 2 * (float)(i + 1) / len - 1.0f;

        if (y1 > 1.0f)
            y1 = 1.0f;
        if (y2 > 1.0f)
            y2 = 1.0f;

        x1 *= k;
        x2 *= k;
        y1 *= k;
        y2 *= k;

        DrawLine(x1, y1, x2, y2, 100);
    }

    // DrawLine( -1.0f, 0.0f, 1.0f, 0.0f, 150);

    return true;
}

//==================================
// Clear()
//==================================
void CPlot::Clear(uint8_t color) {
    for (int y = 0; y < Height; y++) {
        for (int x = 0; x < Width; x++) {
            pFrameBuffer[x + Width * y] = color;
        }
    }
}

//==================================
// SetPixel()
//==================================
void CPlot::SetPixelf(float fx, float fy, uint8_t color) {
    int x = Width / 2.0f + fx * Width / 2.0f;
    int y = Height / 2.0f - fy * Height / 2.0f;

    pFrameBuffer[x + Width * y] = color;
}

void CPlot::SetPixel(int x, int y, uint8_t color) {
    assert(!(x > Width));
    assert(!(y > Height));

    pFrameBuffer[x + Width * y] = color;
}

//==================================
// Initialize()
//==================================
void CPlot::DrawLine(float fx1, float fy1, float fx2, float fy2,
                     uint8_t color) {
    int x1 = Width / 2.0f + fx1 * Width / 2.0f;
    int y1 = Height / 2.0f - fy1 * Height / 2.0f;
    int x2 = Width / 2.0f + fx2 * Width / 2.0f;
    int y2 = Height / 2.0f - fy2 * Height / 2.0f;

    int i;
    int x;
    int y;
    int xinc;
    int yinc;
    int dx;
    int dy;
    int e;

    /* swap x1,y1  with x2,y2 */
    if (x1 > x2) {
        dx = x1;
        x1 = x2;
        x2 = dx;
        dy = y1;
        y1 = y2;
        y2 = dy;
    }

    dx = x2 - x1;
    dy = y2 - y1;

    x = x1;
    y = y1;

    if (dx < 0) {
        xinc = -1;
        dx = -dx;
    } else {
        xinc = 1;
    }

    if (dy < 0) {
        yinc = -1;
        dy = -dy;
    } else {
        yinc = 1;
    }

    if (dx > dy) {
        e = dy - dx;
        for (i = 0; i <= dx; i++) {
            SetPixel(x, y, color);
            if (e >= 0) {
                e -= dx;
                y += yinc;
            }

            e += dy;
            x += xinc;
        }
    } else {
        e = dx - dy;
        for (i = 0; i <= dy; i++) {
            SetPixel(x, y, color);
            if (e >= 0) {
                e -= dy;
                x += xinc;
            }

            e += dx;
            y += yinc;
        }
    }
}

//==================================
// Initialize()
//==================================
int CPlot::Initialize(int width, int height) {

    Width = width;
    Height = height;
    frameWidth = 0.15f;

    pFrameBuffer =
        (uint8_t *)malloc(Width * Height * sizeof(uint8_t) * BYTES_PER_PIXEL);

    Clear(255);

    return true;
}

void CPlot::DrawGrid(int nx, int ny, float xmin, float xmax, float ymin,
                     float ymax) {

    float k = 1.0f - frameWidth;
    float step = k * 2.0 / nx;

    float dy = (ymax - ymin) / ny;
    float dx = (xmax - xmin) / nx;
    char buf[16];
    char xformat[16];
    char yformat[16];

    if (xmax > 10000) {
        strcpy(xformat, "%.2e");
    } else {
        strcpy(xformat, "%.2f");
    }

    if (ymax > 10000) {
        strcpy(yformat, "%.2e");
    } else {
        strcpy(yformat, "%.2f");
    }

    for (float x = -k; x <= k + 0.01; x += step) {

        DrawLine(x, -k, x, k, 200);

        sprintf(buf, xformat, xmin);
        AddLabelf(buf, x, -1.0f + (k / ny));
        xmin += dx;
    }

    step = k * 2.0 / ny;
    for (float y = -k; y <= k + 0.01; y += step) {

        DrawLine(-k, y, k, y, 200);
        sprintf(buf, yformat, ymin);
        AddLabelf(buf, -1, y);
        ymin += dy;
    }
}

void CPlot::RenderChar(uint8_t *c, int index) {
    int x, y;
    uint8_t set;

    uint8_t *bitmap = font8x8_basic[index];

    for (x = 0; x < 8; x++) {
        for (y = 0; y < 8; y++) {
            set = bitmap[x] & 1 << y;

            if (set) {
                set = 255;
            }
            c[(x * 8) + y] = ~set;
            // printf("%c", set ? 'X' : ' ');
        }
        // printf("\n");
    }
}

void CPlot::AddLabelf(char *label, float fx, float fy) {
    int x = Width / 2.0f + fx * Width / 2.0f;
    int y = Height / 2.0f - fy * Height / 2.0f;
    AddLabel(label, x, y);
}
void CPlot::AddLabel(char *label, int x_, int y_) {

    uint8_t *data = pFrameBuffer;
    uint8_t *c;
    int len = strlen(label);
    int x = x_;
    c = (uint8_t *)malloc(8 * 8);

    y_ -= 4; // By char center

    for (int i = 0; i < len; i++) {

        if (label[i] == '\n') {
            data += Width * 8;
            x = x_;
            continue;
        }

        RenderChar(c, (int)label[i]);
        for (int y = 0; y < 8; y++) {
            memcpy(data + (Width * (y + y_)) + x, c + (8 * y), 8);
        }
        x += 8;
    }

    free(c);
}

void CPlot::PlotChart(float *ydata, int len) {
    PlotChart((float)0, (float)len, ydata, len);
}

void CPlot::PlotChart(float xmin, float xmax, float *ydata, int len) {
    Clear(240);

    float ymax = std::numeric_limits<float>::min();
    float ymin = std::numeric_limits<float>::max();

    for (int i = 0; i < len; i++) {
        if (ymax < ydata[i]) {
            ymax = ydata[i];
        }

        if (ymin > ydata[i]) {
            ymin = ydata[i];
        }
    }

    DrawGrid(10, 10, xmin, xmax, ymin, ymax);

    plot2d(ydata, ymin, ymax, len);
}

void CPlot::PlotChart(float *xdata, float *ydata, int len) {

    Clear(240);

    float ymax = std::numeric_limits<float>::min();
    float ymin = std::numeric_limits<float>::max();

    float xmax = std::numeric_limits<float>::min();
    float xmin = std::numeric_limits<float>::max();

    for (int i = 0; i < len; i++) {
        if (ymax < ydata[i]) {
            ymax = ydata[i];
        }

        if (ymin > ydata[i]) {
            ymin = ydata[i];
        }

        if (xmax < xdata[i]) {
            xmax = xdata[i];
        }

        if (xmin > xdata[i]) {
            xmin = xdata[i];
        }
    }

    DrawGrid(10, 10, xmin, xmax, ymin, ymax);

    plot2d(ydata, ymin, ymax, len);
}
