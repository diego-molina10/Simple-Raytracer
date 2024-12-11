#include <vector>
#include <chrono>
#include "raytracing.h"

int main() {
    const int N = 1; 
    ofstream image("output.ppm"); 
    image << "P3\n" << Cw << " " << Ch << "\n255\n";

    auto start = chrono::high_resolution_clock::now();

    vector<vector<Vector>> color_grid((Ch / N) + 1, vector<Vector>((Cw / N) + 1));

    for (int y = -Ch / 2; y < Ch / 2; y += N) {
        for (int x = -Cw / 2; x < Cw / 2; x += N) {
            Vector D = CanvasToViewport(x, y);
            color_grid[(y + Ch / 2) / N][(x + Cw / 2) / N] = TraceRay(O, D, 1.0, numeric_limits<double>::infinity(), 3);
        }
    }

    for (int y = 0; y < Ch; ++y) {
        for (int x = 0; x < Cw; ++x) {
            int grid_x = x / N;
            int grid_y = y / N;
            Vector color = color_grid[grid_y][grid_x];

            image << (int)get<0>(color.v) << " "
                  << (int)get<1>(color.v) << " "
                  << (int)get<2>(color.v) << " ";
        }
        image << "\n";
    }

    image.close();

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Rendering time with subsampling (" << N <<") : " << elapsed.count() << " seconds\n";

    return 0;
}
