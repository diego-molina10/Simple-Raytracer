#include <vector>
#include <chrono>
#include "raytracing.h"

int main() {
    const int N = 1; // Subsampling factor: N determines the size of the grid cells (N=1 means no subsampling)
    ofstream image("output.ppm"); 
    image << "P3\n" << Cw << " " << Ch << "\n255\n";

    auto start = chrono::high_resolution_clock::now();

    // Creates a grid to store colors for subsampled regions.
    // Each cell in the grid represents the color of a block of size N x N pixels in the final image.
    vector<vector<Vector>> color_grid((Ch / N) + 1, vector<Vector>((Cw / N) + 1));

    for (int y = -Ch / 2; y < Ch / 2; y += N) {
        for (int x = -Cw / 2; x < Cw / 2; x += N) {
            Vector D = CanvasToViewport(x, y);

            // Trace the ray and store the computed color in the corresponding grid cell.
            color_grid[(y + Ch / 2) / N][(x + Cw / 2) / N] = TraceRay(O, D, 1.0, numeric_limits<double>::infinity(), 3);
        }
    }

    // Generates the final image by mapping every pixel to the color of its corresponding grid cell
    for (int y = 0; y < Ch; ++y) {
        for (int x = 0; x < Cw; ++x) {

            // Identifies the grid cell corresponding to the current pixel
            int grid_x = x / N;
            int grid_y = y / N;

            // Uses the precomputed color from the grid
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
