#include "as_2d_process_class.h"
#include <QVector>

int main() {
    // Step 1: Generate a random 2D image
    int rows(256);
    int cols(256);
    QVector<QVector<double>> randomImage(rows, QVector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            randomImage[i][j] = static_cast<double>(std::rand()) / RAND_MAX; // Generate random double between 0 and 1
        }
    }

    // Step 2: Instantiate AS_2D_process_class
    int masksize = 10;
    int sf = 1.0;
    int sc = 2.0;
    AS_2D_process_class processor(masksize=10, sf=sf, sc=sc); // Example parameters

    // Step 3: Call processData
    processor.processData(randomImage);

    return 0;
}

QVector<QVector<double>> generateRandomImage(int rows, int cols) {
    QVector<QVector<double>> image(rows, QVector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            image[i][j] = static_cast<double>(std::rand()) / RAND_MAX; // Generate random double between 0 and 1
        }
    }
    return image;
}
