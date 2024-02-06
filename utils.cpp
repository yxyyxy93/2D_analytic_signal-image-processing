#include "utils.h"

#include <QVector>

// To overload the +, -, * operators with the input type QVector<QVector<double>>, point-wise operation
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] + v2[i][j];
        }
    }

    return result;
}
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const double& num){
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] + num;
        }
    }

    return result;
}

QVector<QVector<double>> operator-(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] - v2[i][j];
        }
    }

    return result;
}

QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] * v2[i][j];
        }
    }

    return result;
}
QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const double& num) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] * num;
        }
    }

    return result;
}
