#ifndef UTILS_H
#define UTILS_H

#endif // UTILS_H

#include <QVector>


// To overload the +, -, * operators with the input type QVector<QVector<double>>, point-wise operation
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2);
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const double& num);

QVector<QVector<double>> operator-(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2);

QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2);
QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const double& num);
