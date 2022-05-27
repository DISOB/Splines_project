#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cfloat>
#include <ctime>


using namespace std;

double lower_b;
double upper_b;

// Все переменные для первого сплайна
int N1; // Количество точек первого сплайна
double *xMas1, *yMas1; // Векторы с координатами
double *h1, *l1, *delta1, *lambda1, *c1, *d1, *b1; // Векторы с коэффициентами


// Все переменные для второго сплайна
int N2; // Количество точек первого сплайна
double *xMas2, *yMas2; // Векторы с координатами
double *h2, *l2, *delta2, *lambda2, *c2, *d2, *b2; // Векторы с коэффициентами


void create_vectors_1() {
    xMas1 = new double[N1];
    yMas1 = new double[N1];
    h1 = new double[N1];
    l1 = new double[N1];
    delta1 = new double[N1];
    lambda1 = new double[N1];
    c1 = new double[N1];
    d1 = new double[N1];
    b1 = new double[N1];
}

void create_vectors_2() {
    xMas2 = new double[N2];
    yMas2 = new double[N2];
    h2 = new double[N2];
    l2 = new double[N2];
    delta2 = new double[N2];
    lambda2 = new double[N2];
    c2 = new double[N2];
    d2 = new double[N2];
    b2 = new double[N2];
}

// Функция печати коэффициентов
void printresult(double *yMas, double *b, double *c, double *d, int N) {
    int k = 0;
    printf("\nA[k]\tB[k]\tC[k]\tD[k]\n");
    for (k = 1; k < N; k++) {
        printf("%f\t%f\t%f\t%f\n", yMas[k], b[k], c[k], d[k]);
    }
}

double get_value(double *xMas, double *yMas, double *b, double *c, double *d, int N, double x) {
    int k;
    //find k, where s in [x_k-1; x_k]
    for (k = 1; k < N; k++) {
        if (x >= xMas[k - 1] && x <= xMas[k]) {
            break;
        }
    }
    double F = yMas[k] + b[k] * (x - xMas[k]) + c[k] * pow(x - xMas[k], 2) +
               d[k] * pow(x - xMas[k], 3); // Подсчёт значения сплайна в точке
    return F;
}

void readPoints(double *&xMas, double *&yMas, int N) {
    double x, y;
    for (int i = 0; i < N; i++) {
        cin >> x >> y;
        xMas[i] = x;
        yMas[i] = y;
    }
}

int get_spline_coefficients(double *&xMas, double *&yMas, double *&h, double *&l, double *&lambda, double *&delta,
                            double *&b,
                            double *&c, double *&d, int N) {
    int k = 0;
    // Следующие преобразования деляются для удобства
    for (k = 1; k < N; k++) {
        h[k] = xMas[k] - xMas[k - 1];
        if (h[k] == 0) {
            printf("\nError, x[%d]=x[%d]\n", k, k - 1);
            return 0;
        }
        l[k] = (yMas[k] - yMas[k - 1]) / h[k];
    }

    // Вводим прогоночные коэффициенты
    delta[1] = (-h[2]) / (2 * (h[1] + h[2]));
    lambda[1] = 1.5 * (l[2] - l[1]) / (h[1] + h[2]);
    for (k = 3; k < N; k++) {
        delta[k - 1] = (-h[k]) / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]);
        lambda[k - 1] = (3 * l[k] - 3 * l[k - 1] - h[k - 1] * lambda[k - 2]) /
                        (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]);
    }

    c[0] = 0;
    c[N - 1] = 0;
    // По фурмуле обратной прогонки находим все коэф. С[k]
    for (k = N - 1; k >= 2; k--) {
        c[k - 1] = delta[k - 1] * c[k] + lambda[k - 1];
    }

    // Зная C[k], находим b[k] и d[k]
    for (k = 1; k < N; k++) {
        d[k] = (c[k] - c[k - 1]) / (3 * h[k]);
        b[k] = l[k] + (2 * c[k] * h[k] + h[k] * c[k - 1]) / 3;
    }
    return 1;
}

// Функции для нахождения значений 1 и 2 сплайнов в точках
double F1(double x) {
    return get_value(xMas1, yMas1, b1, c1, d1, N1, x);
}

double F2(double x) {
    return get_value(xMas2, yMas2, b2, c2, d2, N2, x);
}

double dist(double x, double y) {
    if (x < (lower_b)) {
        return -1e11 * (x + y) + (dist(max(x, lower_b), max(y, lower_b)) + 1e11 * max(x, lower_b) * max(y, lower_b));
    }
    if (x > (upper_b)) {
        return 1e11 * (x + y) + (dist(min(x, upper_b), min(y, upper_b)) - 1e11 * min(x, upper_b) * min(y, upper_b));
    }
    return sqrt((F1(x) - F2(y)) * (F1(x) - F2(y)) + (x - y) * (x - y));
}

pair<double, double> derivative(double x, double y, double eps = 1e-10) {
    double f_x = (dist(x + eps, y) - dist(x, y)) / eps;
    double f_y = (dist(x, y + eps) - dist(x, y)) / eps;
    return make_pair(f_x, f_y);
}

pair<pair<double, double>, double>
adam(double init_value_x, double init_value_y, int n_iter, double alpha = 0.0002, double beta_1 = 0.9,
     double beta_2 = 0.999,
     double eps = 1e-8) {
    auto x = make_pair(init_value_x, init_value_y);
    double distance = dist(x.first, x.second);
    auto m = make_pair((double) 0, (double) 0);
    auto v = make_pair((double) 0, (double) 0);
    double beta_1_pow = beta_1;
    double beta_2_pow = beta_2;
    for (int i = 0; i < n_iter; ++i) {
        auto g = derivative(x.first, x.second);
        m.first = beta_1 * m.first + (1 - beta_1) * g.first;
        v.first = beta_2 * v.first + (1 - beta_2) * g.first * g.first;
        double m_hat_x = m.first / (1 - beta_1_pow);
        double v_hat_x = v.first / (1 - beta_2_pow);
        x.first = x.first - alpha * m_hat_x / (sqrt(v_hat_x) + eps);


        m.second = beta_1 * m.second + (1 - beta_1) * g.second;
        v.second = beta_2 * v.second + (1 - beta_2) * g.second * g.second;
        double m_hat_y = m.second / (1 - beta_1_pow);
        double v_hat_y = v.second / (1 - beta_2_pow);
        x.second = x.second - alpha * m_hat_y / (sqrt(v_hat_y) + eps);

        beta_1_pow *= beta_1;
        beta_2_pow *= beta_2;

        distance = dist(x.first, x.second);
//        cout << distance << "\n";
    }
//    cout << x.first << " " << x.second << " " << distance << endl;
    return make_pair(x, distance);
}

double fRand(double fMin, double fMax) {
    double f = (double) rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void find_min(int iter, double eps = 1e-5) {
    pair<double, double> x;
    double best_dist = DBL_MAX;
    for (int i = 0; i < iter; ++i) {
        double init_value_x = fRand(lower_b, upper_b);
        double init_value_y = fRand(lower_b, upper_b);
        pair<pair<double, double>, double> temp = adam(init_value_x, init_value_y, 200000);
        auto x_temp = temp.first;
        double dist_temp = temp.second;
        if (dist_temp < best_dist) {
            x = x_temp;
            best_dist = dist_temp;
        }
        if (i == iter / 2) {
            cout << "Just a little bit left.\n";
        }
    }
    if (best_dist - (int) best_dist < eps) {
        cout << "The splines intersect!" << endl;
        cout << "Point, where splines intersect: " << (x.first + x.second) / 2 << endl;
    } else {
        cout << "Minimum distance between splines: " << best_dist << endl;
        cout << "Point where the distance is minimal: " << endl << x.first << " - point on 1st spline" << endl;
        cout << x.second << " - point on 2nd spline";
    }
}

int main() {
    setlocale(LC_ALL, "Rus");
    srand(time(nullptr));

    cout << "Введите количество точек для первого сплайна: \n";
    cin >> N1;

    create_vectors_1();

    // Получение координат первого сплайна
    readPoints(xMas1, yMas1, N1);

    if (get_spline_coefficients(xMas1, yMas1, h1, l1, lambda1, delta1, b1, c1, d1, N1) == 0) {
        return 0;
    }

    cout << endl << endl << endl;


    cout << "Введите количество точек для второго сплайна: \n";
    cin >> N2;

    create_vectors_2();

    // Получение координат второго сплайна
    readPoints(xMas2, yMas2, N2);

    if (get_spline_coefficients(xMas2, yMas2, h2, l2, lambda2, delta2, b2, c2, d2, N2) == 0) {
        return 0;
    }
    cout << "Splines are constructed, wait." << endl;
    lower_b = max(xMas1[0], xMas2[0]);
    upper_b = min(xMas1[N1 - 1], xMas2[N2 - 1]);
//    cout << F1(0) << " " << F2(0) << "\n";
    find_min(20);
    return 0;
}