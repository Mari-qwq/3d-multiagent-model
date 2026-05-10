#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

struct Segment {
    double t0, t1;

    double x0, y0, z0;

    double vx, vy, vz;
};

struct Agent {
    vector<Segment> segs;
};

static mt19937_64 rng(1234567);

double rnd(double a, double b) {
    uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

Agent generate_agent(
    double a, double b, double c,
    double Vmax,
    double Tmax,
    double T
) {
    Agent A;

    double t = 0.0;

    double x = rnd(0.0, a);
    double y = rnd(0.0, b);
    double z = rnd(0.0, c);

    const double EPS = 1e-12;

    while (t < Tmax - EPS) {

        double dt = rnd(0.0, T);

        if (dt <= EPS)
            dt = EPS;

        // Сферические углы
        double theta = rnd(0.0, 2.0 * M_PI);
        double phi = rnd(-M_PI / 2.0, M_PI / 2.0);

        double v = rnd(0.0, Vmax);

        // Компоненты скорости
        double vx = v * cos(phi) * cos(theta);
        double vy = v * cos(phi) * sin(theta);
        double vz = v * sin(phi);

        double rem = min(dt, Tmax - t);

        while (rem > EPS) {

            // Время до столкновения со стенками
            double tx = INFINITY;
            double ty = INFINITY;
            double tz = INFINITY;

            if (vx > EPS)
                tx = (a - x) / vx;
            else if (vx < -EPS)
                tx = (0.0 - x) / vx;

            if (vy > EPS)
                ty = (b - y) / vy;
            else if (vy < -EPS)
                ty = (0.0 - y) / vy;

            if (vz > EPS)
                tz = (c - z) / vz;
            else if (vz < -EPS)
                tz = (0.0 - z) / vz;

            double thit = min(tx, min(ty, tz));

            // Нет столкновения
            if (!isfinite(thit) || thit > rem - EPS) {

                A.segs.push_back({
                    t,
                    t + rem,
                    x, y, z,
                    vx, vy, vz
                    });

                x += vx * rem;
                y += vy * rem;
                z += vz * rem;

                x = max(0.0, min(a, x));
                y = max(0.0, min(b, y));
                z = max(0.0, min(c, z));

                t += rem;
                rem = 0.0;
            }
            else {

                double dt1 = max(0.0, thit);

                A.segs.push_back({
                    t,
                    t + dt1,
                    x, y, z,
                    vx, vy, vz
                    });

                x += vx * dt1;
                y += vy * dt1;
                z += vz * dt1;

                bool hitX = fabs(thit - tx) < 1e-12;
                bool hitY = fabs(thit - ty) < 1e-12;
                bool hitZ = fabs(thit - tz) < 1e-12;

                // Отражения
                if (hitX) {
                    x = (vx > 0 ? a : 0.0);
                    vx = -vx;
                }

                if (hitY) {
                    y = (vy > 0 ? b : 0.0);
                    vy = -vy;
                }

                if (hitZ) {
                    z = (vz > 0 ? c : 0.0);
                    vz = -vz;
                }

                t += dt1;
                rem -= dt1;
            }
        }
    }

    return A;
}

struct Point3D {
    double x, y, z;
};

Point3D position(const Agent& A, double t) {

    for (const auto& s : A.segs) {

        if (t >= s.t0 && t <= s.t1) {

            double u = t - s.t0;

            return {
                s.x0 + s.vx * u,
                s.y0 + s.vy * u,
                s.z0 + s.vz * u
            };
        }
    }

    const auto& s = A.segs.back();

    double u = s.t1 - s.t0;

    return {
        s.x0 + s.vx * u,
        s.y0 + s.vy * u,
        s.z0 + s.vz * u
    };
}

vector<double> events_for_pair(
    const Agent& A,
    const Agent& B,
    double r,
    double Tmax
) {
    vector<double> ev;

    size_t ia = 0;
    size_t ib = 0;

    double r2 = r * r;

    const double EPS = 1e-12;

    while (ia < A.segs.size() && ib < B.segs.size()) {

        const auto& sa = A.segs[ia];
        const auto& sb = B.segs[ib];

        double t0 = max(sa.t0, sb.t0);
        double t1 = min(sa.t1, sb.t1);

        if (t1 > Tmax)
            t1 = Tmax;

        if (t0 < t1 - EPS) {

            double xi = sa.x0 + sa.vx * (t0 - sa.t0);
            double yi = sa.y0 + sa.vy * (t0 - sa.t0);
            double zi = sa.z0 + sa.vz * (t0 - sa.t0);

            double xj = sb.x0 + sb.vx * (t0 - sb.t0);
            double yj = sb.y0 + sb.vy * (t0 - sb.t0);
            double zj = sb.z0 + sb.vz * (t0 - sb.t0);

            double rx = xi - xj;
            double ry = yi - yj;
            double rz = zi - zj;

            double ux = sa.vx - sb.vx;
            double uy = sa.vy - sb.vy;
            double uz = sa.vz - sb.vz;

            double alpha =
                ux * ux +
                uy * uy +
                uz * uz;

            double beta =
                2.0 * (
                    rx * ux +
                    ry * uy +
                    rz * uz
                    );

            double gamma =
                rx * rx +
                ry * ry +
                rz * rz -
                r2;

            double L = t1 - t0;

            if (alpha > EPS) {

                double D =
                    beta * beta -
                    4.0 * alpha * gamma;

                if (D >= -EPS) {

                    D = max(0.0, D);

                    double root = sqrt(D);

                    double s1 =
                        (-beta - root) /
                        (2.0 * alpha);

                    double s2 =
                        (-beta + root) /
                        (2.0 * alpha);

                    if (s1 >= -EPS && s1 <= L + EPS)
                        ev.push_back(t0 + s1);

                    if (s2 >= -EPS && s2 <= L + EPS)
                        ev.push_back(t0 + s2);
                }
            }
        }

        if (sa.t1 < sb.t1 - EPS)
            ia++;
        else if (sb.t1 < sa.t1 - EPS)
            ib++;
        else {
            ia++;
            ib++;
        }
    }

    sort(ev.begin(), ev.end());

    vector<double> clean;

    for (double t : ev) {

        if (
            clean.empty() ||
            fabs(clean.back() - t) > 1e-9
            )
        {
            clean.push_back(t);
        }
    }

    return clean;
}

struct PairDurations {
    vector<double> inside;
    vector<double> outside;
};

PairDurations durations_for_pair(
    const Agent& A,
    const Agent& B,
    double r,
    double Tmax
) {
    PairDurations pd;

    auto ev = events_for_pair(A, B, r, Tmax);

    double r2 = r * r;

    double prev = 0.0;

    auto pA0 = position(A, 0.0);
    auto pB0 = position(B, 0.0);

    double dx = pA0.x - pB0.x;
    double dy = pA0.y - pB0.y;
    double dz = pA0.z - pB0.z;

    bool inside =
        (dx * dx + dy * dy + dz * dz)
        <= r2 + 1e-12;

    for (double t : ev) {

        if (t > Tmax)
            break;

        double dt = t - prev;

        if (dt > 1e-12) {

            if (inside)
                pd.inside.push_back(dt);
            else
                pd.outside.push_back(dt);
        }

        inside = !inside;

        prev = t;
    }

    if (prev < Tmax - 1e-12) {

        double dt = Tmax - prev;

        if (inside)
            pd.inside.push_back(dt);
        else
            pd.outside.push_back(dt);
    }

    return pd;
}

double mean_of(const vector<double>& v) {

    if (v.empty())
        return 0.0;

    long double s = 0.0;

    for (double x : v)
        s += x;

    return (double)(s / v.size());
}

int main() {

    // Размеры параллелепипеда
    double a = 30.0;
    double b = 30.0;
    double c = 30.0;

    double r = 0.5;

    int n = 5;

    double Tmax = 200.0;

    int K = 10000;

    // Более плотная сетка параметров
    vector<double> Vmax_grid = {
        0.5, 1.0, 1.5,
        2.0, 2.5, 3.0,
        3.5, 4.0, 4.5, 5.0
    };

    vector<double> T_grid = {
        0.5, 1.0, 1.5,
        2.0, 2.5, 3.0,
        3.5, 4.0, 4.5, 5.0
    };

    ofstream summary("summary3D.csv");

    summary <<
        "Vmax,T,"
        "mean_inside,"
        "mean_outside,"
        "count_inside,"
        "count_outside\n";

    for (double Vmax : Vmax_grid) {

        for (double T : T_grid) {

            vector<double> all_inside;
            vector<double> all_outside;

            for (int k = 0; k < K; ++k) {

                vector<Agent> agents;

                agents.reserve(n);

                for (int i = 0; i < n; ++i) {

                    agents.push_back(
                        generate_agent(
                            a, b, c,
                            Vmax,
                            Tmax,
                            T
                        )
                    );
                }

                for (int i = 0; i < n; ++i) {

                    for (int j = i + 1; j < n; ++j) {

                        auto pd =
                            durations_for_pair(
                                agents[i],
                                agents[j],
                                r,
                                Tmax
                            );

                        all_inside.insert(
                            all_inside.end(),
                            pd.inside.begin(),
                            pd.inside.end()
                        );

                        all_outside.insert(
                            all_outside.end(),
                            pd.outside.begin(),
                            pd.outside.end()
                        );
                    }
                }
            }

            double meanInside =
                mean_of(all_inside);

            double meanOutside =
                mean_of(all_outside);

            summary
                << Vmax << ","
                << T << ","
                << meanInside << ","
                << meanOutside << ","
                << all_inside.size() << ","
                << all_outside.size()
                << "\n";

            // Сохранение сырых данных inside/outside

            string fname_in =
                "inside_Vmax_" +
                to_string(Vmax) +
                "_T_" +
                to_string(T) +
                ".csv";

            string fname_out =
                "outside_Vmax_" +
                to_string(Vmax) +
                "_T_" +
                to_string(T) +
                ".csv";

            ofstream fin(fname_in);
            ofstream fout(fname_out);

            for (double x : all_inside)
                fin << x << "\n";

            for (double x : all_outside)
                fout << x << "\n";

            fin.close();
            fout.close();

            cout
                << "Done Vmax="
                << Vmax
                << " T="
                << T
                << endl;
        }
    }

    summary.close();

    cout << "3D simulation completed.\n";

    return 0;
}