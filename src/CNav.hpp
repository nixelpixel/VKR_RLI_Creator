/*
 * CNav.h
 *
 *  Created on: 3 ���. 2023
 *      Author: User
 */

#ifndef CNav_H_
#define CNav_H_

#include "CLog.hpp"
#include "helpers.hpp"
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <tuple>
#include <vector>

class Velocity;
class Elevation;
class Dim3;
class CNav {
  public:
    CNav();
    CNav(CLog *_log, const char *fn);
    virtual ~CNav();
    int read_msec_from(uint32_t from_msec, uint32_t msec);
    int read_msec(uint32_t msec);
    int read_sec_from(float sec_from, float sec);
    int read_sec(float sec);
    size_t size();
    Velocity velocity();
    Elevation elevation();

    double lat(size_t n = 0);
    double lon(size_t n = 0);
    double course(size_t n = 0);

    std::vector<Dim3> gyro();
    std::vector<Dim3> accel();
    std::vector<Dim3> compass();

    friend class Velocity;
    friend class Elevation;
    friend class ParameterBase;

  private:
    struct nav_pack_t {
        uint32_t _marker = 0;
        uint8_t version = 0;
        double _lat = 0;
        double _lon = 0;
        double _ele = 0;
        double _velocity_course = 0;
        float _pressure = 0;
        int16_t _temperature = 0;
        int16_t _gx = 0;
        int16_t _gy = 0;
        int16_t _gz = 0;
        int16_t _ax = 0;
        int16_t _ay = 0;
        int16_t _az = 0;
        int16_t _cx = 0;
        int16_t _cy = 0;
        int16_t _cz = 0;
        double pitch = 0;
        double roll = 0;
        double yaw = 0;
        double _course = 0;
        uint64_t _time = 0;
        uint8_t valid = 0;
        uint16_t crc = 0;

        float pressure() { return swap_endian(_pressure); }
        float temperature() {
            float t = swap_endian(_temperature);
            return t / 10.0;
        }

        int16_t gx() { return swap_endian(_gx); }
        int16_t gy() { return swap_endian(_gy); }
        int16_t gz() { return swap_endian(_gz); }
        int16_t ax() { return swap_endian(_ax); }
        int16_t ay() { return swap_endian(_ay); }
        int16_t az() { return swap_endian(_az); }
        int16_t cx() { return swap_endian(_cx); }
        int16_t cy() { return swap_endian(_cy); }
        int16_t cz() { return swap_endian(_cz); }

        double lat() {
            // return rad2deg(swap_endian<double>(_lat));
            return swap_endian<double>(_lat);
        }

        double lon() {
            return swap_endian<double>(_lon);
            // return rad2deg(swap_endian<double>(_lon));
        }

        double course() {
            return swap_endian<double>(_course);
            // return rad2deg(swap_endian<double>(_course));
        }

        double velocity_course() {
            return swap_endian<double>(_velocity_course);
        }

        double ele() { return swap_endian<double>(_ele); }
        uint32_t marker() { return swap_endian<uint32_t>(_marker); }
        uint64_t time() { return swap_endian<uint64_t>(_time); }

    } __attribute__((packed));

    struct nav_pack_t nav_pack;
    FILE *nav_file;
    CLog *log;
    std::vector<struct nav_pack_t> nav_vec;
};

class ParameterBase {
  public:
    ParameterBase(CNav *_nav) : Nav(_nav) {}

  protected:
    CNav *Nav;
};

class Velocity : public ParameterBase {
  public:
    Velocity(CNav *_nav) : ParameterBase(_nav) {}

    float mean() { return (this->*_mean)(); }
    float get(int n) { return (this->*_get)(n); }

    Velocity as_kmh() {
        _mean = &Velocity::mean_kmh;
        _get = &Velocity::get_kmh;
        return *this;
    }

    Velocity as_si() {
        _mean = &Velocity::mean_si;
        _get = &Velocity::get_si;
        return *this;
    }

  private:
    float (Velocity::*_mean)() = &Velocity::mean_si;
    float (Velocity::*_get)(int n) = &Velocity::get_si;

    float mean_si() {
        float acc = 0;
        for (auto nav : Nav->nav_vec) {
            acc += nav.velocity_course();
        }
        acc /= Nav->nav_vec.size();
        return acc;
    }

    float mean_kmh() { return mean_si() * 3.6; }

    float get_si(int n) {
        if (n < Nav->nav_vec.size()) {
            return Nav->nav_vec[n].velocity_course();
        }
        return NAN;
    }

    float get_kmh(int n) { return get_si(n) * 3.6; }
};

class Elevation : public ParameterBase {
  public:
    Elevation(CNav *_nav) : ParameterBase(_nav) {}
    float mean() {
        float acc = 0;
        for (auto nav : Nav->nav_vec) {
            acc += nav.ele();
        }
        acc /= Nav->nav_vec.size();
        return acc;
    }
    float mean_above(float above) { return mean() - above; }
};

class Dim3 {
  public:
    Dim3() {}
    Dim3(int16_t _x, int16_t _y, int16_t _z) {
        x = _x;
        y = _y;
        z = _z;
    }
    int16_t x = 0;
    int16_t y = 0;
    int16_t z = 0;
};

#endif /* CNav_H_ */
