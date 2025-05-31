/*
 * CNav.cpp
 *
 *  Created on: 3 ���. 2023
 *      Author: User
 */

#include "CNav.hpp"

CNav::CNav(CLog *_log, const char *fn) : log(_log) {
    nav_file = fopen(fn, "r");
    if (!nav_file) {
        printf("Ошибка чтения файла: %s!\n", fn);
    }
}

CNav::CNav() {}

CNav::~CNav() {
    if (nav_file) {
        fclose(nav_file);
        nav_file = NULL;
    }
    nav_vec.clear();
}

int CNav::read_msec_from(uint32_t from_msec, uint32_t msec) {

    if (!nav_file) {
        return -1;
    }

    rewind(nav_file);

    // printf("sizeof(nav_pack_t): %ld\n", sizeof(nav_pack_t) );
    int n = fread(&nav_pack, sizeof(nav_pack_t), 1, nav_file);
    if (!n || nav_pack.marker() != 0x55AA55AA) {
        log->LogPrintf("Отсутствуют навигационные данные!\n");
        return -2;
    }

    const uint64_t start_time = nav_pack.time();
    while ((nav_pack.time() - start_time) < from_msec) {
        n = fread(&nav_pack, sizeof(nav_pack_t), 1, nav_file);

        if (!n) {
            log->LogPrintf("Ошибка: не найдена временная метка начала "
                           "синтезирования ( %d мсек )\n",
                           from_msec);
            return -3;
        }
    }

    nav_vec.clear();
    nav_vec.push_back(nav_pack);

    while ((nav_pack.time() - start_time - from_msec) < msec) {
        n = fread(&nav_pack, sizeof(nav_pack_t), 1, nav_file);
        if (!n) {
            log->LogPrintf("Внимание: конец файла достигнут раньше задуманного "
                           "( %d отсчетов загружено )\n",
                           nav_vec.size());
            return -4;
        }

        nav_vec.push_back(nav_pack);
    }

    return 0;
}

int CNav::read_msec(uint32_t msec) { return read_msec_from(0, msec); }

int CNav::read_sec(float sec) { return read_msec_from(0, sec * 1000); }

int CNav::read_sec_from(float sec_from, float sec) {
    return read_msec_from(sec_from, sec * 1000);
}

size_t CNav::size() { return nav_vec.size(); }

Velocity CNav::velocity() { return Velocity(this); }

Elevation CNav::elevation() { return Elevation(this); }

double CNav::lat(size_t n) {
    if (n < nav_vec.size()) {
        return nav_vec[n].lat();
    }
    return NAN;
}

double CNav::lon(size_t n) {
    if (n < nav_vec.size()) {
        return nav_vec[n].lon();
    }
    return NAN;
}

double CNav::course(size_t n) {
    if (n < nav_vec.size()) {
        return nav_vec[n].course();
    }
    return NAN;
}

std::vector<Dim3> CNav::gyro() {
    std::vector<Dim3> vec;

    for (auto nav : this->nav_vec) {
        vec.push_back(Dim3(nav.gx(), nav.gy(), nav.gz()));
    }

    return vec;
}

std::vector<Dim3> CNav::accel() {
    std::vector<Dim3> vec;

    for (auto nav : this->nav_vec) {
        vec.push_back(Dim3(nav.ax(), nav.ay(), nav.az()));
    }

    return vec;
}

std::vector<Dim3> CNav::compass() {
    std::vector<Dim3> vec;

    for (auto nav : this->nav_vec) {
        vec.push_back(Dim3(nav.cx(), nav.cy(), nav.cz()));
    }

    return vec;
}
