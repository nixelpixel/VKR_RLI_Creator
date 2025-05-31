/*
 * CRadarSTT.h
 *
 *  Created on: 3 июл. 2017 г.
 *      Author: user
 */

#ifndef CRADARSTT_H_
#define CRADARSTT_H_

#include "CLog.hpp"
#include "CRadar.hpp"
#include "CTCPServer.hpp"

class CRadar_STT : public CRadar {
  public:
    CRadar_STT(CLog *log, const char *address, const char *port);
    CRadar_STT(CLog *log);
    CRadar_STT();
    virtual ~CRadar_STT();

    int ReadDataFromFile(char *filename, CImage *image, int Np, float tshift,
                         int syncpos, int kr, int comp2, int rmdist);
    int ReadBlockFromFile(char *filename, CImage *image, int np, int kr);
    int ReadBlockSetPos(int syncpos, float tshift);
    int ReadBlockResetPos(void);
    // int               ReadDataFromSocket(void);
    int FindSyncPos(CImage *image, int *sync_begin, int *sync_width);
    pybind11::tuple FindSyncPosPy(CImage *image);
    float GetFileDuration(char *filename);

    int RemoveTransients(CImage *image);

    int OpenHologram(char *fn, CImage *image, int np, int syncpos, int kr,
                     int comp2);
    int ReadBlock(CImage *image, int block);
    void FreeHologram();

    int OpenHologramSocket(char *ip_addr, CImage *image, int np, int syncpos,
                           int kr, int comp2);
    int ReadBlockFromSocket(CImage *image);

    int InitDataFile(char *_fn, CImage *image, int np, int kr, int comp2);
    int ReadBlockFromInitializedFile(CImage *image, int syncpos, int block);

  private:
    CTCPServer *tcp_srv = nullptr;
};

#endif /* CRADARSTT_H_ */
