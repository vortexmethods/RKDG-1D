#ifndef FLUXHLLC_H_
#define FLUXHLLC_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "Flux.h"
#include "ProblemGas1D.h"

using namespace std;

#define HLLC FluxHLLC

class FluxHLLC :
    public Flux
{
private:

    //РЎРѕР±СЃС‚РІРµРЅРЅС‹Рµ С‡РёСЃР»Р° РЅР° РІСЃРµР№ СЃРµС‚РєРµ (РјРµР¶РґСѓ i-Р№ Рё (i-1)-Р№ СЏС‡РµР№РєР°РјРё) 
    vector<vector<double>> L;

	//РЎРїРѕСЃРѕР± РІС‹С‡РёСЃР»РµРЅРёСЏ СЃРєРѕСЂРѕСЃС‚Рё РЅР° РіСЂР°РЅРёС†Рµ СЏС‡РµРµРє
	SoundVelType SoundVel;

    //РџРѕС‚РѕРєРё HLLРЎ РЅР° Р»РµРІРѕРј Рё РїСЂР°РІРѕРј РєРѕРЅС†Рµ СЏС‡РµР№РєРё
    vector<double> hllcL, hllcR;

    //HLLC-СЂРµС€РµРЅРёСЏ "СЃРѕ Р·РІРµР·РґРѕС‡РєРѕР№"
    vector<double> UastL, UastR;

    //Р’СЂРµРјРµРЅРЅС‹Рµ РїРµСЂРµРјРµРЅРЅС‹Рµ: РІРµРєС‚РѕСЂС‹ - HLLC-СЃРєР°С‡РѕРє СЂРµС€РµРЅРёСЏ РЅР° РіСЂР°РЅРёС†Р°С… i-Р№ СЏС‡РµР№РєРё
    vector<double> LdUast, RdUast;

    void getUast(const int cell, const double Sast, \
		const vector<vector<double>>& leftsol, \
		const vector<vector<double>>& mysol, \
        vector<double>& Uast);

    double getSast(const int cell, \
		const vector<vector<double>>& leftsol, \
		const vector<vector<double>>& mysol);
public:

    //РљРѕРЅСЃС‚СЂСѓРєС‚РѕСЂ  (prm - РѕР±СЉРµРєС‚, СЃРѕРґРµСЂР¶Р°С‰РёР№ РїР°СЂР°РјРµС‚СЂС‹ Р·Р°РґР°С‡Рё)
	FluxHLLC(const BaseParams& prm, Problem& prb, SoundVelType soundvel);

    //Р”РµСЃС‚СЂСѓРєС‚РѕСЂ
    ~FluxHLLC();

	void step(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& DSOL, \
		const double cft);

};

#endif