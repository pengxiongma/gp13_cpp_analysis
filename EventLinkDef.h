/*************************************************************************
    > File Name: EventLinkDef.h
    > Author: Ma PengXiong
    > Mail: 1182905413@qq.com 
    > Created Time: Thu 01 Sep 2022 04:19:38 PM CST
 ************************************************************************/

#include<vector>
using namespace std;
#ifdef __CINT__
#pragma link C++ class std::vector<int>;
#pragma link C++ class std::vector<unsigned short>;
#pragma link C++ class std::vector<unsigned int>;
#pragma link C++ class std::vector< std::vector<float> >;
#pragma link C++ class std::vector< std::vector<short> >;
#pragma link C++ class std::vector< std::vector<unsigned short> >;
#endif
