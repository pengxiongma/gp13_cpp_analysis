# gp13_cpp_analysis

# use conda to install ROOT, https://root.cern/install/#conda

#Just call once to make dictionary
rootcling -f EventDict.cxx -c EventLinkDef.h


#g++ -w -o read.out _Read_file.cpp  EventDict.cxx    $(root-config --cflags --libs)   -lMinuit
g++ -w -o exe.out  _Analysis_gp13_ByPengxiong.cpp  EventDict.cxx    $(root-config --cflags --libs)   -lMinuit
./exe.out /home/mapx/mapx/DunhuangData/ROOTFile/TD/Calibration_20dB_XY_GRAND.TEST-RAW-10s-ch2andch3-20dB-du78_81_82.20230417200013.013.dat.root DH_20dB_bp_GausFit_galacticnoise_Calibration_20dB_XY_GRAND.TEST-RAW-10s-ch2andch3-20dB-du78_81_82.20230417200013.013.dat.root_1078_sigma0.root 1 1078 20230417200013

# below is an example to generate shell script to run data, to change to your style and definition
for x in {1078,1081,1082};do ls -1 /home/mapx/mapx/DunhuangData/ROOTFile/TD/Calibration_20dB_XY_GRAND.TEST-RAW-10s-ch2andch3-20dB-du78_81_82.2023*.dat.root  | while read file ; do time=$(echo ${file} | gawk -F "." '{print $3}');  echo -e "${PWD}/exe.out ${file} DH_20dB_bp_GausFit_galacticnoise_$(basename ${file})_${x}_sigma0.root 1 $x ${time}  > ${PWD}/Parameters_GP13/DH_20dB_after20Hz_bp_GausFit_galacticnoise_${x}_sigma0_${time}.LOG  " > _Run_20dB_DU${x}_${time}.sh ;done ;done
