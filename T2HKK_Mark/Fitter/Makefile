all:
	make -C ThrowParams/src
	make -C OscParams/src
	make -C SystParams/src
	make -C MakeHistograms/src
	make -C src

links:
	@ [ -d include ] || mkdir include
	@ [ -d lib ]     || mkdir lib
	@ [ -d bin ]     || mkdir bin
	ln -s ../../Prob3++/BargerPropagator.h include/BargerPropagator.h
	ln -s ../../Prob3++/NeutrinoPropagator.h include/NeutrinoPropagator.h
	ln -s ../../Prob3++/EarthDensity.h include/EarthDensity.h
	ln -s ../../Prob3++/mosc3.h include/mosc3.h
	ln -s ../../Prob3++/mosc.h include/mosc.h
	ln -s ../../Prob3++/probWrapper.h include/probWrapper.h
	ln -s ../ThrowParams/src/ThrowParams.h include/ThrowParams.h
	ln -s ../OscParams/src/OscParams.h include/OscParams.h
	ln -s ../SystParams/src/SystParams.h include/SystParams.h
	ln -s ../MakeHistograms/src/MakeHistograms.h include/MakeHistograms.h

clean:
	make clean -C ThrowParams/src
	make clean -C OscParams/src
	make clean -C SystParams/src
	make clean -C MakeHistograms/src
	make clean -C src



