default:

	# first make the 2 libraries
	cd simfold; make libsimfold.a
	cd LE; make
	
	# make the HotKnots and computeEnergy
	cd hotspot; make
	
	# copy them into the bin directory
	cp hotspot/HotKnots bin
	cp hotspot/computeEnergy bin

clean:
	cd simfold; make clean
	cd LE; make clean
	cd hotspot; make clean
	rm -rf bin/HotKnots bin/computeEnergy 
	rm -rf bin/output/*

