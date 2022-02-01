default:

	# first make the 2 libraries
	cd simfold; make libsimfold.a
	cd LE; make
	
	# make the HotKnots and computeEnergy
	cd src; make
	
	# copy them into the bin directory
	mv src/HotKnots .
	mv src/computeEnergy .

clean:
	cd simfold; make clean
	cd LE; make clean
	cd src; make clean
	rm -rf HotKnots computeEnergy 
	rm -fr build
	rm -fr dist
