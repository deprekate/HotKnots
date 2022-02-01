default:

	# first make the 2 libraries
	cd simfold; make libsimfold.a
	cd LE; make
	
	# make the HotKnots and computeEnergy
	cd src; make
	
	# copy them into the bin directory
	mv src/HotKnots .

clean:
	cd simfold; make clean
	cd LE; make clean
	cd src; make clean
	rm -rf HotKnots
	rm -fr build
	rm -fr dist
