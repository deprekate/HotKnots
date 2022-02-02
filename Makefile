default:
	python3 setup.py install --user

clean:
	rm -fr build
	rm -fr dist
	rm -fr HotKnots.egg-info
	pip uninstall HotKnots -y

