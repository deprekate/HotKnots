default:
	pip install ../hotknots/ --user

clean:
	rm -fr build
	rm -fr dist
	rm -fr hotknots.egg-info
	pip uninstall hotknots -y

