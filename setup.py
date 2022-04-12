#!/usr/bin/env python3
# encoding: utf-8

import os
#from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils import sysconfig


#os.environ["CC"] = "g++"

sources = "temp.cpp Bands.cpp Input.cpp Loop.cpp LoopList.cpp Stack.cpp common.cpp commonPK.cpp goodStem.cpp hotspot.cpp init.cpp initPK.cpp params.cpp paramsPK.cpp python.cpp s_energy_matrix.cpp s_hairpin_loop.cpp s_internal_loop.cpp s_min_folding.cpp s_multi_loop.cpp s_multi_loop_sub.cpp s_partition_function.cpp s_specific_functions.cpp s_stacked_pair.cpp s_sub_folding.cpp sc.cpp score.cpp timer.cpp utils.cpp"

module = Extension('hotknots.hotknots',
			language='c++',
			extra_compile_args=['-g'],
			extra_link_args=['-lm'],
			include_dirs=[
						 '.',
						 '...',
						 os.path.join(os.getcwd(), 'include'),
			],
			sources = ["./src/"+item for item in sources.split()]
			)

def readme():
	with open("README.md", "r") as fh:
		long_desc = fh.read()
	return long_desc

def get_version():
	with open("VERSION.md", 'r') as f:
		v = f.readline().strip()
		return v

def main():
	setup (
		name = 'hotknots',
		version = get_version(),
		author = "Katelyn McNair, Jihong Ren, Baharak Rastegari, Cristina Pop, Mirela Andronescu ",
		author_email = "deprekate@gmail.com",
		description = 'A a tool to predict RNA secondary structure, including pseudoknots',
		long_description = readme(),
		long_description_content_type="text/markdown",
		url =  "https://github.com/deprekate/HotKnots",
		scripts=['hotknots.py'],
		classifiers=[
			"Programming Language :: Python :: 3",
			"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
			"Operating System :: OS Independent",
		],
		python_requires='>3.5.2',
		packages=find_packages(),
		package_data = {'': ['*.dat','*.dh','*.conf','*.dg','*.txt'],},
		#install_requires=[''],
		ext_modules = [module],
		#include_package_data = True,
	)


if __name__ == "__main__":
	main()
