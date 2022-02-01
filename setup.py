#!/usr/bin/env python3
# encoding: utf-8

import os
#from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils import sysconfig

class custom_build_ext(build_ext):
	def build_extensions(self):
		# Override the compiler executables. Importantly, this
		# removes the "default" compiler flags that would
		# otherwise get passed on to to the compiler, i.e.,
		# distutils.sysconfig.get_var("CFLAGS").
		#print(distutils.sysconfig._config_vars)
		compiler = sysconfig.get_config_vars("CC")
		self.compiler.set_executable("compiler_so", compiler)
		self.compiler.set_executable("compiler_cxx", compiler)
		self.compiler.set_executable("linker_so", compiler)
		build_ext.build_extensions(self)

#os.environ["CC"] = "gcc"
compile_args = ["-g"] # -Wall -O2"]
link_args	= ["-lsimfold", "-lLEModel"]
#deps = "src/LinearFoldEval.cpp src/LinearFold.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h".split(' ') 

sources = "HotKnot.cpp HotKnotEnergy.cpp common.cpp computeEnergy.cpp goodStem.cpp hotspot.cpp init.cpp params.cpp python.cpp s_energy_matrix.cpp s_hairpin_loop.cpp s_internal_loop.cpp s_min_folding.cpp s_multi_loop.cpp s_multi_loop_sub.cpp s_partition_function.cpp s_specific_functions.cpp s_stacked_pair.cpp s_sub_folding.cpp sc.cpp score.cpp timer.cpp utils.cpp"

mod= Extension('simfold',
			language='c++',
			extra_compile_args=['-g', '-Wunused-variable'],
			extra_link_args=['-lm'],
			include_dirs=[
						 '.',
						 '...',
						 os.path.join(os.getcwd(), 'include'),
			],
			sources = ["./src/"+item for item in sources.split()]
			)
module = Extension('HotKnots',
			language='c++',
			extra_compile_args=compile_args,
			extra_link_args = ["-lLEModel", "-lsimfold"],
			include_dirs=[
						 '.',
						 '...',
						 os.path.join(os.getcwd(), 'src'),
						 os.path.join(os.getcwd(), 'simfold'),
						 os.path.join(os.getcwd(), 'LE'),
			],
			library_dirs = [os.getcwd(),'./LE', './simfold'],
			sources = ['src/python.cpp'] ) #, 'src/LinearFold.cpp']  )

mod_sim = Extension('simfold',
			language='c++',
			extra_compile_args=['-g','-Wno-deprecated', '-Wunused-variable'],
			extra_link_args=['-lm'],
			include_dirs=[
						 '.',
						 '...',
						 os.path.join(os.getcwd(), 'LE/H'),
						 os.path.join(os.getcwd(), 'simfold/include'),
						 os.path.join(os.getcwd(), 'simfold/src/common'),
						 os.path.join(os.getcwd(), 'simfold/src/simfold'),
			],
			library_dirs = [os.getcwd(),'./LE', './simfold'],
			sources = ["./LE/"+item for item in "Loop.cpp Stack.cpp Input.cpp Bands.cpp LoopList.cpp commonPK.cpp paramsPK.cpp initPK.cpp".split()]
			)
mod_le = Extension('LE',
			language='c++',
			extra_compile_args=['-g','-Wno-deprecated', '-Wunused-variable'],
			extra_link_args=['-lm'],
			include_dirs=[
						 '.',
						 '...',
						 os.path.join(os.getcwd(), 'LE/H'),
						 os.path.join(os.getcwd(), 'simfold/include'),
						 os.path.join(os.getcwd(), 'simfold/src/common'),
						 os.path.join(os.getcwd(), 'simfold/src/simfold'),
			],
			library_dirs = [os.getcwd(),'./LE', './simfold'],
			sources = ["./LE/"+item for item in "Loop.cpp Stack.cpp Input.cpp Bands.cpp LoopList.cpp commonPK.cpp paramsPK.cpp initPK.cpp".split()]
			)

def readme():
	with open("README.md", "r") as fh:
		long_desc = fh.read()
	return long_desc

def get_version():
	with open("VERSION", 'r') as f:
		v = f.readline().strip()
		return v

def main():
	setup (
		name = 'HotKnots',
		version = get_version(),
		author = "Katelyn McNair, Jihong Ren, Baharak Rastegari, Cristina Pop, Mirela Andronescu ",
		author_email = "deprekate@gmail.com",
		description = 'A a tool to predict RNA secondary structure, including pseudoknots',
		long_description = readme(),
		long_description_content_type="text/markdown",
		url =  "https://github.com/deprekate/HotKnots",
		scripts=[],
		classifiers=[
			"Programming Language :: Python :: 3",
			"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
			"Operating System :: OS Independent",
		],
		python_requires='>3.5.2',
		packages=find_packages(),
		#install_requires=[''],
		ext_modules = [mod],
		#cmdclass={"build_ext":custom_build_ext}
		#cmdclass={"build_ext":build_ext}
	)


if __name__ == "__main__":
	main()
