//#include <stdio.h>
//#include <limits.h>
#include <Python.h>

#include "hotspot.h"
#include "utils.h"

extern int initiate(char *currentModel, char *paramsFile, char *cfile, char *cfilepk);
extern struct Fold* best( char *sequence, char *currentModel);

typedef struct {
	PyObject_HEAD
	const unsigned char* dna;
	unsigned int len;
	unsigned int i;
	unsigned int f;
} windows_Iterator;

PyObject* windows_Iterator_iter(PyObject *self){
	Py_INCREF(self);
	return self;
}

PyObject* windows_Iterator_iternext(PyObject *self){
	windows_Iterator *p = (windows_Iterator *)self;
	if( (p->i)  <  (p->len - 2) ){
		PyObject *py_list = Py_BuildValue("[f]", -1.0);
		return py_list;
	}else{
		PyErr_SetNone(PyExc_StopIteration);
		return NULL;
	}
}

static void Iter_dealloc(windows_Iterator *self){ PyObject_Del(self); }

static PyTypeObject IterableType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "Iter",
	.tp_basicsize = sizeof(windows_Iterator),
	.tp_itemsize = 0,
	.tp_dealloc = (destructor) Iter_dealloc,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	.tp_doc = "Custom objects",
	.tp_iter	  = windows_Iterator_iter,
	.tp_iternext  = windows_Iterator_iternext
};

static PyObject * get_windows(PyObject *self, PyObject *args){
	windows_Iterator *p;
	p = PyObject_New(windows_Iterator, &IterableType);
	if (!p) return NULL;
	if (PyType_Ready(&IterableType) < 0) {
		return NULL;
	}
	if (!PyArg_ParseTuple(args, "s", &p->dna)) {
		return NULL;
	}
	
	p->i = 0;
	p->f = 1;
	p->len = strlen( (const char*) p->dna);
	for (int i=0; p->dna[i] ; i++){
	}

	/* I'm not sure if it's strictly necessary. */
	if (!PyObject_Init((PyObject *)p, &IterableType)) {
		Py_DECREF(p);
		return NULL;
	}

	return (PyObject *)p;
}

// Module method definitions
static PyObject* initialize(PyObject *self, PyObject *args, PyObject *kwargs){
	char *model;
	char *params;
	char *cfile;
	char *cfilepk;
	static char *kwlist[] = { (char *) "model", (char *) "params", (char *) "cfile", (char *) "cfilepk",  NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "ssss", kwlist, &model, &params, &cfile, &cfilepk)){
		return NULL;
	}
	initiate(model , params, cfile, cfilepk);	
	Py_RETURN_NONE;
}

static PyObject* fold(PyObject *self, PyObject *args, PyObject *kwargs){
	char *sequence;
	char *model;
	PyObject *retval;
	static char *kwlist[] = { (char *) "sequence", (char *) "model", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s|s", kwlist, &sequence, &model)){
		return NULL;
	}
	//struct Node* bestNode = best( sequence , model);	
	//bpseq2dp( (int) strlen(sequence), bestNode->secStructure, outputStructure);
	struct Fold* bestFold = best( sequence , model);	
	retval = Py_BuildValue("[sf]", bestFold->structure, bestFold->score);
	free(bestFold);
	return retval;
}

// Method definition object for this extension, these argumens mean:
static PyMethodDef hotknots_methods[] = { 
	{"get_windows",         get_windows,            METH_VARARGS,                 "Empty for now, can be used to yield a python iterator."},  
	{"initialize",        (PyCFunction) initialize, METH_VARARGS | METH_KEYWORDS, "Calculates the minimum free energy of the sequence."},  
	{"fold",              (PyCFunction) fold,       METH_VARARGS | METH_KEYWORDS,                  "do it"},  
	{NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef hotknots_definition = { 
	PyModuleDef_HEAD_INIT,
	"HotKnots",
	"A Python module that does RNA folding.",
	-1, 
	hotknots_methods
};

// Module initialization
// Python calls this function when importing your extension. It is important
// that this function is named PyInit_[[your_module_name]] exactly, and matches
// the name keyword argument in setup.py's setup() call.
PyMODINIT_FUNC PyInit_hotknots(void) {
	//Py_Initialize();
	return PyModule_Create(&hotknots_definition);
}





