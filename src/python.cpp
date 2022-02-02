//#include <stdio.h>
//#include <limits.h>
#include <Python.h>

extern int run( const char *string, char *currentModel, char *paramsFile);

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

static PyObject* fold(PyObject *self, PyObject *args, PyObject *kwargs){
	char *string;
	char *model;
	char *params;
	static char *kwlist[] = {(char *) "string", (char *) "model", (char *) "params", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "s|ss", kwlist, &string, &model, &params)){
		return NULL;
	}
	run("UUUGCCCUGAAACUGGCGCGUGAGAUGGGGCGACCCGACUGGCGUGCCAU", "CC", "params/parameters_DP09.txt");	
	return Py_BuildValue("[sf]", "string", -100.0);
}

// Module method definitions
static PyObject* no_args(PyObject *self, PyObject *args) {
	Py_RETURN_NONE;
}

// Method definition object for this extension, these argumens mean:
static PyMethodDef HotKnots_methods[] = { 
	{"get_windows",         get_windows, METH_VARARGS,                 "Empty for now, can be used to yield a python iterator."},  
	{"fold",        (PyCFunction)  fold, METH_VARARGS | METH_KEYWORDS, "Calculates the minimum free energy of the sequence."},  
	{"no_args",                 no_args, METH_NOARGS,                  "Empty for now."},  
	{NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef HotKnots_definition = { 
	PyModuleDef_HEAD_INIT,
	"HotKnots",
	"A Python module that does RNA folding.",
	-1, 
	HotKnots_methods
};

// Module initialization
// Python calls this function when importing your extension. It is important
// that this function is named PyInit_[[your_module_name]] exactly, and matches
// the name keyword argument in setup.py's setup() call.
PyMODINIT_FUNC PyInit_HotKnots(void) {
	//Py_Initialize();
	return PyModule_Create(&HotKnots_definition);
}





