#include <Python.h>
#include "wrapper.h"

/* Docstrings */
static char module_docstring[] =
    "This library is a wrapper ";

/* Available functions */
static PyObject *_Init_ADC(PyObject *self, PyObject *args);
static PyObject *_Read_Single_Channel(PyObject *self, PyObject *args);
static PyObject *_Init_Single_Channel(PyObject *self, PyObject *args);
static PyObject *_ADC_Stop(PyObject *self, PyObject *args);
static PyObject *_thread1(PyObject *self, PyObject *args);
static PyObject *_kill_flag(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
 //   {"chi2", chi2_chi2, METH_VARARGS, chi2_docstring},
    {"init", _Init_ADC, METH_VARARGS, {"Initializes the ADC"}},
    {"read_channel", _Read_Single_Channel, METH_VARARGS, {"Reads a given channel"}},
    {"init_channel", _Init_Single_Channel, METH_VARARGS, {"Initializes a given channel"}},
    {"stop", _ADC_Stop, 0, {"Closes the ADC"}},
    {"run", _thread1, METH_VARARGS, {"Runs the main program"}},
    {"kill", _kill_flag, METH_VARARGS, {"Softly kills the program"}},
    {NULL, NULL, 0, NULL}
};

/* Initialize the module */
//~ PyMODINIT_FUNC initads1256(void)
//~ {
    //~ PyObject *m = Py_InitModule3("ads1256", module_methods, module_docstring);
    //~ if (m == NULL)
        //~ return;

//~ }

static struct PyModuleDef initads1256 =
{
    PyModuleDef_HEAD_INIT,
    "ads1256", /* name of module */
    module_docstring, /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};

PyMODINIT_FUNC PyInit_ads1256(void)
{
    return PyModule_Create(&initads1256);
}

static PyObject *_Init_ADC(PyObject *self, PyObject *args)
{
    double gain, sps;
    unsigned char mode;
    PyObject *yerr_obj;
    
    printf("ERROR");
    
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args,"ddb",&gain,&sps,&mode,&yerr_obj))
        return NULL;

    /* execute the code */ 
    Init_ADC(gain,sps,mode);
    
    Py_RETURN_NONE;
}

static PyObject *_Read_Single_Channel(PyObject *self, PyObject *args)
{
    unsigned char  ch;
    int value;
    PyObject *yerr_obj;
    
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "b", &ch,&yerr_obj))
        return NULL;
                                       
    value = Read_Single_Channel(ch);
    return Py_BuildValue("l", value);
}


static PyObject *_Init_Single_Channel(PyObject *self, PyObject *args)
{
    unsigned char ch;
    PyObject *yerr_obj;
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "b", &ch,&yerr_obj))
        return NULL;
    
    Init_Single_Channel(ch);
    Py_RETURN_NONE;
}

static PyObject *_ADC_Stop(PyObject *self, PyObject *args)
{
    /* execute the code */ 
    int value = ADC_Stop();
    
    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("i",value);
    return ret;
}

static PyObject *_thread1(PyObject *self, PyObject *args)
{
    double duration, latitude, longitude, elevation;
    char *stationcode, *instrumentstring, *stationchannel,*path;
    int end;
    unsigned char mode;
    
    
    PyObject *yerr_obj;
    
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args,"dbssdddss",&duration,&mode,&stationcode,&stationchannel,&latitude,&longitude,&elevation,&instrumentstring,&path,&yerr_obj))
        return NULL;
    
    /* execute the code */ 
    end = thread1(duration,mode,stationcode,stationchannel,latitude,longitude,elevation,instrumentstring,path);
    
    return Py_BuildValue("i", end);
}

static PyObject *_kill_flag(PyObject *self, PyObject *args)
{
    /* execute the code */ 
    killProgram();
    
    Py_RETURN_NONE;
}
