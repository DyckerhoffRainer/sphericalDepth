#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void aHD_double_val(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void aHD_int(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void aHD_int_val(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"aHD_double_val", (DL_FUNC) &aHD_double_val, 13},
  {"aHD_int",        (DL_FUNC) &aHD_int,        12},
  {"aHD_int_val",    (DL_FUNC) &aHD_int_val,    13},
  {NULL, NULL, 0}
};

void R_init_sphericalDepth(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
