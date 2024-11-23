#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void aHD_C_double_val(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void aHD_C_int(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void aHD_C_int_val(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void aHD_R_double_val(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void aHD_R_int(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void aHD_R_int_val(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"aHD_C_double_val", (DL_FUNC) &aHD_C_double_val, 12},
  {"aHD_C_int",        (DL_FUNC) &aHD_C_int,        11},
  {"aHD_C_int_val",    (DL_FUNC) &aHD_C_int_val,    12},
  {"aHD_R_double_val", (DL_FUNC) &aHD_R_double_val, 12},
  {"aHD_R_int",        (DL_FUNC) &aHD_R_int,        11},
  {"aHD_R_int_val",    (DL_FUNC) &aHD_R_int_val,    12},
  {NULL, NULL, 0}
};

void R_init_sphericalDepth(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
