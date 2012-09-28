

#pragma once

#include <R.h>
#include <Rdefines.h>

#define printf Rprintf

SEXP slice_dc(SEXP gen_ind, SEXP l_ind,SEXP scores,SEXP start,SEXP end);
