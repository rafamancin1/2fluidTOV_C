#ifndef GSL_FUNCTION_PP_HPP
#define GSL_FUNCTION_PP_HPP

#include <gsl/gsl_math.h>

// Wraps any callable (e.g. a capturing lambda) as a gsl_function.
// The callable must outlive this wrapper.
template< typename F >  class gsl_function_pp : public gsl_function {
 public:
 gsl_function_pp(const F& func) : _func(func) {
   function = &gsl_function_pp::invoke;
   params=this;
 }
 private:
 const F& _func;
 static double invoke(double x, void *params) {
 return static_cast<gsl_function_pp*>(params)->_func(x);
 }
 };

#endif
