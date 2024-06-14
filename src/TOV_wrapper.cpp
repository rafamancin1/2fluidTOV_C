#include "../include/TOV_family.hpp"
//#include <boost/python.hpp>
#include <pybind11/pybind11.h>
/*
using namespace boost::python;

BOOST_PYTHON_MODULE(TOV_ext) {
    class_<TOV_Family>("TOV_Family", init<>())
    .def("calc_lambda_parallel", &TOV_Family::calc_lambda_parallel)
    .def("calc_lambda", &TOV_Family::calc_lambda)
    ;
}
*/
namespace py = pybind11;
PYBIND11_MODULE(twofluidTOV, m) {
    py::class_<EOS_Tabular>(m, "EOS_Tabular")
        .def(py::init<std::string>());
    py::class_<EOS_Poly>(m, "EOS_Poly")
        .def(py::init<std::string, const double, const double>());
    py::class_<TOV_Family>(m, "TOV_Family")
        .def(py::init<EOS_Tabular&, EOS_Poly&>())
        .def("calc_lambda_parallel", &TOV_Family::calc_lambda_parallel)
        .def("calc_lambda", &TOV_Family::calc_lambda);
}


