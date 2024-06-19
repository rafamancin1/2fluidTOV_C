// cppimport
#include "../include/TOV_family.hpp"
#include <pybind11/pybind11.h>

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


/*
<%
cfg['dependencies'] = ['conversions.hpp', 'eos_poly.hpp', 'eos_tabular.hpp', 'eos.hpp', 'spline.h', 'TOV_family.hpp', 'twofluid_TOV.hpp']
cfg['extra_link_args'] = ['`gsl-config --cflags --libs`']
cfg['libraries'] = ['gsl']
cfg['include_dirs'] = ['libInterpolate', 'Eigen', 'boost']
setup_pybind11(cfg)
%>
*/

