// cppimport
#include "../include/TOV_family.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
PYBIND11_MODULE(twofluidTOV, m) {
    py::class_<EOS_Tabular>(m, "EOS_Tabular")
        .def(py::init<std::string>());
    py::class_<EOS_Poly>(m, "EOS_Poly")
        .def(py::init<std::string, const double, const double>());
    py::class_<TOV_Family>(m, "TOV_Family")
        .def(py::init<EOS_Tabular&, EOS_Poly&>())
        .def(py::init<EOS_Tabular&, EOS_Poly&, const double>())
        .def(py::init<EOS_Tabular&>())
        .def("calc_lambda_parallel", &TOV_Family::calc_lambda_parallel)
        .def("calc_lambda", &TOV_Family::calc_lambda)
        .def("calc_lambda_and_mass_directly", &TOV_Family::calc_lambda_and_mass_directly)
        .def_readonly("Rs", &TOV_Family::Rs)
        .def_readonly("Ms", &TOV_Family::Ms)
        .def_readonly("k2s", &TOV_Family::k2s)
        .def_readonly("e1s", &TOV_Family::e1s)
        .def_readonly("e2s", &TOV_Family::e2s)
        .def_readonly("RBs", &TOV_Family::RBs)
        .def_readonly("RDs", &TOV_Family::RDs);
    py::class_<TwoFluid_TOV>(m, "TwoFluid_TOV")
        .def(py::init<EOS_Tabular&, EOS_Poly&>())
        .def("integrate_two_fluid_tov", &TwoFluid_TOV::integrate_two_fluid_tov);
    py::class_<TOV_result>(m, "TOV_result")
        .def_readonly("M", &TOV_result::M)
        .def_readonly("R", &TOV_result::R)
        .def_readonly("F_chi", &TOV_result::F_chi)
        .def_readonly("lambda_param", &TOV_result::lambda);
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

