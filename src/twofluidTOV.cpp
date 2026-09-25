// cppimport
#include "../include/TOV_family.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
PYBIND11_MODULE(twofluidTOV, m) {
    py::class_<EOS_Tabular>(m, "EOS_Tabular")
        .def(py::init<std::string>())
        .def_readonly("eos_name", &EOS_Tabular::eos_name)
        .def("energy_from_pressure", &EOS_Tabular::energy_from_pressure)
        .def("pc_from_ec", &EOS_Tabular::pc_from_ec)
        .def("dedp", &EOS_Tabular::dedp);
    // Registering the abstract base lets any analytic EOS be passed where EOS_Analytic& is expected
    py::class_<EOS_Analytic>(m, "EOS_Analytic")
        .def_readonly("eos_name", &EOS_Analytic::eos_name)
        .def_readonly("p_surface", &EOS_Analytic::p_surface);
    py::class_<EOS_Poly, EOS_Analytic>(m, "EOS_Poly")
        .def(py::init<std::string, const double, const double>())
        .def("energy_from_pressure", &EOS_Poly::energy_from_pressure)
        .def("pc_from_ec", &EOS_Poly::pc_from_ec)
        .def("dedp", &EOS_Poly::dedp)
        .def_readwrite("K", &EOS_Poly::K);
    py::class_<EOS_SIDM, EOS_Analytic>(m, "EOS_SIDM")
        .def(py::init<const double, const double>())
        .def("energy_from_pressure", &EOS_SIDM::energy_from_pressure)
        .def("pc_from_ec", &EOS_SIDM::pc_from_ec)
        .def("dedp", &EOS_SIDM::dedp)
        .def_readwrite("m_chi", &EOS_SIDM::m_chi)
        .def_readwrite("lambda_chi", &EOS_SIDM::lambda_chi)
        .def_readwrite("lambda", &EOS_SIDM::lambda_chi); // old name, reachable only via getattr
    // TOV_Family and TwoFluid_TOV store references to their EOSs, so keep_alive
    // stops Python from freeing an EOS while a solver still points at it
    py::class_<TOV_Family>(m, "TOV_Family")
        .def(py::init<EOS_Tabular&, EOS_Poly&>(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def(py::init<EOS_Tabular&, EOS_Poly&, const double, const double>(),
             py::arg("eos1"), py::arg("eos2"), py::arg("F_chi"), py::arg("e2_max") = E2_MAX_DEFAULT,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def(py::init<EOS_Tabular&, EOS_SIDM&>(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def(py::init<EOS_Tabular&, EOS_SIDM&, const double, const double>(),
             py::arg("eos1"), py::arg("eos2"), py::arg("F_chi"), py::arg("e2_max") = E2_MAX_DEFAULT,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def(py::init<EOS_Tabular&>(), py::keep_alive<1, 2>())
        .def("calc_lambda", &TOV_Family::calc_lambda)
        .def("calc_lambda_and_mass_directly", &TOV_Family::calc_lambda_and_mass_directly)
        .def("calc_lambda_and_mass_directly_v2", &TOV_Family::calc_lambda_and_mass_directly_v2)
        .def("lambda_from_mass", &TOV_Family::lambda_from_mass)
        .def("set_analytic_eos", &TOV_Family::set_analytic_eos, py::keep_alive<1, 2>())
        .def_readonly("Rs", &TOV_Family::Rs)
        .def_readonly("Ms", &TOV_Family::Ms)
        .def_readonly("k2s", &TOV_Family::k2s)
        .def_readonly("lambdas", &TOV_Family::lambdas)
        .def_readonly("e1s", &TOV_Family::e1s)
        .def_readonly("e2s", &TOV_Family::e2s)
        .def_readonly("RBs", &TOV_Family::RBs)
        .def_readonly("RDs", &TOV_Family::RDs)
        .def_readwrite("e2_max", &TOV_Family::e2_max);
    py::class_<TwoFluid_TOV>(m, "TwoFluid_TOV")
        .def(py::init<EOS_Tabular&, EOS_Poly&>(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def(py::init<EOS_Tabular&, EOS_SIDM&>(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def("integrate_two_fluid_tov", &TwoFluid_TOV::integrate_two_fluid_tov)
        .def("reset_state", &TwoFluid_TOV::reset_state);
    py::class_<TOV_result>(m, "TOV_result")
        .def_readonly("M", &TOV_result::M)
        .def_readonly("R", &TOV_result::R)
        .def_readonly("F_chi", &TOV_result::F_chi)
        .def_readonly("lambda_param", &TOV_result::lambda)
        .def_readonly("R_B", &TOV_result::R_B)
        .def_readonly("R_D", &TOV_result::R_D);
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

