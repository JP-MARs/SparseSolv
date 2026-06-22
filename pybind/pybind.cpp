
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/complex.h>

#include <SparseSolv/SparseBuilder.hpp>
#include <SparseSolv/SparseMat.hpp>
#include <SparseSolv/MatSolverBase.hpp>
#include <SparseSolv/MatSolverICCG.hpp>
#include <SparseSolv/MatSolverMRTR.hpp>

namespace py = pybind11;
using namespace JPMRspace::SparseSolv;

PYBIND11_MODULE(SparseSolvPy, m) {
    m.doc() = "SparseSolv bindings";

    py::class_<SparseBuilder>(m, "SparseBuilder")
        .def(py::init<slv_int>(), py::arg("size"))
        .def("getSize", &SparseBuilder::getSize)
        .def("initialize", py::overload_cast<slv_int>(&SparseBuilder::initialize), py::arg("size"))
        .def("reInitialize", &SparseBuilder::reInitialize, py::arg("size"))
        .def("add", &SparseBuilder::add, py::arg("row"), py::arg("col"), py::arg("val"))
        .def("resetMat", &SparseBuilder::resetMat)
    ;

    py::class_<SparseBuilderC>(m, "SparseBuilderC")
        .def(py::init<slv_int>(), py::arg("size"))
        .def("get_size", &SparseBuilderC::getSize)
        .def("initialize", py::overload_cast<slv_int>(&SparseBuilderC::initialize), py::arg("size"))
        .def("reInitialize", &SparseBuilderC::reInitialize, py::arg("size"))
        .def("add", &SparseBuilderC::add, py::arg("row"), py::arg("col"), py::arg("val"))
        .def("resetMat", &SparseBuilderC::resetMat)
    ;

    py::class_<SparseMat>(m, "SparseMat")
        .def(py::init<>())
        .def(py::init<SparseBuilder&, bool>(), py::arg("builder"), py::arg("to_square") = true)
        .def(py::init<slv_int,
                      const std::vector<slv_int>&,
                      const std::vector<slv_int>&,
                      const std::vector<double>&>(),
             py::arg("n_nonzero"),
             py::arg("rows"),
             py::arg("cols"),
             py::arg("vals"))
        .def("getSize", &SparseMat::getSize)
        .def("getMaxCol", &SparseMat::getMaxCol)
        .def("print_mat", &SparseMat::printMat, py::arg("filename") = "Mat.csv")
    ;

    py::class_<SparseMatC>(m, "SparseMatC")
        .def(py::init<>())
        .def(py::init<SparseBuilderC&, bool>(), py::arg("builder"), py::arg("to_square") = true)
        .def(py::init<slv_int,
                      const std::vector<slv_int>&,
                      const std::vector<slv_int>&,
                      const std::vector<JPMRspace::dcomplex>&>(),
             py::arg("n_nonzero"),
             py::arg("rows"),
             py::arg("cols"),
             py::arg("vals"))
        .def("getSize", &SparseMatC::getSize)
        .def("getMaxCol", &SparseMatC::getMaxCol)
        .def("print_mat", &SparseMatC::printMat, py::arg("filename") = "Mat.csv")
    ;

	py::class_<MatSolverBase>(m, "MatSolverBase")
        .def("setInitSolution", &MatSolverBase::setInitSolution)
        .def("setDiagScale", &MatSolverBase::setDiagScale)
        .def("setSaveBest", &MatSolverBase::setSaveBest)
        .def("setSaveLog", &MatSolverBase::setSaveLog)
        .def("setDirvegeType", &MatSolverBase::setDirvegeType)
        .def("setBadDivVal", &MatSolverBase::setBadDivVal)
        .def("setBadDivCount", &MatSolverBase::setBadDivCount)
        .def("setConvNormalizeType", &MatSolverBase::setConvNormalizeType)
        .def("setConvNormalizeConst", &MatSolverBase::setConvNormalizeConst)
        .def("getResidualLog", [](MatSolverBase& self) {
            std::vector<double> log;
            self.getResidualLog(log);
            return log;
        })
    ;

    py::class_<MatSolverICCG, MatSolverBase>(m, "MatSolverICCG")
        .def(py::init<>())
        .def("setAccelFactor", &MatSolverICCG::setAccelFactor)
        .def("solve_real",
             [](MatSolverICCG& self,
                double conv_cri,
                int max_ite,
                const SparseMat& A,
                const Eigen::VectorXd& b,
                Eigen::VectorXd x0) {
                 bool ok = self.solve(conv_cri, max_ite, A, b, x0);
                 return py::make_tuple(ok, x0);
             },
             py::arg("conv_cri"), py::arg("max_ite"), py::arg("A"), py::arg("b"), py::arg("x0"))
        .def("solve_complex",
             [](MatSolverICCG& self,
                double conv_cri,
                int max_ite,
                const SparseMatC& A,
                const Eigen::VectorXcd& b,
                Eigen::VectorXcd x0) {
                 bool ok = self.solve(conv_cri, max_ite, A, b, x0);
                 return py::make_tuple(ok, x0);
             },
             py::arg("conv_cri"), py::arg("max_ite"), py::arg("A"), py::arg("b"), py::arg("x0"))
    ;

    py::class_<MatSolverMRTR, MatSolverBase>(m, "MatSolverMRTR")
        .def(py::init<>())
        .def("setAccelFactor", &MatSolverMRTR::setAccelFactor)
        .def("setPreconditionType", &MatSolverMRTR::setPreconditionType)
        .def("solve_real",
             [](MatSolverMRTR& self,
                double conv_cri,
                int max_ite,
                const SparseMat& A,
                const Eigen::VectorXd& b,
                Eigen::VectorXd x0) {
                 bool ok = self.solve(conv_cri, max_ite, A, b, x0);
                 return py::make_tuple(ok, x0);
             },
             py::arg("conv_cri"), py::arg("max_ite"), py::arg("A"), py::arg("b"), py::arg("x0"))
        .def("solve_complex",
             [](MatSolverMRTR& self,
                double conv_cri,
                int max_ite,
                const SparseMatC& A,
                const Eigen::VectorXcd& b,
                Eigen::VectorXcd x0) {
                 bool ok = self.solve(conv_cri, max_ite, A, b, x0);
                 return py::make_tuple(ok, x0);
             },
             py::arg("conv_cri"), py::arg("max_ite"), py::arg("A"), py::arg("b"), py::arg("x0"))
    ;

}

