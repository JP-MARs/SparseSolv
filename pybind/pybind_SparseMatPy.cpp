
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <memory>

#include <SparseSolve/SparseMat.hpp>
#include <SparseSolve/MatSolvers.hpp>

using namespace SRLfem;
namespace py = pybind11;

PYBIND11_MODULE(SparseSolvPy, m)
{
	// publish SparseMat to Python
	py::class_<SparseMat>(m, "SparseMat")
		.def(py::init<int>())
		.def(py::init<int, const std::vector<int>&,
						const std::vector<int>&,
						const std::vector<double>&>())
		.def("add", &SparseMat::add)
		.def("resetMat", &SparseMat::resetMat)
		.def("fix", &SparseMat::fix)
		.def("delFlagPosition", &SparseMat::delFlagPosition)
		.def("DiagFlagPosition", &SparseMat::DiagFlagPosition)
		.def("print", &SparseMat::print)
		.def("printMat", &SparseMat::printMat)
	;
	// publish SparseMatC to Python
	py::class_<SparseMatC>(m, "SparseMatC")
		.def(py::init<int>())
		.def(py::init<int, const std::vector<int>&,
						const std::vector<int>&,
						const std::vector<dcomplex>&>())
		.def("add", &SparseMatC::add)
		.def("resetMat", &SparseMatC::resetMat)
		.def("fix", &SparseMatC::fix)
		.def("delFlagPosition", &SparseMatC::delFlagPosition)
		.def("DiagFlagPosition", &SparseMatC::DiagFlagPosition)
		.def("print", &SparseMatC::print)
		.def("printMat", &SparseMatC::printMat)
	;
	// publish SparseMatC to Python
	py::class_<MatSolvers>(m, "MatSolvers")
		.def(py::init<>())
		.def("setDiagScale", &MatSolvers::setDiagScale)
		.def("setSaveBest", &MatSolvers::setSaveBest)
		.def("setSaveLog", &MatSolvers::setSaveLog)
		.def("getResidualLog", &MatSolvers::getResidualLog)
		.def("setDirvegeType", &MatSolvers::setDirvegeType)
		.def("setBadDivVal", &MatSolvers::setBadDivVal)
		.def("setBadDivCount", &MatSolvers::setBadDivCount)
		.def("setConvNormalizeType", &MatSolvers::setConvNormalizeType)
		.def("setConvNormalizeConst", &MatSolvers::setConvNormalizeConst)

		.def("getResidualLog_py", &MatSolvers::getResidualLog_py)
		.def("solveICCG_py", py::overload_cast<const int, const double, const int, const double, const SparseMat&, py::array_t<double>, py::array_t<double>, bool>(&MatSolvers::solveICCG_py) )
		.def("solveICCG_py", py::overload_cast<const int, const double, const int, const double, const SparseMatC&, py::array_t<std::complex<double>>, py::array_t<std::complex<double>>, bool>(&MatSolvers::solveICCG_py) )
		.def("solveICCGwithABMC_py", py::overload_cast<const int, const double, const int, const double, const SparseMat&, py::array_t<double>, py::array_t<double>, int, int, bool>(&MatSolvers::solveICCGwithABMC_py) )
		.def("solveICMRTR_py", py::overload_cast<const int, const double, const int, const double, const SparseMat&, py::array_t<double>, py::array_t<double>, bool>(&MatSolvers::solveICMRTR_py) )
		.def("solveICMRTR_py", py::overload_cast<const int, const double, const int, const double, const SparseMatC&, py::array_t<std::complex<double>>, py::array_t<std::complex<double>>, bool>(&MatSolvers::solveICMRTR_py) )
		.def("solveSGSMRTR_py", py::overload_cast<const int, const double, const int, const SparseMat&, py::array_t<double>, py::array_t<double>, bool>(&MatSolvers::solveSGSMRTR_py) )
		.def("solveEigenBiCGstab_accel_py", py::overload_cast<const int, const double, const int, const SparseMatC&, py::array_t<std::complex<double>>, py::array_t<std::complex<double>>, bool>(&MatSolvers::solveEigenBiCGstab_accel_py) )
		.def("solveEigenBiCGstab_accel_py", py::overload_cast<const int, const double, const int, const SparseMat&, py::array_t<double>, py::array_t<double>, bool>(&MatSolvers::solveEigenBiCGstab_accel_py) )
		.def("solveEigenGMRES_py", py::overload_cast<const int, const double, const int, const SparseMatC&, py::array_t<std::complex<double>>, py::array_t<std::complex<double>>, bool>(&MatSolvers::solveEigenGMRES_py) )
		.def("solveEigenGMRES_py", py::overload_cast<const int, const double, const int, const SparseMat&, py::array_t<double>, py::array_t<double>, bool>(&MatSolvers::solveEigenGMRES_py) )
	;
}



