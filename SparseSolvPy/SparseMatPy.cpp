
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <memory>

#include "SparseMat.hpp"
#include "MatSolversICCG.hpp"

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
	py::class_<MatSolversICCG>(m, "MatSolversICCG")
		.def(py::init<>())
		.def("solveICCGPy", [](int size, const double conv_cri, const int max_ite, const double accera, const SparseMat& A, const py::array_t<double>& b, py::array_t<double>& x, bool is_log){
			MatSolversICCG iccg;
			iccg.setSaveLog(is_log);
			bool bl = iccg.solveICCG( size, conv_cri, max_ite, accera, A, b.data(), const_cast<double*>(x.data()) );
			std::vector<double> log;
			iccg.getResidualLog(log);
			return log;
		})
		//.def("solveICCG", py::overload_cast<int, const double, const int, const double, const SparseMat&, const std::vector<double>&, std::vector<double>&, bool>(&MatSolversICCG::solveICCG))
		//.def("solveICCGwithABMC", py::overload_cast<int, const double, const int, const double, const SparseMat&, const std::vector<double>&, std::vector<double>&, int, int, bool>(&MatSolversICCG::solveICCGwithABMC))
	;
}

