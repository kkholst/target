#include <string>
#include <complex>
#include <memory>     // smart pointers (unique_ptr)
#include <vector>
#include <pybind11/stl.h>
#include <target/glm.hpp>
#include "armapy.hpp"
#include "py_test.hpp"
#include "loop.hpp"
#include <algorithm>
#include <string>
#include <random>
#include <iterator>
#include <torch/torch.h>

void torch_test() {
  torch::Tensor tensor = torch::rand({2, 3});
  std::cout << tensor << std::endl;
}

pyarray expit(pyarray &x) {
  arma::mat res = target::expit(pymat(x));
  return matpy(res);
}

std::string samp() {
  std::string in = "Sangmin", out;
  std::sample(in.begin(), in.end(), std::back_inserter(out),
    5, std::mt19937{std::random_device{}()});
  return out;
}

PYBIND11_MODULE(__targeted_template_c__, m) {
  m.doc() = "Python bindings for the target C++ library";
  m.def("expit", &expit, "Sigmoid function (inverse logit)");
  //m.#define ("expit", [](arma::mat &x) { return target::expit(x); }, "Sigmoid function (inverse logit)");
  m.def("myloop", &myloop, "test loop");
  m.def("scale2", &scale2);
  m.def("add", &add);
  m.def("samp", &samp);
  m.def("torch_test", &torch_test);

  py::class_<Test>(m, "testclass")
    .def(py::init<std::vector<double> &, std::vector<double> &,  // x, y
       std::string>())                   // model-type: 'default', 'advanced'
    .def_readwrite("X", &Test::x, "x variable")
    .def_readwrite("Y", &Test::y, "y variable")
    .def("calc", &Test::calc);
}
