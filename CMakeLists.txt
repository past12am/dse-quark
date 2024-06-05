cmake_minimum_required(VERSION 3.27)
project(QuarkDSE_v3)

set(CMAKE_CXX_STANDARD 17)

add_executable(QuarkDSE_v3 main.cpp
        src/grid/MomentumGrid.cpp
        include/grid/MomentumGrid.hpp
        include/Definitions.hpp
        src/dse/DSE.cpp
        include/dse/DSE.hpp
        src/qcd/QuarkPropagator.cpp
        include/qcd/QuarkPropagator.hpp
        src/qcd/QuarkSelfEnergy.cpp
        include/qcd/QuarkSelfEnergy.hpp
        include/numerics/roots/NewtonRootFinder.hpp
        include/numerics/integration/GaussChebyshev.hpp
        include/numerics/integration/GaussLegendre.hpp
        include/numerics/polynomials/LegendrePolynomials.hpp
        src/numerics/roots/NewtonRootFinder.cpp
        src/numerics/integration/GaussChebyshev.cpp
        src/numerics/integration/GaussLegendre.cpp
        src/numerics/polynomials/LegendrePolynomials.cpp
)

find_package(GSL REQUIRED)
target_link_libraries(QuarkDSE_v3 GSL::gsl)