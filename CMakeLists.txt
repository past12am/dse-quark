cmake_minimum_required(VERSION 3.27)
project(QuarkDSE_v3)

set(CMAKE_CXX_STANDARD 17)

find_package(GSL REQUIRED)

add_executable(QuarkDSE_v3 main.cpp
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
        src/grid/ComplexMomentumGrid.cpp
        include/grid/ComplexMomentumGrid.hpp
        include/helpers/InterpAcceleratorsStruct.hpp
)

target_link_libraries(QuarkDSE_v3 GSL::gsl)