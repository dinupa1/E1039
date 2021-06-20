cmake_minimum_required(VERSION 3.19)
project(analysis)

set(CMAKE_CXX_STANDARD 11)

find_package(ROOT)

add_executable(analysis main.cpp)

target_link_libraries(analysis ROOT::Core ROOT::Imt ROOT::RIO ROOT::Net ROOT::Hist
        ROOT::Graf ROOT::Graf3d ROOT::Gpad ROOT::ROOTVecOps ROOT::Tree
        ROOT::TreePlayer ROOT::Rint ROOT::Postscript ROOT::Matrix ROOT::Physics
        ROOT::MathCore ROOT::Thread ROOT::MultiProc ROOT::ROOTDataFrame)