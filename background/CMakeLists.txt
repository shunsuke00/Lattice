# CMakeのバージョンを設定
cmake_minimum_required(VERSION 3.13)

set ( CMAKE_CXX_COMPILER g++)
set ( CMAKE_CXX_FLAGS "-fopenmp" )

# プロジェクト名と使用する言語を設定
project(test_cmake CXX)

# a.outという実行ファイルをadd.cppから作成
add_executable(a.out main.cpp)
