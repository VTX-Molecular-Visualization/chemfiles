import os
from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout
from conan.tools.scm import Git

class VTXChemfilesRecipe(ConanFile):
    name = "chemfiles"
    version = "1.0"
    package_type = "library"
    
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}
    
    generators = "CMakeDeps", "CMakeToolchain"
    
    #exports_sources = "CMakeLists.txt", "src/*", "include/*"
        
    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def layout(self):
        cmake_layout(self)    
        # Add generated include dir.        
        self.cpp.source.includedirs = ["include", os.path.join(self.folders.build, "include")]
        #self.cpp.package.includedirs = ["include", os.path.join(self.folders.build, "include")]
        
    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["chemfiles"]
        if self.settings.os == "Windows":
            self.cpp_info.system_libs.append('ws2_32')

