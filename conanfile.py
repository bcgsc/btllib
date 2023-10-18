from conan import ConanFile
from conan.tools.meson import Meson, MesonToolchain


class BtllibRecipe(ConanFile):
    name = "btllib"
    version = "1.6.2"
    url = "https://github.com/bcgsc/btllib"
    description = "Bioinformatics common code library"
    settings = "os", "compiler", "build_type", "arch"
    tool_requires = "meson/1.2.2", "cmake/3.25.2"

    exports_sources = (
        "include/*",
        "recipes/*",
        "scripts/*",
        "src/*",
        "subprojects/*",
        "wrappers/*",
        "meson.build",
    )

    def layout(self):
        self.folders.build = "install"

    def generate(self):
        tc = MesonToolchain(self)
        tc.generate()

    def build(self):
        ms = Meson(self)
        ms.configure()
        ms.build()

    def package(self):
        ms = Meson(self)
        ms.install()

    def package_info(self):
        self.cpp_info.libs = ["btllib"]
