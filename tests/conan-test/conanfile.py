import os

from conan import ConanFile
from conan.tools.meson import Meson
from conan.tools.build import can_run


class BtllibTestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "PkgConfigDeps", "MesonToolchain"

    def requirements(self):
        self.requires(self.tested_reference_str)

    def build(self):
        ms = Meson(self)
        ms.configure()
        ms.build()

    def test(self):
        if can_run(self):
            cmd = os.path.join(self.cpp.build.bindir, "btllib-conan-test")
            self.run(cmd, env="conanrun")
