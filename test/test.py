#!/usr/bin/env python
import json
from unittest import TestCase, main
from subprocess import check_output
from tempfile import NamedTemporaryFile

import numpy
from numpy import testing
from scipy.io import loadmat
import h5py

class TestHamiltomnianAgainstReference(TestCase):
    @classmethod
    def setUpClass(cls):
        with open("reference.json", 'r') as f:
            cls.reference = json.load(f)
        for field in ("vectors", "basis_atom", "basis_orbital", "basis_spin", "H", "S"):
            cls.reference[field] = numpy.array(cls.reference[field])
        for field in ("H", "S"):
            cls.reference[field] = cls.reference[field][..., 0] + 1.j * cls.reference[field][..., 1]

    def __test_it__(self, data, convert_complex=False, convert_map=False):
        if convert_map:
            data = {k: numpy.array(v) for k, v in data.items()}

        if convert_complex:
            for field in ("H", "S"):
                data[field] = numpy.array(data[field])
                data[field] = data[field][..., 0] + 1.j * data[field][..., 1]

        self.assertEqual(self.reference["fermi"], data["fermi"])
        testing.assert_array_equal(self.reference["vectors"], data["vectors"])

        testing.assert_array_equal(self.reference["basis_atom"], data["basis_atom"])
        testing.assert_array_equal(self.reference["basis_orbital"], data["basis_orbital"])
        testing.assert_array_equal(self.reference["basis_spin"], data["basis_spin"])

        testing.assert_array_equal(self.reference["H"], data["H"])
        testing.assert_array_equal(self.reference["S"], data["S"])

    def test_json(self):
        f = NamedTemporaryFile('w+', suffix='.json')
        print(check_output(['../build/openmx-hks', 'extract-hamiltonian', 'data.hks', f.name]).decode('utf-8'))
        f.seek(0)
        data = json.load(f)
        self.__test_it__(data, convert_complex=True)

    def test_mat(self):
        f = NamedTemporaryFile('w+', suffix='.mat')
        print(check_output(['../build/openmx-hks', 'extract-hamiltonian', 'data.hks', f.name]).decode('utf-8'))
        f.seek(0)
        data = loadmat(f.name, squeeze_me=True)
        self.__test_it__(data)

    def test_h5(self):
        f = NamedTemporaryFile('w+', suffix='.h5')
        print(check_output(['../build/openmx-hks', 'extract-hamiltonian', 'data.hks', f.name]).decode('utf-8'))
        f.seek(0)
        data = h5py.File(f.name, 'r')
        self.__test_it__(data, convert_complex=True, convert_map=True)

if __name__ == "__main__":
    main()
