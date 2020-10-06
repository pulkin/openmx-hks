#!/usr/bin/env python
import json
from unittest import TestCase, main
from subprocess import check_output
from tempfile import NamedTemporaryFile

import numpy
from numpy import testing
from scipy.io import loadmat
from scipy.sparse import csr_matrix
import h5py

def load_mat(f):
    return loadmat(f.name, squeeze_me=True)

def load_h5(f):
    return h5py.File(f.name, 'r')

class TestHamiltomnianAgainstReference(TestCase):
    @classmethod
    def setUpClass(cls):
        with open("reference.json", 'r') as f:
            cls.reference = json.load(f)
        for field in ("vectors", "basis_atom", "basis_orbital", "basis_spin", "H", "S"):
            cls.reference[field] = numpy.array(cls.reference[field])
        for field in ("H", "S"):
            cls.reference[field] = cls.reference[field][..., 0] + 1.j * cls.reference[field][..., 1]

    def __test_it__(self, data, convert_complex=False, convert_map=False, sparse=False):
        if convert_map:
            data = {k: numpy.array(v) for k, v in data.items()}

        if convert_complex:
            for field in ("H", "S"):
                data[field] = numpy.array(data[field])
                data[field] = data[field][..., 0] + 1.j * data[field][..., 1]

        if sparse:
            n = len(data["basis_atom"])
            for field in ("H", "S"):
                sparse = csr_matrix((data[field], data[field+"_indices"], data[field+"_indptr"]), (n * len(data["vectors"]), n))
                data[field] = sparse.toarray().reshape(-1, n, n)

        self.assertEqual(self.reference["fermi"], data["fermi"])
        testing.assert_array_equal(self.reference["vectors"], data["vectors"])

        testing.assert_array_equal(self.reference["basis_atom"], data["basis_atom"])
        testing.assert_array_equal(self.reference["basis_orbital"], data["basis_orbital"])
        testing.assert_array_equal(self.reference["basis_spin"], data["basis_spin"])

        testing.assert_array_equal(self.reference["H"], data["H"])
        testing.assert_array_equal(self.reference["S"], data["S"])

    def __test_options__(self, suffix, driver, args=None, **kwargs):
        if args is None:
            args = []
        with NamedTemporaryFile('w+', suffix=suffix) as f:
            print(check_output(['../build/openmx-hks', 'extract-hamiltonian', 'data.hks', f.name] + args).decode('utf-8'))
            f.seek(0)
            data = driver(f)
            self.__test_it__(data, **kwargs)

    def test_json(self):
        self.__test_options__(".json", json.load, convert_complex=True)

    def test_sparse_json(self):
        self.__test_options__(".json", json.load, ["--sparse"], convert_complex=True, sparse=True)

    def test_mat(self):
        self.__test_options__(".mat", load_mat)

    def test_sparse_mat(self):
        self.__test_options__(".mat", load_mat, ["--sparse"], sparse=True)

    def test_h5(self):
        self.__test_options__(".h5", load_h5, convert_complex=True, convert_map=True)

    def test_sparse_h5(self):
        self.__test_options__(".h5", load_h5, ["--sparse"], convert_complex=True, convert_map=True, sparse=True)

if __name__ == "__main__":
    main()
