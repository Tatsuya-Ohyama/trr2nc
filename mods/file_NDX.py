#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pickle

from molecule_topology import MoleculeTopology


# =============== class =============== #
class FileNDX:
	""" NDXFile class """
	def __init__(self, obj_mol):
		# member variables
		self._obj_mol = obj_mol
		self._defs = []

		# initiation
		self.add_def("System", "*")


	def save_pickle(self, output_file):
		"""
		Pickle ファイルに保存するメソッド
		@param output_file: 出力する pickle ファイルのパス
		@return: 自身を返す (チェーンメソッドのため)"""
		with open(output_file, "wb") as obj_output:
			pickle.dump(self, obj_output)
		return self


	def restore_pickle(self, input_file):
		"""
		Pickle ファイルから復元するメソッド
		@param input_file: pickle ファイルのパス
		@return: 自身を返す (チェーンメソッドのため)
		"""
		with open(input_file, "rb") as obj_input:
			self = pickle.load(obj_input)
		return self


	def add_def(self, name, amber_mask):
		"""
		原子群を定義を追加するメソッド
		@param name: 定義名
		@param amber_mask: 定義する原子の AmberMask
		@return self
		"""
		self._obj_mol.set_mask(amber_mask)
		self._defs.append([name, self._obj_mol.get_info("atom", "atom_index")])
		return self


	def del_def(self, index):
		"""
		原子群の定義を削除するメソッド
		@param index: 原子群の定義のインデックス
		@return self
		"""
		del(self._defs[index])
		return self


	def get_def(self, name = None):
		"""
		定義された原子群を表示するメソッド
		@param name: 定義名 (Default: None)
		@return 定義名と原子群のインデックスリスト
		"""
		if name is not None:
			# 定義名がある場合
			list_check = [i for i, x in enumerate(self._def) if x[0] == name]
			if len(list_check) != 0:
				return self._defs[list_check[0]]
			else:
				sys.stderr.write("ERROR: Invalid define name.\n")
				return False
		else:
			# 定義名がない場合
			return self._defs


	def output_ndx(self, output_file):
		"""
		ndx ファイルを出力するメソッド
		@param output_file: 出力する ndx のファイルパス
		@return: self
		"""
		with open(output_file, "w") as obj_output:
			for list_def in self._defs:
				obj_output.write("[ {0} ]\n".format(list_def[0]))
				list_def[1] = ["{0:>4d}".format(x + 1) for x in list_def[1]]
				list_range = [0, 15]
				while list_range[0] < len(list_def[1]):
					obj_output.write("{0}\n".format(" ".join([str(x) for x in list_def[1][list_range[0]:list_range[1]]])))
					list_range = [x + 15 for x in list_range]

		return self


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
