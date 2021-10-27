#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import parmed



# =============== constant =============== #
VAL_PER_ROW = 15
VAL_LENGTH = 4
VAL_FORMAT = "{0:>" + str(VAL_LENGTH) + "}"



# =============== class =============== #
class FileNDX:
	""" NDXFile class """
	def __init__(self, obj_mol):
		# member variables
		self._obj_mol = obj_mol
		self._def_names = []
		self._def_list = {}

		# initiation
		self.add_def("System", "*")


	def add_def(self, name, amber_mask):
		"""
		Method to define atom group

		Args:
			name(str): define name
			amber_mask(str): AmberMask

		Returns:
			self
		"""
		obj_mask = parmed.amber.AmberMask(self._obj_mol, amber_mask)
		self._def_names.append(name)
		self._def_list[name] = [i + 1 for i in obj_mask.Selected()]
		return self


	def del_def(self, name):
		"""
		Method to delete defined atom group

		Args:
			name(str or int): name or index for atom group

		Returns:
			self
		"""
		if isinstance(name) == str:
			del(self._def_list[name])
			self._def_names.remove(name)

		elif isinstance(name) == int:
			del(self._def_list[self._def_names[name]])
			del(self._def_names[name])

		else:
			sys.stderr.write("ERROR: undefined data type (str or int).\n")
			sys.exit(1)

		return self


	def get_def(self, name=None):
		"""
		Method to return defined atom groups

		Args:
			name(str or int): name or index for atom group

		Returns:
			dict: defined name and list for atomic index
		"""
		if isinstance(name) == str:
			if name in self._def_list.keys():
				return self._def_list[name]
			else:
				sys.stderr.write("ERROR: Invalid define name.\n")
				return False

		elif isinstance(name) == int:
			return self._def_list[self._def_names[name]]

		else:
			# 定義名がない場合
			return self._def_list


	def output_ndx(self, output_file):
		"""
		ndx ファイルを出力するメソッド
		@param output_file: 出力する ndx のファイルパス
		@return: self
		"""
		with open(output_file, "w") as obj_output:
			for name in self._def_names:
				obj_output.write("[ {0} ]\n".format(name))
				for i in range(0, len(self._def_list[name]), VAL_PER_ROW):
					row = " ".join([VAL_FORMAT.format(v) for v in self._def_list[name][i:i+VAL_PER_ROW]])
					obj_output.write("{0}\n".format(row))

		return self
