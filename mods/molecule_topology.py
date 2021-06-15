#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Molecule Topology class
"""

import sys
import os
import parmed
import numpy as np
from pathlib import Path
import copy
import pickle



# =============== classes =============== #
# --------------- MoleculeTopology --------------- #
class MoleculeTopology:
	"""
	分子トポロジークラス
	@param topology_file: トポロジーファイル
	@mask: 有効化マスク
	"""
	def __init__(self, topology_file = None, mask = "*", flag_destructive = False):
		# メンバ変数
		self._obj_topology = None
		self._mask = None
		self._file_type = None

		self._list_target_atom = None
		self._residue_info = []

		# 初期化処理
		self.load_topology(topology_file)
		self.set_mask(mask, flag_destructive)



	def load_topology(self, topology_file):
		"""
		ファイルを読み込むメソッド
		@param topology_file: トポロジーファイル
		@return self
		"""
		if topology_file is not None:
			if ".top" in topology_file:
				# Gromacs topology
				self._obj_topology = parmed.gromacs.GromacsTopologyFile(topology_file)
				self._file_type = "top"
			elif ".gro" in topology_file:
				# Gro file
				self._obj_topology = parmed.gromacs.GromacsGroFile.parse(topology_file)
				self._file_type = "gro"
			elif ".prmtop" in topology_file:
				# Amber topology
				self._obj_topology = parmed.amber.LoadParm(topology_file)
				self._file_type = "prmtop"
			elif ".pdb" in topology_file or ".mol2" in topology_file or ".sdf" in topology_file:
				# Other molecule files
				self._obj_topology = parmed.formats.load_file(topology_file)
				self._file_type = "pdb"
			else:
				sys.stderr.write("ERROR: Unsupported format for topology. ({0})\n".format(topology_file))
				sys.exit(1)
			self._list_target_atom = [idx for idx in range(len(self._obj_topology.atoms))]
		return self


	def set_mask(self, mask = None, flag_destructive = False):
		"""
		マスクを設定するメソッド
		@param mask: AmberMask for enable
		@param flag_destructive: マスク内容をオブジェクトにも適用する (破壊的メソッド) (Default: False)
		@return: self
		"""
		if self._obj_topology is not None:
			self._mask = mask
			self._clear_residue_info()
			if mask is not None:
				self._list_target_atom = list(parmed.amber.AmberMask(self._obj_topology, mask).Selected())
				self._make_residue_info()

				if flag_destructive:
					self._obj_topology.strip("!(" + mask + ")")

		return self


	def _make_residue_info(self):
		"""
		キャッシュデータを作成するメソッド
		"""
		if len(self._residue_info) == 0:
			# 残基情報を登録
			previous_info = ""
			for idx in self._list_target_atom:
				obj_atom = self._obj_topology.atoms[idx]
				info = obj_atom.residue.name + "-" + str(obj_atom.residue.number)
				if previous_info != info:
					# 前回の残基と異なる残基に所属する場合
					self._residue_info.append([obj_atom.residue.name, obj_atom.residue.number, [], [], [], 0.0, 0.0])
					previous_info = info
				self._residue_info[-1][2].append(obj_atom.element)
				self._residue_info[-1][3].append(obj_atom.name)
				self._residue_info[-1][4].append(idx)
				self._residue_info[-1][5] += obj_atom.charge
				self._residue_info[-1][6] += obj_atom.mass

		return self


	def _clear_residue_info(self):
		""" キャッシュデータを削除するメソッド """
		self._residue_info = []
		return self


	def get_file_type(self):
		"""
		読み込んだファイルタイプを返すメソッド
		@return: file_type
		"""
		return self._file_type


	@property
	def file_type(self):
		"""
		読み込んだファイルタイプを返すメソッド
		@return: file_type
		"""
		return self.get_file_type()


	def get_info(self, target, info_type):
		"""
		登録情報を返すメソッド
		@param target: 返す情報の対象 ("atom" or "residue")
		@param info_type: 返す情報 ("atom_index", "atom_name", "element", "charge", "mass", "residue_index", or "residue_name")
		@return 指定された情報はリストで返される
		"""
		if target == "atom":
			# 原子情報を返す
			if info_type == "atom_index":
				# 原子インデックスを返す
				return self._list_target_atom

			elif info_type == "atom_name":
				# 原子名を返す
				return [self._obj_topology.atoms[idx].name for idx in self._list_target_atom]

			elif info_type == "element":
				# 原子番号を返す
				return [self._obj_topology.atoms[idx].element for idx in self._list_target_atom]

			elif info_type == "charge":
				# 電荷を返す
				if self._file_type in ["pdb", "gro"]:
					sys.stderr.write("WARNING: PDB and GRO files do not have any charge information.\n")
				return [self._obj_topology.atoms[idx].charge for idx in self._list_target_atom]

			elif info_type == "mass":
				# 質量を返す
				return [self._obj_topology.atoms[idx].mass for idx in self._list_target_atom]

			elif info_type == "residue_index":
				# 残基インデックスを返す
				return [self._obj_topology.atoms[idx].residue.number for idx in self._list_target_atom]

			elif info_type == "residue_name":
				# 残基名を返す
				return [self._obj_topology.atoms[idx].residue.name for idx in self._list_target_atom]

			else:
				sys.stderr.write("ERROR: Undefined info_type in Molecule.get_info()\n")
				sys.exit(1)

		elif target == "residue":
			# 残基情報を返す
			if info_type == "residue_index":
				# 残基インデックスを返す
				return [x[1] for x in self._residue_info]

			elif info_type == "residue_name":
				# 残基名を返す
				return [x[0] for x in self._residue_info]

			elif info_type == "element":
				# 属する原子インデックスリストを返す
				return [x[2] for x in self._residue_info]

			elif info_type == "charge":
				# 残基の電荷を返す
				if self._file_type in ["pdb", "gro"]:
					sys.stderr.write("WARNING: PDB and GRO files do not have any charge information.\n")
				return [x[5] for x in self._residue_info]

			elif info_type == "mass":
				# 質量を返す
				return [x[6] for x in self._residue_info]

			elif info_type == "atom_index":
				# 属する原子インデックスリストを返す
				return [x[4] for x in self._residue_info]

			elif info_type == "atom_name":
				# 属する原子名リストを返す
				return [x[3] for x in self._residue_info]

			else:
				sys.stderr.write("ERROR: Undefined info_type in Molecule.get_info()\n")
				sys.exit(1)

		else:
			sys.stderr.write("ERROR: Undefined target in Molecule.get_info().\n")
			sys.exit(1)


	def convert_index2name(self, target_type_src, info_type_dst, info_value):
		"""
		原子および残基のインデックスをそれぞれの名前や元素で返すメソッド
		@param target_type_src: 変換したい情報タイプ ("atom_index" or "residue_name")
		@param info_type_dst: 変換後の情報タイプ ("atom_name", "element", or "residue_name")
		@param info_value[]: 変換したいインデックス
		@return 変換後の情報[]
		"""
		if not isinstance(info_value, list):
			# データをすべてリストに変換
			info_value = [info_value]

		results = []
		if target_type_src == "atom_index":
			# 原子インデックスの場合
			if info_type_dst == "atom_name":
				results = [residue_info[3][residue_info[4].index(info_idx)] for residue_info in self._residue_info for info_idx in info_value if info_idx in residue_info[4]]
			elif info_type_dst == "element":
				results = [residue_info[2][residue_info[4].index(info_idx)] for residue_info in self._residue_info for info_idx in info_value if info_idx in residue_info[4]]
			elif info_type_dst == "residue_name":
				info_type_dst = 0
				results = [x[0] for info_idx in info_value for i, x in enumerate(self._residue_info) if info_idx in x[4]]
			else:
				sys.stderr.write("ERROR: Undefined info_type_dst in Molecule.convert_index2name().\n")
				sys.exit(1)

		elif target_type_src == "residue_index":
			# 残基インデックスの場合
			if info_type_dst == "residue_name":
				info_type_dst = 2
			else:
				sys.stderr.write("ERROR: Undefined info_type_dst in Molecule.convert_index2name().\n")
				sys.exit(1)

			for info_idx in info_value:
				# 検索
				results = [residue_info[0] for residue_info in self._residue_info for info_idx in info_value if residue_info[1] == info_idx]

		return results


	def save_file(self, output_file):
		"""
		ファイルとして保存するメソッド
		@param output_file: 出力ファイルパス
		@return self
		"""
		output_path = Path(output_file)
		if output_path.exists():
			output_path.unlink()

		if self._mask is not None:
			# マスクがある場合
			obj_topology = copy.copy(self._obj_topology)
			strip_mask = parmed.amber.AmberMask(obj_topology, "!(" + self._mask + ")")
			obj_topology.strip(strip_mask)
			obj_topology.save(output_file)

		else:
			# マスクがない場合
			self._obj_topology.save(output_file)

		return self


	@property
	def molecule(self):
		"""
		parmed オブジェクトを返すメソッド
		@retrun obj_topology
		"""
		return self._obj_topology


	@property
	def atoms(self):
		"""
		parmed の原子オブジェクトリストを返すメソッド
		@return obj_atoms(list)
		"""
		return self._obj_topology.atoms


	@property
	def residues(self):
		"""
		parmed の残基オブジェクトを返すメソッド
		@return obj_residues(list)
		"""
		return self._obj_topology.residues


	@property
	def mask(self):
		"""
		mask を返すメソッド
		@return mask
		"""
		return self._mask


	@mask.deleter
	def mask(self):
		"""
		mask をリセットするメソッド
		@return self
		"""
		self._mask = None
		return self



# =============== main =============== #
# if __name__ == '__main__':
#	pass
