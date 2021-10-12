#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import copy


# =============== variable =============== #
re_atom = re.compile(r"^\s*\d{,5}.{5}.{5}\s*\d{,5}(?:\s*-?\d+\.\d+){3}")
re_vel = re.compile(r"^(?:\s*-?\d+\.\d{4}){3}")
re_cell = re.compile(r"^(\s*\d+\.\d+){3,}")


# =============== FileGRO class =============== #
class FileGRO:
	""" GRO ファイルクラス """
	def __init__(self, input_file = None):
		self._input_file = ""
		self._obj_frame = []

		# initiation
		self.load_file(input_file)


	def load_file(self, input_file):
		"""
		ファイルを読み込むメソッド
		@param input_file: 入力ファイル
		@return: self
		"""
		if input_file is not None:
			self._input_file = input_file
			flag_title = True
			flag_atom = False
			with open(input_file, "r") as obj_input:
				for line_val in obj_input:
					if flag_title:
						# タイトル
						self._obj_frame.append(FrameGRO())
						self._obj_frame[-1].set_title(line_val.rstrip("\r\n"))
						flag_title = False
						flag_atom = True

					elif flag_atom:
						# 原子数
						self._obj_frame[-1].set_atom(int(line_val.strip()))
						flag_atom = False

					elif re_atom.search(line_val):
						# 原子レコード
						data = [int(line_val[0:5].strip()), line_val[5:10], line_val[10:15], int(line_val[15:20].strip()), float(line_val[21:29].strip()), float(line_val[29:37].strip()), float(line_val[37:45].strip())]

						vel_data = line_val[44:].rstrip("\r\n")
						if re_vel.search(vel_data):
							# 速度データの追加
							data.extend([float(vel_data[0:8]), float(vel_data[8:16]), float(vel_data[16:24])])
						self._obj_frame[-1].add_data(data)

					elif re_cell.search(line_val):
						# セルサイズ
						self._obj_frame[-1].set_cell([float(x) for x in line_val.strip().split()])
						flag_title = True

		return self


	def get_input_file(self):
		"""
		入力ファイルパスを返すメソッド
		@return 入力ファイルパス
		"""
		return self._input_file


	def set_title(self, frame_idx, title):
		"""
		タイトルを設定するメソッド
		@param frame_idx: フレームインデックス
		@param title: タイトル文字列
		@return self
		"""
		self._obj_frame[frame_idx].set_title(title)
		return self


	def get_title(self, frame_idx = None):
		"""
		タイトルを返すメソッド
		@param frame_idx: フレームインデックス
		@return タイトルリスト
		"""
		if frame_idx is not None:
			return [self._obj_frame[frame_idx].get_title()]
		else:
			return [frame.get_title() for frame in self._obj_frame]


	def set_atom(self, frame_idx, atom):
		"""
		原子数を設定するメソッド
		@param frame_idx: フレームインデックス
		@param atom: 原子数
		@return self
		"""
		self._obj_frame[frame_idx].set_atom(atom)
		return self


	def get_atom(self, frame_idx = None):
		"""
		原子数を返すメソッド
		@param frame_idx: フレームインデックス
		@return 原子数リスト
		"""
		if frame_idx is not None:
			return [self._obj_frame[frame_idx].get_atom()]
		else:
			return [frame.get_atom() for frame in self._obj_frame]


	def add_frame(self, coord):
		"""
		フレームデータを追加するメソッド
		@param coord: 座標
		@return self
		"""
		self._obj_frame.append(copy.deepcopy(self._obj_frame[0]))
		frame_idx = len(self._obj_frame) - 1
		self._obj_frame[-1].set_data("coord", coord)
		self._obj_frame[-1].set_title("{0} copied by file_GRO.add_frame() -> {1}".format(self._obj_frame[-1].get_title().rstrip("\r\n"), frame_idx))
		return self


	def get_frame_index(self):
		"""
		フレームインデックスを返すメソッド
		@return frame_index
		"""
		return list(range(len(self._obj_frame)))


	def set_data(self, frame_idx, data_type, data):
		"""
		データを上書きするメソッド
		@param frame_idx: フレームインデックス
		@param data_type: 情報の種類 ("residue_index", "residue_name", "atom_name", "atom_index", "coord", "coord_x", "coord_y", "coord_z", "vel", "vel_x", "vel_y", "vel_z", "*")
		@param data: 原子情報
		@return self
		"""
		self._obj_frame[frame_idx].set_data(data_type, data)
		return self


	def get_data(self, frame_idx, data_type, flag_unique = True):
		"""
		原子情報を返すメソッド
		@param frame_idx: フレームインデックス
		@param data_type: 情報の種類 ("residue_index", "residue_name", "atom_name", "atom_index", "coord", "coord_x", "coord_y", "coord_z", "vel", "vel_x", "vel_y", "vel_z", "*")
		@param data: 原子情報
		@param flag_unique: 残基情報をまとめるかどうか (True: まとめる / False: まとめない)
		@return self
		"""
		if frame_idx is not None:
			return self._obj_frame[frame_idx].get_data(data_type, flag_unique = True)
		else:
			return [frame.get_data(data_type, flag_unique) for frame in self._obj_frame]


	def set_cell(self, frame_idx, cell):
		"""
		セル情報を設定するメソッド
		@param frame_idx: フレームインデックス
		@param cell: セル情報
		@return self
		"""
		self._obj_frame[frame_idx].set_cell(cell)
		return self


	def get_cell(self, frame_idx = None):
		"""
		セル情報を返すメソッド
		@param frame_idx: フレームインデックス
		@return セル情報リスト
		"""
		if frame_idx is not None:
			return [self._obj_frame[frame_idx].get_cell()]
		else:
			return [frame.get_cell() for frame in self._obj_frame]


	def output_gro(self, output_file):
		"""
		GRO ファイルを出力するメソッド
		@param output_file: 出力パス (指定されていない場合、テキストを返す)
		@return: self or .gro formatted text リスト
		"""
		if output_file:
			with open(output_file, "w") as obj_output:
				for obj_frame in self._obj_frame:
					obj_output.write(obj_frame.output_gro())

		else:
			return [frame.output_gro() for frame in self._obj_frame]



# =============== DataGRO class =============== #
class FrameGRO:
	""" GRO フレームクラス """
	def __init__(self):
		# member variables
		self._title = ""
		self._atom = 0
		self._atom_info = []
		self._cell = []


	def set_title(self, title):
		"""
		タイトルを設定するメソッド
		@param title: 設定するタイトル
		@return self
		"""
		self._title = title
		return self


	def set_atom(self, atom):
		"""
		原子数を設定するメソッド
		@param atom:  原子数
		@return self
		"""
		self._atom = atom
		return self


	def add_data(self, data):
		"""
		データを追加するメソッド
		@param data: 追加情報
		@return self
		"""
		self._atom_info.append(data)
		return self


	def set_data(self, data_type, data):
		"""
		データを上書きするメソッド
		@param data_type: データのタイプ
		@param data: データ
		@return: self
		"""
		if data_type == "*":
			# 全データ
			self._atom_info = data

		elif data_type == "residue_index":
			# 残基インデックス
			self._atom_info = [y + x[1:] for x, y in zip(self._atom_info, data)]

		elif data_type == "resdue_name":
			#  残基名
			self._atom_info = [x[0] + y + x[2:] for x, y in zip(self._atom_info, data)]

		elif data_type == "atom_name":
			#  原子名
			self._atom_info = [x[0:2] + y + x[3:] for x, y in zip(self._atom_info, data)]

		elif data_type == "atom_index":
			#  原子インデックス
			self._atom_info = [x[0:3] + y + x[4:] for x, y in zip(self._atom_info, data)]

		elif data_type == "coord":
			#  座標
			self._atom_info = [x[0:4] + y + x[7:] for x, y in zip(self._atom_info, data)]

		elif data_type == "coord_x":
			#  座標 x
			self._atom_info = [x[0:3] + y + x[5:] for x, y in zip(self._atom_info, data)]

		elif data_type == "coord_y":
			#  座標 y
			self._atom_info = [x[0:4] + y + x[6:] for x, y in zip(self._atom_info, data)]

		elif data_type == "coord_z":
			#  座標 z
			self._atom_info = [x[0:5] + y + x[7:] for x, y in zip(self._atom_info, data)]

		elif data_type == "vel":
			#  速度
			self._atom_info = [x[0:6] + y for x, y in zip(self._atom_info, data)]

		elif data_type == "vel_x":
			#  速度 x
			if 8 <= len(self._atom_info[0]):
				self._atom_info = [x[0:6] + y + x[8:] for x, y in zip(self._atom_info, data)]
			else:
				sys.stderr.write("ERROR: adding velocity needs to three dimensional values.\n")
				sys.exit(1)

		elif data_type == "vel_y":
			#  速度 y
			if 8 <= len(self._atom_info[0]):
				self._atom_info = [x[0:7] + y + x[9:] for x, y in zip(self._atom_info, data)]
			else:
				sys.stderr.write("ERROR: adding velocity needs to three dimensional values.\n")
				sys.exit(1)

		elif data_type == "vel_z":
			#  速度 z
			if 8 <= len(self._atom_info[0]):
				self._atom_info = [x[0:8] + y for x, y in zip(self._atom_info, data)]
			else:
				sys.stderr.write("ERROR: adding velocity needs to three dimensional values.\n")
				sys.exit(1)


	def set_cell(self, cell):
		"""
		セル情報を設定するメソッド
		@param cell: セル情報
		@return self
		"""
		self._cell = cell
		return self


	def get_title(self):
		"""
		タイトルを返すメソッド
		@return: title
		"""
		return self._title


	def get_atom(self):
		"""
		原子数を返すメソッド
		@return 原子数
		"""
		return self._atom


	def get_data(self, data_type = "*", flag_unique = True):
		"""
		原子情報を返すメソッド
		@param data_type: 情報の種類 ("residue_index", "residue_name", "atom_name", "atom_index", "coord", "coord_x", "coord_y", "coord_z", "vel", "vel_x", "vel_y", "vel_z", "*")
		@param flag_unique: 残基情報をまとめるかどうか (True: まとめる / False: まとめない)
		@return 原子情報
		"""
		if data_type == "*":
			# 全データ
			return self._atom_info

		elif data_type == "residue_index" or data_type == "reisude_name":
			# 残基番号
			if flag_unique:
				# 残基情報をまとめる
				residue_info_saved = ""
				new_info = []
				for info in self._atom_info:
					residue_info = info[1] + "." + info[0]
					if residue_info_saved != residue_info:
						if data_type == "residue_index":
							# 残基インデックス
							new_info.append(info[0])
						elif data_type == "residue_name":
							# 残基名
							new_info.append(info[1])
						residue_info_saved = residue_info
			else:
				# 残基情報をまとめない
				if data_type == "residue_index":
					return [x[0] for x in self._atom_info]
				elif data_type == "residue_name":
					return [x[1] for x in self._atom_info]

		elif data_type == "atom_name":
			# 原子名
			return [x[2] for x in self._atom_info]

		elif data_type == "atom_index":
			# 原子インデックス
			return [x[3] for x in self._atom_info]

		elif data_type == "coord":
			# 座標
			return [[x[4], x[5], x[6]] for x in self._atom_info]

		elif data_type == "coord_x":
			# 座標 x
			return [x[4]for x in self._atom_info]

		elif data_type == "coord_y":
			# 座標 y
			return [x[5]for x in self._atom_info]

		elif data_type == "coord_z":
			# 座標 z
			return [x[6]for x in self._atom_info]

		elif data_type == "vel":
			# 速度
			if 8 <= len(self._atom_info[0]):
				return [[x[7], x[8], x[9]] for x in self._atom_info]

		elif data_type == "vel_x":
			# 速度 x
			if 8 <= len(self._atom_info[0]):
				return [x[7]for x in self._atom_info]

		elif data_type == "vel_y":
			# 速度 y
			if 8 <= len(self._atom_info[0]):
				return [x[8]for x in self._atom_info]

		elif data_type == "vel_z":
			# 速度 z
			if 8 <= len(self._atom_info[0]):
				return [x[9]for x in self._atom_info]

		else:
			# 未定義
			sys.stderr.write("ERROR: Undefined data_type in fileGRO class.\n")
			sys.exit(1)


	def get_cell(self):
		"""
		セルサイズを返すメソッド
		@return: cell
		"""
		return self._cell


	def output_gro(self, output_file = None):
		"""
		GRO ファイルを出力するメソッド
		@param output_file: 出力パス (指定されていない場合、テキストを返す)
		@return: self or .gro formatted text
		"""
		if output_file:
			with open(output_file, "w") as obj_output:
				obj_output.write("{0}\n".format(self._title))
				obj_output.write("{0}\n".format(self._atom))
				if 8 <= len(self._atom_info[0]):
					# 速度がある場合
					for info in self._atom_info:
						obj_output.write("{0[0]:>5}{0[1]:<5}{0[2]:<5}{0[3]:>5}{0[4]:>8.3f}{0[5]:>8.3f}{0[6]:>8.3f}{0[7]:>8.4f}{0[8]:>8.4f}{0[9]:>8.4f}\n".format(info))

				else:
					for info in self._atom_info:
						obj_output.write("{0[0]:>5}{0[1]:<5}{0[2]:<5}{0[3]:>5}{0[4]:>8.3f}{0[5]:>8.3f}{0[6]:>8.3f}\n".format(info))

				obj_output.write("{0}\n".format("".join(["{0:>10.5f}".format(x) for x in self._cell])))

		else:
			# テキスト形式で返す
			output = ""
			output += "{0}\n".format(self._title)
			output += "{0}\n".format(self._atom)
			if 8 <= len(self._atom_info[0]):
				# 速度がある場合
				for info in self._atom_info:
					output += "{0[0]:>5}{0[1]:<5}{0[2]:<5}{0[3]:>5}{0[4]:>8.3f}{0[5]:>8.3f}{0[6]:>8.3f}{0[7]:>8.4f}{0[8]:>8.4f}{0[9]:>8.4f}\n".format(info)

			else:
				for info in self._atom_info:
					output += "{0[0]:>5}{0[1]:<5}{0[2]:<5}{0[3]:>5}{0[4]:>8.3f}{0[5]:>8.3f}{0[6]:>8.3f}\n".format(info)

			output += "{0}\n".format("".join(["{0:>10.5f}".format(x) for x in self._cell]))
			return output



# =============== main =============== #
# if __name__ == '__main__':
# 	pass
