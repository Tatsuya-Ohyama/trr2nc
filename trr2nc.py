#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
trr2nc
Convert Gromacs trajectory (.trr) to AMBER trajectory (.nc)
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import subprocess
import tempfile

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), "modules"))
from py_module_basic import basic
import parmed

# =============== functions =============== #
def check_command(command_name):
	""" コマンドの有無を確認してコマンドの絶対パスを返す関数 """
	process = subprocess.Popen("which %s" % command_name, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	(stdout, stderr) = process.communicate()

	stdout = stdout.decode("utf-8")
	if stdout == "":
		sys.stderr.write("ERROR: Invalid command (%s).\n       Check PATH or installed in the system.\n" % command_name)
		sys.exit(1)
	command_path = stdout.rstrip("\r\n")

	return command_path


def exec_sp(command, operation = False):
	""" subprocess で外部プログラムを実行する関数 """
	if operation:
		process = subprocess.Popen(command, shell = True)
	else:
		process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	streamdata = process.communicate()

	if process.returncode == 1:
		sys.stderr.write("ERROR: subprocess failed\n    '{0}'.\n".format(command))
		sys.exit(1)


def make_prmtop(top, prmtop, strip_mask):
	""" prmtop を作成する関数 """
	sys.stderr.write("{start}Creating prmtop ({file}){end}\n".format(file = prmtop, start = basic.color.LRED + basic.color.BOLD, end = basic.color.END))

	gromacs_top = parmed.gromacs.GromacsTopologyFile(top)
	if strip_mask is not None:
		gromacs_top.strip(strip_mask)
	amber_top = parmed.amber.AmberParm.from_structure(gromacs_top)
	amber_top.write_parm(prmtop)


def make_ndx(top, strip_mask, center_mask):
	""" ndx  ファイルを作成する関数 """
	tempfile_ndx = ""

	with tempfile.NamedTemporaryFile(mode = "w", prefix = ".trr2nc_", dir = args.temp_dir) as obj_temp:
		tempfile_name = obj_temp.name
	tempfile_ndx = tempfile_name + ".ndx"
	sys.stderr.write("{start}Creating ndx ({file}){end}\n".format(file = tempfile_ndx, start = basic.color.LRED + basic.color.BOLD, end = basic.color.END))

	import copy

	with open(tempfile_ndx, "w") as obj_output:
		# 全体の ndx の出力
		obj_output.write("[ System ]\n")
		whole_structure = parmed.gromacs.GromacsTopologyFile(top)

		line_elem = 0
		for atom_idx in [x.idx for x in whole_structure.atoms]:
			line_elem += 1
			if line_elem == 1:
				obj_output.write("{0:>5}".format(atom_idx + 1))
			else:
				obj_output.write(" {0:>5}".format(atom_idx + 1))

			if line_elem % 15 == 0:
				line_elem = 0
				obj_output.write("\n")

		if line_elem % 15 != 0:
			line_elem = 0
			obj_output.write("\n")

		obj_output.write("\n")

		# 中心構造の出力
		obj_output.write("[ Center ]\n")
		if center_mask is None:
			# center_mask が指定されていない場合 (全体と同じ)
			line_elem = 0
			for atom_idx in [x.idx for x in whole_structure.atoms]:
				line_elem += 1
				if line_elem == 1:
					obj_output.write("{0:>5}".format(atom_idx + 1))
				else:
					obj_output.write(" {0:>5}".format(atom_idx + 1))

				if line_elem % 15 == 0:
					line_elem = 0
					obj_output.write("\n")

			if line_elem % 15 != 0:
				line_elem = 0
				obj_output.write("\n")

		else:
			# center_mask がある場合
			line_elem = 0
			for atom_idx in parmed.amber.AmberMask(whole_structure, center_mask).Selected():
				line_elem += 1
				if line_elem == 1:
					obj_output.write("{0:>5}".format(atom_idx + 1))
				else:
					obj_output.write(" {0:>5}".format(atom_idx + 1))

				if line_elem % 15 == 0:
					line_elem = 0
					obj_output.write("\n")

			if line_elem % 15 != 0:
				line_elem = 0
				obj_output.write("\n")

		obj_output.write("\n")

		# strip 構造の出力
		obj_output.write("[ Strip ]\n")
		if strip_mask is None:
			# strip_mask が指定されていない場合 (全体と同じ)
			line_elem = 0
			for atom_idx in [x.idx for x in whole_structure.atoms]:
				line_elem += 1
				if line_elem == 1:
					obj_output.write("{0:>5}".format(atom_idx + 1))
				else:
					obj_output.write(" {0:>5}".format(atom_idx + 1))

				if line_elem % 15 == 0:
					line_elem = 0
					obj_output.write("\n")

			if line_elem % 15 != 0:
				line_elem = 0
				obj_output.write("\n")

		else:
			# strip_mask がある場合
			line_elem = 0
			for atom_idx in parmed.amber.AmberMask(whole_structure, "!" + strip_mask).Selected():
				line_elem += 1
				if line_elem == 1:
					obj_output.write("{0:>5}".format(atom_idx + 1))
				else:
					obj_output.write(" {0:>5}".format(atom_idx + 1))

				if line_elem % 15 == 0:
					line_elem = 0
					obj_output.write("\n")

			if line_elem % 15 != 0:
				line_elem = 0
				obj_output.write("\n")

		obj_output.write("\n")

	return tempfile_ndx


def convert_trajectory(top, tpr, trr, ndx, begin, end, prmtop, strip_mask, output):
	""" トラジェクトリを trr から nc に変換する関数 """
	command_gmx = check_command("gmx")

	tempfile_name = ""
	with tempfile.NamedTemporaryFile(mode = "w", prefix = ".trr2nc_", dir = args.temp_dir) as obj_temp:
		tempfile_name = obj_temp.name

	# 周期境界でジャンプしないトラジェクトリの作成
	temp_traj1 = tempfile_name + "1.trr"
	sys.stderr.write("{start}Creating nojump trajectory ({file}){end}\n".format(file = temp_traj1, start = basic.color.LRED + basic.color.BOLD, end = basic.color.END))
	trajectories = " ".join(trr)
	command = "{command} trjconv -s {tpr} -f {trajectory} -o {output} -pbc nojump -center".format(command = command_gmx, tpr = tpr, trajectory = trajectories, output = temp_traj1)
	if begin is not None:
		command += " -b {0}".format(begin)
	if end is not None:
		command += " -e {0}".format(end)
	command += " -n {0} << 'EOF'\n1\n2\nEOF".format(ndx)
	exec_sp(command, False)

	# strip_mask が指定されていた場合
	if strip_mask is not None:
		sys.stderr.write("{start}Creating stripped tpr file ({file})\n{end}".format(file = tpr, start = basic.color.LRED + basic.color.BOLD, end = basic.color.END))

		# grompp 用の構造ファイルの作成 (strip 済み)
		ref_coord = tempfile_name + "_ref.gro"
		command = "{command} trjconv -s {input_tpr} -f {trajectory} -o {ref_coord} -n {ndx} -b 1 -e 1 << 'EOF'\n2\nEOF".format(command = command_gmx, trajectory = trajectories, input_tpr = tpr, ndx = ndx, ref_coord = ref_coord)
		exec_sp(command, False)

		# strip したトポロジーの作成
		temp_top = tempfile_name + ".top"
		obj_top = parmed.gromacs.GromacsTopologyFile(top)
		obj_top.strip(strip_mask)
		obj_top.write(temp_top)

		# grompp 用の mdp ファイルの作成
		temp_mdp1 = tempfile_name + "1.mdp"
		temp_mdp2 = tempfile_name + "2.mdp"
		with open(temp_mdp1, "w") as obj_output:
			obj_output.write("integrator = steep\n emtol = 1000.0\n emstep = 0.01\n nsteps = 50000\n nstlist = 100\n ns_type = grid\n rlist = 1.0\n coulombtype = PME\n rcoulomb = 1.0\n nstlog = 1 \n pbc = xyz\n vdwtype = cut-off\n constraints = none\n cutoff-scheme = Verlet\n")

		# tpr の作成
		temp_tpr = tempfile_name + ".tpr"
		command = "{command} grompp -f {temp_mdp1} -c {ref_coord} -p {temp_top} -o {temp_tpr} -po {temp_mdp2}".format(command = command_gmx, temp_mdp1 = temp_mdp1, ref_coord = ref_coord, temp_top = temp_top, temp_tpr = temp_tpr, temp_mdp2 = temp_mdp2)
		print(command)
		exec_sp(command, True)

		# 後処理
		tpr = temp_tpr
		os.remove(ref_coord)
		os.remove(temp_top)
		os.remove(temp_mdp1)
		os.remove(temp_mdp2)

	# 分子を中央に配置したトラジェクトリの作成
	temp_traj2 = tempfile_name + "2.trr"
	sys.stderr.write("{start}Creating centered trajectory ({file}){end}\n".format(file = temp_traj2, start = basic.color.LRED + basic.color.BOLD, end = basic.color.END))
	trajectories = " ".join(trr)
	command = "{command} trjconv -s {tpr} -f {trajectory} -o {output} -pbc res -ur compact".format(command = command_gmx, tpr = tpr, trajectory = temp_traj1, output = temp_traj2)
	if begin is not None:
		command += " -b {0}".format(begin)
	if end is not None:
		command += " -e {0}".format(end)
	command += " << 'EOF'\n0\nEOF"
	exec_sp(command, True)
	os.remove(temp_traj1)
	if strip_mask is not None:
		os.remove(tpr)

	# nc ファイルに変換
	temp_in = tempfile_name + ".in"
	sys.stderr.write("{start}Converting AMBER trajectory ({file}){end}\n".format(file = output, start = basic.color.LRED + basic.color.BOLD, end = basic.color.END))
	with open(temp_in, "w") as obj_output:
		obj_output.write("parm {0}\n".format(prmtop))
		obj_output.write("trajin {0}\n".format(temp_traj2))
		obj_output.write("trajout {0}\n".format(output))
		obj_output.write("go\n")

	command_cpptraj = check_command("cpptraj")
	exec_sp("{0} -i {1}".format(command_cpptraj, temp_in), True)

	os.remove(temp_traj2)
	os.remove(temp_in)


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "trr2nc.py - Convert trr to nc with treating PBC", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-s", dest = "tpr", metavar = "INPUT.tpr", required = True, help = "Gromacs run input file")
	parser.add_argument("-x", dest = "trr", metavar = "INPUT.<trr|xtc|gro>", required = True, nargs = "+", help = "Gromacs trajectory file")
	parser.add_argument("-o", dest = "nc", metavar = "OUTPUT.<nc|mdcrd>", required = True, help = "output for Amber trajectory (.nc)")
	parser.add_argument("-p", dest = "prmtop", metavar = "INPUT.prmtop", required = True, help = "Amber topology file")
	parser.add_argument("-t", dest = "top", metavar = "INPUT.top", required = True, help = "Gromacs topology file when prmtop does not exist")
	parser.add_argument("-sc", dest = "temp_dir", metavar = "TEMP_DIR", default = ".", help = "temporary directory (Default: current dir)")

	gmx_option = parser.add_argument_group("gromacs option")
	gmx_option.add_argument("-b", dest = "begin", metavar = "START_TIME", type = int, help = "First frame (ps) to read from trajectory")
	gmx_option.add_argument("-e", dest = "end", metavar = "END_TIME", type = int, help = "Last frame (ps) to read from trajectory")

	cpptraj_option = parser.add_argument_group("cpptraj option")
	cpptraj_option.add_argument("-mc", dest = "center_mask", metavar = "CENTER_MASK", help = "center mask for cpptraj")
	cpptraj_option.add_argument("-ms", dest = "strip_mask", metavar = "STRIP_MASK", help = "strip mask for cpptraj")

	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")

	args = parser.parse_args()

	basic.check_exist(args.tpr, 2)
	for trj_file in args.trr:
		basic.check_exist(trj_file, 2)
	basic.check_exist(args.top, 2)

	if args.flag_overwrite == False:
		basic.check_overwrite(args.prmtop)
	if os.path.exists(args.prmtop):
		os.remove(args.prmtop)
	make_prmtop(args.top, args.prmtop, args.strip_mask)
	ndx = make_ndx(args.top, args.strip_mask, args.center_mask)

	convert_trajectory(args.top, args.tpr, args.trr, ndx, args.begin, args.end, args.prmtop, args.strip_mask, args.nc)

	if args.strip_mask is not None:
		os.remove(ndx)
