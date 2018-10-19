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
from termcolor import colored

from basic_func import check_exist, check_overwrite
from molecule_topology import MoleculeTopology
from file_NDX import NDXFile


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


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "trr2nc.py - Convert trr to nc with treating PBC", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-s", dest = "tpr", metavar = "INPUT.tpr", required = True, help = "Gromacs run input file")
	parser.add_argument("-x", dest = "trr", metavar = "INPUT.<trr|xtc|gro>", required = True, nargs = "+", help = "Gromacs trajectory file")
	parser.add_argument("-o", dest = "nc", metavar = "OUTPUT.<nc|mdcrd>", required = True, help = "output for Amber trajectory (.nc)")
	parser.add_argument("-t", dest = "top", metavar = "INPUT.top", required = True, help = "Gromacs topology file when prmtop does not exist")
	parser.add_argument("-p", dest = "prmtop", metavar = "OUTPUT.prmtop", required = True, help = "Amber topology file")
	parser.add_argument("-sc", dest = "temp_dir", metavar = "TEMP_DIR", default = ".", help = "Temporary directory (Default: current dir)")

	gmx_option = parser.add_argument_group("gromacs option")
	gmx_option.add_argument("-b", dest = "begin", metavar = "START_TIME", type = int, help = "First frame (ps) to read from trajectory (start from 0)")
	gmx_option.add_argument("-e", dest = "end", metavar = "END_TIME", type = int, help = "Last frame (ps) to read from trajectory (start from 0)")
	gmx_option.add_argument("--offset", dest = "offset", metavar = "OFFSET", type = int, default = 1, help = "Output interval (Default: 1)")

	cpptraj_option = parser.add_argument_group("cpptraj option")
	cpptraj_option.add_argument("-mc", dest = "center_mask", metavar = "CENTER_MASK", help = "center mask for cpptraj")
	cpptraj_option.add_argument("-ms", dest = "strip_mask", metavar = "STRIP_MASK", help = "strip mask for cpptraj")

	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")

	args = parser.parse_args()

	check_exist(args.tpr, 2)
	for trj_file in args.trr:
		check_exist(trj_file, 2)
	check_exist(args.top, 2)


	command_gmx = check_command("gmx")

	# 一時ファイル名の決定
	tempfile_name = ""
	with tempfile.NamedTemporaryFile(mode = "w", prefix = ".trr2nc_", dir = args.temp_dir) as obj_temp:
		tempfile_name = obj_temp.name


	# トポロジーファイルの読み込み
	obj_top = MoleculeTopology(args.top)


	# 周期境界でジャンプしないトラジェクトリの作成
	temp_traj1 = tempfile_name + "1.trr"
	temp_ndx = tempfile_name + ".ndx"
	sys.stderr.write(colored("Creating nojump trajectory ({file})\n".format(file = temp_traj1), "red", attrs = ["bold"]))

	tpr = args.tpr
	trajectories = " ".join(args.trr)
	command = "{command} trjconv -s {tpr} -f {trajectory} -o {output} -pbc nojump -skip {offset}".format(command = command_gmx, tpr = tpr, trajectory = trajectories, output = temp_traj1, offset = args.offset)
	if args.begin is not None:
		command += " -b {0}".format(args.begin)
	if args.end is not None:
		command += " -e {0}".format(args.end)

	obj_ndx = NDXFile(obj_top)
	if args.strip_mask is None:
		# strip_mask が指定されていない場合
		obj_ndx.output_ndx(temp_ndx)
		command += " -n {0} << 'EOF'\n0\nEOF".format(temp_ndx)
	else:
		# strip_mask が指定されていた場合
		# ndx ファイルの作成
		obj_ndx.add_def("Strip", "!" + args.strip_mask).output_ndx(temp_ndx)
		obj_ndx.del_def(1)
		command += " -n {0} << 'EOF'\n1\nEOF".format(temp_ndx)


		# 次の中央配置のための準備
		# grompp 用の構造ファイルの作成 (strip 済み)
		ref_coord = tempfile_name + "_ref.gro"
		command_sub = "{command} trjconv -s {input_tpr} -f {trajectory} -o {ref_coord} -n {ndx} -b 1 -e 1 << 'EOF'\n1\nEOF".format(command = command_gmx, trajectory = trajectories, input_tpr = args.tpr, ndx = temp_ndx, ref_coord = ref_coord)
		exec_sp(command_sub, False)

		# strip したトポロジーの作成
		temp_top = tempfile_name + ".top"
		obj_top.set_mask("!" + args.strip_mask, True)
		obj_top.save_file(temp_top)

		# ndx ファイルの更新
		obj_ndx = NDXFile(obj_top)

		# grompp 用の mdp ファイルの作成
		temp_mdp1 = tempfile_name + "1.mdp"
		temp_mdp2 = tempfile_name + "2.mdp"
		with open(temp_mdp1, "w") as obj_output:
			obj_output.write("integrator = steep\n emtol = 1000.0\n emstep = 0.01\n nsteps = 50000\n nstlist = 100\n ns_type = grid\n rlist = 1.0\n coulombtype = PME\n rcoulomb = 1.0\n nstlog = 1 \n pbc = xyz\n vdwtype = cut-off\n constraints = none\n cutoff-scheme = Verlet\n")

		# tpr の作成
		temp_tpr = tempfile_name + ".tpr"
		command_sub = "{command} grompp -f {temp_mdp1} -c {ref_coord} -p {temp_top} -o {temp_tpr} -po {temp_mdp2} -maxwarn 10".format(command = command_gmx, temp_mdp1 = temp_mdp1, ref_coord = ref_coord, temp_top = temp_top, temp_tpr = temp_tpr, temp_mdp2 = temp_mdp2)
		exec_sp(command_sub, False)

		# 後処理
		tpr = temp_tpr
		os.remove(ref_coord)
		os.remove(temp_top)
		os.remove(temp_mdp1)
		os.remove(temp_mdp2)

	exec_sp(command, True)


	# 分子を中央に配置したトラジェクトリの作成
	temp_traj2 = tempfile_name + "2.trr"
	sys.stderr.write(colored("INFO: Creating centered trajectory ({file})\n".format(file = temp_traj2), "red", attrs = ["bold"]))
	trajectories = " ".join(args.trr)
	command = "{command} trjconv -s {tpr} -f {trajectory} -o {output} -pbc res -ur compact".format(command = command_gmx, tpr = tpr, trajectory = temp_traj1, output = temp_traj2)


	if args.center_mask is None:
		command += " << 'EOF'\n0\nEOF"
	else:
		obj_ndx = NDXFile(obj_top).add_def("Center", args.center_mask).output_ndx(temp_ndx)
		command += " -center -n {ndx} << 'EOF'\n1\n0\nEOF".format(ndx = temp_ndx)
	exec_sp(command, True)

	os.remove(temp_traj1)
	if args.strip_mask is not None:
		os.remove(tpr)
	os.remove(temp_ndx)


	# prmtop の作成
	sys.stderr.write(colored("INFO: Creating prmtop ({file})\n".format(file = args.prmtop), "red", attrs = ["bold"]))
	if args.flag_overwrite == False:
		check_overwrite(prmtop)
	obj_top.set_mask("*").save_file(args.prmtop)


	# nc ファイルに変換
	if args.flag_overwrite == False:
		check_overwrite(args.nc)

	temp_in = tempfile_name + ".in"
	sys.stderr.write(colored("INFO: Converting AMBER trajectory ({file})\n".format(file = args.nc), "red", attrs = ["bold"]))
	with open(temp_in, "w") as obj_output:
		obj_output.write("parm {0}\n".format(args.prmtop))
		obj_output.write("trajin {0}\n".format(temp_traj2))
		obj_output.write("trajout {0}\n".format(args.nc))
		obj_output.write("go\n")

	command_cpptraj = check_command("cpptraj")
	exec_sp("{0} -i {1}".format(command_cpptraj, temp_in), True)

	os.remove(temp_traj2)
	os.remove(temp_in)
