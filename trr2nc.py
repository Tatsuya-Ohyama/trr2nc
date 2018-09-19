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
import parmed

from basic_func import check_exist, check_overwrite
from classes.ndx_file import NDXFile

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


def convert_trajectory(top, tpr, trr, prmtop, output, flag_overwrite, center_mask = None, strip_mask = None, begin = None, end = None):
	""" トラジェクトリを trr から nc に変換する関数 """
	command_gmx = check_command("gmx")

	# 一時ファイル名の決定
	tempfile_name = ""
	with tempfile.NamedTemporaryFile(mode = "w", prefix = ".trr2nc_", dir = args.temp_dir) as obj_temp:
		tempfile_name = obj_temp.name

	# トポロジーファイルの読み込み
	obj_top = parmed.gromacs.GromacsTopologyFile(top)


	# 周期境界でジャンプしないトラジェクトリの作成
	temp_traj1 = tempfile_name + "1.trr"
	temp_ndx = tempfile_name + ".ndx"
	sys.stderr.write(colored("Creating nojump trajectory ({file})\n".format(file = temp_traj1), "red", attrs = ["bold"]))

	trajectories = " ".join(trr)
	command = "{command} trjconv -s {tpr} -f {trajectory} -o {output} -pbc nojump".format(command = command_gmx, tpr = tpr, trajectory = trajectories, output = temp_traj1)
	if begin is not None:
		command += " -b {0}".format(begin)
	if end is not None:
		command += " -e {0}".format(end)
	if strip_mask is None:
		# strip_mask が指定されていない場合
		make_ndx(temp_ndx, [x.idx for x in obj_top.atoms])
		command += " -n {0} << 'EOF'\n0\nEOF".format(temp_ndx)
	else:
		# strip_mask が指定されていた場合
		amber_mask = parmed.amber.AmberMask(obj_top, strip_mask)
		make_ndx(temp_ndx, [x.idx for x in obj_top.atoms], None, [x for x in amber_mask.Selected(invert = True)])
		command += " -n {0} << 'EOF'\n1\nEOF".format(temp_ndx)

		# grompp 用の構造ファイルの作成 (strip 済み)
		ref_coord = tempfile_name + "_ref.gro"
		command_sub = "{command} trjconv -s {input_tpr} -f {trajectory} -o {ref_coord} -n {ndx} -b 1 -e 1 << 'EOF'\n1\nEOF".format(command = command_gmx, trajectory = trajectories, input_tpr = tpr, ndx = temp_ndx, ref_coord = ref_coord)
		exec_sp(command_sub, False)

		# strip したトポロジーの作成
		temp_top = tempfile_name + ".top"
		obj_top.strip(strip_mask)
		obj_top.write(temp_top)

		# grompp 用の mdp ファイルの作成
		temp_mdp1 = tempfile_name + "1.mdp"
		temp_mdp2 = tempfile_name + "2.mdp"
		with open(temp_mdp1, "w") as obj_output:
			obj_output.write("integrator = steep\n emtol = 1000.0\n emstep = 0.01\n nsteps = 50000\n nstlist = 100\n ns_type = grid\n rlist = 1.0\n coulombtype = PME\n rcoulomb = 1.0\n nstlog = 1 \n pbc = xyz\n vdwtype = cut-off\n constraints = none\n cutoff-scheme = Verlet\n")

		# tpr の作成
		temp_tpr = tempfile_name + ".tpr"
		command_sub = "{command} grompp -f {temp_mdp1} -c {ref_coord} -p {temp_top} -o {temp_tpr} -po {temp_mdp2} -maxwarn 10".format(command = command_gmx, temp_mdp1 = temp_mdp1, ref_coord = ref_coord, temp_top = temp_top, temp_tpr = temp_tpr, temp_mdp2 = temp_mdp2)
		exec_sp(command_sub, True)

		# 後処理
		tpr = temp_tpr
		os.remove(ref_coord)
		os.remove(temp_top)
		os.remove(temp_mdp1)
		os.remove(temp_mdp2)

	exec_sp(command, False)

	# 分子を中央に配置したトラジェクトリの作成
	temp_traj2 = tempfile_name + "2.trr"
	sys.stderr.write(colored("Creating centered trajectory ({file})\n".format(file = temp_traj2), "red", attrs = ["bold"]))
	trajectories = " ".join(trr)
	command = "{command} trjconv -s {tpr} -f {trajectory} -o {output} -pbc res -ur compact".format(command = command_gmx, tpr = tpr, trajectory = temp_traj1, output = temp_traj2)
	if begin is not None:
		command += " -b {0}".format(begin)
	if end is not None:
		command += " -e {0}".format(end)
	if center_mask is None:
		command += " << 'EOF'\n0\nEOF"
	else:
		amber_mask = parmed.amber.AmberMask(obj_top, center_mask)
		make_ndx(temp_ndx, [x.idx for x in obj_top.atoms], [x for x in amber_mask.Selected()], None)
		command += " -center -n {ndx} << 'EOF'\n1\n0\nEOF".format(ndx = temp_ndx)
	exec_sp(command, True)
	os.remove(temp_traj1)
	if strip_mask is not None:
		os.remove(tpr)
	os.remove(temp_ndx)

	# prmtop の作成
	sys.stderr.write(colored("Creating prmtop ({file})\n".format(file = prmtop), "red", attrs = ["bold"]))
	if flag_overwrite == False:
		check_overwrite(prmtop)

	amber_top = parmed.amber.AmberParm.from_structure(obj_top)
	amber_top.write_parm(prmtop)

	# nc ファイルに変換
	temp_in = tempfile_name + ".in"
	sys.stderr.write(colored("{start}Converting AMBER trajectory ({file}){end}\n".format(file = output), "red", attrs = ["bold"]))
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

	check_exist(args.tpr, 2)
	for trj_file in args.trr:
		check_exist(trj_file, 2)
	check_exist(args.top, 2)

	if args.flag_overwrite == False:
		check_overwrite(args.prmtop)

	convert_trajectory(args.top, args.tpr, args.trr, args.prmtop, args.nc, args.flag_overwrite, args.center_mask, args.strip_mask, args.begin, args.end)
