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
	tempfile_name = ""
	with tempfile.NamedTemporaryFile(mode = "w", prefix = ".trr2nc_", dir = ".") as obj_temp:
		tempfile_name = obj_temp.name

	command_parmed = check_command("parmed")
	temp_in = tempfile_name + ".in"
	with open(temp_in, "w") as obj_output:
		obj_output.write("gromber {0}\n".format(top))
		obj_output.write("outparm {0}\n".format(prmtop))
		obj_output.write("quit\n")
	exec_sp("{0} -n < {1}".format(command_parmed, temp_in), True)
	os.remove(temp_in)


def convert_trajectory(tpr, trr, begin, end, prmtop, strip_mask, fitting_mask, output):
	command_gmx = check_command("gmx")

	tempfile_name = ""
	with tempfile.NamedTemporaryFile(mode = "w", prefix = ".trr2nc_", dir = ".") as obj_temp:
		tempfile_name = obj_temp.name

	# 周期境界条件でジャンプしないトラジェクトリの作成
	temp_traj1 = tempfile_name + "1.trr"
	sys.stderr.write("{start}Creating nojumped trajectory ({file}){end}\n".format(file = temp_traj1, start = basic.color.LRED + basic.color.BOLD, end = basic.color.END))
	trajectories = " ".join(trr)
	command = "{0} trjconv -s {1} -f {2} -o {3} -pbc nojump".format(command_gmx, tpr, trajectories, temp_traj1)
	if begin is not None:
		command += " -b {0}".format(begin)
	if end is not None:
		command += " -e {0}".format(end)
	command += " << 'EOF'\n0\nEOF"
	exec_sp(command, True)

	# 特定分子を中央に配置したトラジェクトリの作成
	temp_traj2 = tempfile_name + "2.trr"
	sys.stderr.write("{start}Creating solute centered trajectory ({file}){end}\n".format(file = temp_traj2, start = basic.color.LRED + basic.color.BOLD, end = basic.color.END))
	exec_sp("{0} trjconv -s {1} -f {2} -o {3} -pbc mol -center -ur compact << 'EOF'\n0\n0\nEOF".format(command_gmx, tpr, temp_traj1, temp_traj2), True)
	os.remove(temp_traj1)

	# nc ファイルに変換
	temp_in = tempfile_name + ".in"
	sys.stderr.write("{start}Converting AMBER trajectory ({file}){end}\n".format(file = output, start = basic.color.LRED + basic.color.BOLD, end = basic.color.END))
	with open(temp_in, "w") as obj_output:
		obj_output.write("parm {0}\n".format(prmtop))
		obj_output.write("trajin {0}\n".format(temp_traj2))
		if strip_mask is not None:
			obj_output.write("strip {0}\n".format(strip_mask))
		if fitting_mask is not None:
			obj_output.write("rms {0} first mass\n".format(fitting_mask))
		obj_output.write("trajout {0}\n".format(output))
		obj_output.write("go\n")
		if strip_mask is not None:
			obj_output.write("clear trajin\n")
			obj_output.write("parmstrip {0}\n".format(strip_mask))
			obj_output.write("parmwrite out {0}\n".format(prmtop))
			obj_output.write("go\n")

	command_cpptraj = check_command("cpptraj")
	exec_sp("{0} -i {1}".format(command_cpptraj, temp_in), True)

	os.remove(temp_traj2)
	os.remove(temp_in)


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "trr2nc.py - Convert trr to nc with treating PBC", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-s", dest = "tpr", metavar = "INPUT.tpr", required = True, help = "Gromacs run input file")
	parser.add_argument("-x", dest = "trr", metavar = "INPUT.<trr|xtc|gro>", nargs = "+", help = "Gromacs trajectory file")
	parser.add_argument("-o", dest = "nc", metavar = "OUTPUT.<nc|mdcrd>", required = True, help = "output for Amber trajectory (.nc)")
	parser.add_argument("-p", dest = "prmtop", metavar = "INPUT.prmtop", required = True, help = "Amber topology file")

	gmx_option = parser.add_argument_group("gromacs option")
	gmx_option.add_argument("-b", dest = "begin", metavar = "START_TIME", type = int, help = "First frame (ps) to read from trajectory")
	gmx_option.add_argument("-e", dest = "end", metavar = "END_TIME", type = int, help = "Last frame (ps) to read from trajectory")
	gmx_option.add_argument("-t", dest = "top", metavar = "INPUT.top", help = "Gromacs topology file when prmtop does not exist")

	cpptraj_option = parser.add_argument_group("cpptraj option")
	cpptraj_option.add_argument("-mf", dest = "fitting_mask", metavar = "FITTING_MASK", help = "fitting mask for cpptraj")
	cpptraj_option.add_argument("-ms", dest = "strip_mask", metavar = "STRIP_MASK", help = "strip mask for cpptraj")

	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")

	args = parser.parse_args()

	basic.check_exist(args.tpr, 2)
	for trj_file in args.trr:
		basic.check_exist(trj_file, 2)

	if os.path.exists(args.prmtop) == False:
		# prmtop がない場合
		if args.top is None:
			# 力場ファイルが指定されていない場合
			sys.stderr.write(" ERROR: -t option is not specified when prmtop does not exist\n")
			sys.exit(1)

		make_prmtop(args.top, args.prmtop, args.strip_mask)

	convert_trajectory(args.tpr, args.trr, args.begin, args.end, args.prmtop, args.strip_mask, args.fitting_mask, args.nc)
