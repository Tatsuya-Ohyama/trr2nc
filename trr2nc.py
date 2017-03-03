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
# check_command (コマンドの有無を確認)
def check_command(command_name):
	import subprocess

	process = subprocess.Popen("which %s" % command_name, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	(stdout, stderr) = process.communicate()

	stdout = stdout.decode("utf-8")
	if stdout == "":
		sys.stderr.write("ERROR: Invalid command (%s).\n       Check PATH or installed in the system.\n" % command_name)
		sys.exit(1)
	command_path = stdout.rstrip("\r\n")

	return command_path


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "trr2nc.py - Convert trr to nc with treating PBC", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-s", dest = "tpr", metavar = "tpr", required = True, help = "Gromacs topology file (.tpr)")

	parser.add_argument("-x", dest = "trr", metavar = "trajectory", required = True, help = "Gromacs trajectory file (.trr, .xtc)")

	parser.add_argument("-p", dest = "prmtop", metavar = "prmtop", help = "Amber topology file (.prmtop)")

	parser.add_argument("-m", dest = "mask", metavar = "mask", required = True, help = "fitting mask for cpptraj")
	parser.add_argument("-gc", dest = "group_center", metavar = "center_of_group" , required = True, type = int, help = "center of group")
	parser.add_argument("-sm", dest = "strip_mask", metavar = "strip_mask", help = "strip mask for cpptraj")
	parser.add_argument("-o", dest = "nc", metavar = "nc", required = True, help = "output for Amber trajectory (.nc)")
	parser.add_argument("-n", dest = "ndx", metavar = "ndx", help = "index file for Gromacs (.ndx)")
	parser.add_argument("-b", dest = "begin", metavar = "startTime", type = int, help = "First frame (ps) to read from trajectory")
	parser.add_argument("-e", dest = "end", metavar = "startTime", type = int, help = "Last frame (ps) to read from trajectory")
	parser.add_argument("-skip", dest = "skip", metavar = "skipFrame", type = int, help = "Only write every nr-th frame")
	args = parser.parse_args()


# メモ 2016/12/12
# 1. cpptraj で prmtop を読み込み
# 2. アウトプットで指定された prmtop を作成
# 3. 中心に合わせた prmtop を作成
# 4. prmtop を読み込み、原子順序番号が列挙された ndx ファイルをプログラム内で作成
# 5. 作成した ndx ファイルを gmx に与えて実行


"""

	# 必須ファイルの確認
	check_file(args.tpr)
	check_file(args.trr)
	check_file(args.prmtop)
	if args.ndx != None:
		check_file(args.ndx)

	# コマンドの確認
	path_gmx = check_command("gmx")
	path_cpptraj = check_command("cpptraj")

	# 一時ファイル名取得
	obj_temp = tempfile.NamedTemporaryFile(mode = "w", prefix = ".trr2nc_", delete = True, dir = ".")
	random_name = obj_temp.name
	obj_temp.close()
	tempfile_prmtop = random_name + "_prmtop.prmtop"
	tempfile_nojump = random_name + "_nojump.trr"
	tempfile_center = random_name + "_center.trr"

	# strip された prmtop 作成
	prmtop = args.prmtop
	if args.strip_mask != None:
		tempfile_in = random_name + "_prmtop.prmtop"
		cpptraj_infile = [
			"parm %s" % args.prmtop,
			"parmstrip %s" % args.strip_mask,
			"parmwrite out %s" % tempfile_prmtop,
			"go"
		]
		command = "%s <<'INPUT'\n%s\nINPUT" % (path_cpptraj, "\n".join(cpptraj_infile))
		process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		while True:
			# エラーログ取得
			line = process.stderr.readline().decode("utf-8").rstrip("\r\n")
			if line:
				log.append(line)
			if not line and process.poll() is not None:
				break
		process.wait()
		if process.returncode == 0:
			# コマンド成功
			sys.stderr.write("done\n")
		else:
			# コマンド失敗
			message = "\nSomething wrong in follow command\n %s\n" % command
			sys.stderr.write(message)
			sys.stderr.write("Check follow STDERR:\n")
			for line in log:
				line = line + "\n"
				sys.stderr.write(line)
			sys.exit(1)
		prmtop = tempfile_prmtop


	# 周期境界でジャンプしないトラジェクトリの作成
	log = []
	sys.stderr.write("Create nojump trajectory ... ")
	sys.stderr.flush()
	command = "%s trjconv -s %s -f %s -o %s -pbc nojump" % (path_gmx, args.tpr, args.trr, tempfile_nojump)
	if args.begin != None:
		command += " -b %d" % args.begin
	if args.end != None:
		command += " -e %d" % args.end
	if args.skip != None:
		command += " -skip %d" % args.skip
	command += " <<'INPUT'\n%d\nINPUT" % args.group_center
	process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	while True:
		# エラーログ取得
		line = process.stderr.readline().decode("utf-8").rstrip("\r\n")
		if line:
			log.append(line)
		if not line and process.poll() is not None:
			break
	process.wait()
	if process.returncode == 0:
		# コマンド成功
		sys.stderr.write("done\n")
	else:
		# コマンド失敗
		message = "\nSomething wrong in follow command\n %s\n" % command
		sys.stderr.write(message)
		sys.stderr.write("Check follow STDERR:\n")
		for line in log:
			line = line + "\n"
			sys.stderr.write(line)
		sys.exit(1)


	# ターゲット分子を中心としたトラジェクトリの作成
	log = []
	sys.stderr.write("Create molecular centered trajectory ... ")
	sys.stderr.flush()
	command = "%s trjconv -s %s -f %s -o %s -center -pbc mol -ur compact " % (path_gmx, args.tpr, tempfile_nojump, tempfile_center)
	if args.ndx != None:
		command += "-n %s " % args.ndx
	command += "<<'INPUT'\n%d\n0\nINPUT" % (args.group_center)
	process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	while True:
		# エラーログ取得
		line = process.stderr.readline().decode("utf-8").rstrip("\r\n")
		if line:
			log.append(line)
		if not line and process.poll() is not None:
			break
	process.wait()
	if process.returncode == 0:
		# コマンド成功
		sys.stderr.write("done\n")
	else:
		# コマンド失敗
		message = "\nSomething wrong in follow command\n %s\n" % command
		sys.stderr.write(message)
		sys.stderr.write("Check follow STDERR:\n")
		for line in log:
			line = line + "\n"
			sys.stderr.write(line)
		sys.exit(1)

	# AMBER 用トラジェクトリ作成
	log = []
	sys.stderr.write("Create AMBER trajectory ... ")
	sys.stderr.flush()
	cpptraj_infile = [
		"parm %s" % args.prmtop,
		"trajin %s" % tempfile_center,
		"autoimage",
		"rms %s first mass" % args.mask,
		"trajout %s" % args.nc,
		"go"
	]
	command = "%s <<'INPUT'\n%s\nINPUT" % (path_cpptraj, "\n".join(cpptraj_infile))
	process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	while True:
		# エラーログ取得
		line = process.stderr.readline().decode("utf-8").rstrip("\r\n")
		if line:
			log.append(line)
		if not line and process.poll() is not None:
			break
	process.wait()
	if process.returncode == 0:
		# コマンド成功
		sys.stderr.write("done\n")
	else:
		# コマンド失敗
		message = "\nSomething wrong in follow command\n %s\n" % command
		sys.stderr.write(message)
		sys.stderr.write("Check follow STDERR:\n")
		for line in log:
			line = line + "\n"
			sys.stderr.write(line)
		sys.exit(1)

	os.remove(tempfile_nojump)
	os.remove(tempfile_center)
"""
