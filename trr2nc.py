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


def pdbfix(input, output):
	""" trjconv によって生成された PDB を tleap で読み込めるようにする関数 """
	re_atom = re.compile(r"^(?:(?:ATOM)|(?:HETATM))")
	with tempfile.TemporaryFile(mode = "r+") as obj_temp:
		with open(input, "r") as obj_input:
			for line in obj_input:
				if re_atom.search(line):
					line = line[0:54]
					atomname = line[12:14].strip()
					atomtype = line[12:16]
					residuename = line[17:20]

					if residuename == "SOL":
						# 水分子
						if "OW" in atomtype:
							line = line[0 : 12] + " O  " + line[16 : 17] + "WAT" + line[20:]
						elif "HW" in atomtype:
							if "HW1" in atomtype:
								line = line[0 : 12] + " H1 " + line[16 : 17] + "WAT" + line[20:]
							elif "HW2" in atomtype:
								line = line[0 : 12] + " H2 " + line[16 : 17] + "WAT" + line[20:]

					elif "K" in atomname:
						line = line[0 : 12] + " K+ " + line[16 : 17] + " K+" + line[20:]

					elif "CL" in atomname:
						line = line[0 : 12] + "Cl- " + line[16 : 17] + "Cl-" + line[20:]

					line = "{0}\n".format(line)

				obj_temp.write(line)

		obj_temp.seek(0)
		with open(output, "w") as obj_output:
			for line in obj_temp:
				obj_output.write(line)


def exec_sp(command):
	""" subprocess で外部プログラムを実行する関数 """
	process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	streamdata = process.communicate()

	if process.returncode == 1:
		sys.stderr.write("ERROR: subprocess failed in make_prmtop():\n    '{0}'.\n".format(command))
		sys.exit(1)


def make_prmtop(tpr, trr, forcefields, output):
	""" prmtop を作成する関数 """
	command_tleap = check_command("tleap")
	command_gmx = check_command("gmx")

	tempfile_name = ""
	with tempfile.NamedTemporaryFile(mode = "w", prefix = ".trr2nc_", delete = True, dir = ".") as obj_temp:
		tempfile_name = obj_temp.name

	tmp_gro = tempfile_name + ".gro"
	tmp_pdb = tempfile_name + ".pdb"
	tmp_in = tempfile_name + ".in"
	tmp_inpcrd = tempfile_name + ".inpcrd"
	exec_sp("{0} trjconv -s {1} -f {2} -o {3} -pbc nojump <<'EOF'\n0\nEOF".format(command_gmx, tpr, trr, tmp_gro))
	exec_sp("{0} trjconv -s {1} -f {2} -o {3} -pbc mol -center -ur compact << 'EOF'\n1\n0\nEOF".format(command_gmx, tpr, tmp_gro, tmp_pdb))
	pdbfix(tmp_pdb, tmp_pdb)
	with open(tmp_in, "w") as obj_output:
		for forcefield in forcefields:
			obj_output.write("source {0}\n".format(forcefield))
		obj_output.write("complex = loadpdb {0}\n".format(tmp_pdb))
		obj_output.write("saveamberparm complex {0} {1}\n".format(output, tmp_inpcrd))
		obj_output.write("quit\n")
	exec_sp("{0} -f {1} << 'EOF'\nquit\nEOF".format(command_tleap, tmp_in))

	os.remove(tmp_gro)
	os.remove(tmp_pdb)
	os.remove(tmp_in)
	os.remove(tmp_inpcrd)


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "trr2nc.py - Convert trr to nc with treating PBC", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-s", dest = "tpr", metavar = "INPUT.tpr", required = True, help = "Gromacs topology file")
	parser.add_argument("-x", dest = "trr", metavar = "INPUT.<trr|xtc>", required = True, help = "Gromacs trajectory file for input")

	parser.add_argument("-pi", dest = "prmtop_input", metavar = "INPUT.prmtop", help = "Amber topology file")

	parser.add_argument("-po", dest = "prmtop_output", metavar = "OUTPUT.prmtop", help = "Amber topology file for input")
	parser.add_argument("-ff", dest = "forcefield", metavar = "FF_FILE", nargs = "*", default = [], help = "Force field files of AMBER when not specify -pi")
	# parser.add_argument("-m", dest = "mask", metavar = "mask", required = True, help = "fitting mask for cpptraj")
	# parser.add_argument("-gc", dest = "group_center", metavar = "center_of_group" , required = True, type = int, help = "center of group")
	# parser.add_argument("-sm", dest = "strip_mask", metavar = "strip_mask", help = "strip mask for cpptraj")
	# parser.add_argument("-o", dest = "nc", metavar = "nc", required = True, help = "output for Amber trajectory (.nc)")
	# parser.add_argument("-n", dest = "ndx", metavar = "ndx", help = "index file for Gromacs (.ndx)")
	# parser.add_argument("-b", dest = "begin", metavar = "startTime", type = int, help = "First frame (ps) to read from trajectory")
	# parser.add_argument("-e", dest = "end", metavar = "startTime", type = int, help = "Last frame (ps) to read from trajectory")
	# parser.add_argument("-skip", dest = "skip", metavar = "skipFrame", type = int, help = "Only write every nr-th frame")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	args = parser.parse_args()

	basic.check_exist(args.tpr, 2)
	basic.check_exist(args.trr, 2)

	if args.prmtop_input is None:
		# prmtop がない場合

		if len(args.forcefield) == 0:
			# 力場ファイルが指定されていない場合
			sys.stderr.write(" ERROR: -ff option is not specified when no -pi option\n")
			sys.exit(1)

		if args.prmtop_output is None:
			sys.stderr.write(" ERROR: -po option does not specified when no -pi option\n")
			sys.exit(1)

		if args.flag_overwrite == False:
			basic.check_overwrite(args.prmtop_output, args.forcefield)
		make_prmtop(args.tpr, args.trr, args.forcefield, args.prmtop_output)


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
