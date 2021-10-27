#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
trr2nc
Program to convert Gromacs trajectory to AMBER trajectory
"""

import sys
import signal
import os
import argparse
import subprocess
import tempfile
from termcolor import colored
import re
import parmed

from mods.func_prompt_io import *
from mods.file_NDX import FileNDX


global delete_files
def delete_all(signal=0):
	"""
	Function to delete temporary files

	Args:
		signal (int, optional): signal number (Default: 0)
	"""
	for file in delete_files:
		os.remove(file)

	if signal in [2, 15]:
		# ctrl-c
		sys.exit(1)
	else:
		sys.exit(0)
signal.signal(signal.SIGINT, delete_all)



# =============== constant =============== #
COMMAND_NAME_GMX = "gmx"
COMMAND_NAME_CPPTRAJ = "cpptraj"
LOG_COLOR = "yellow"
RE_CPPTRAJ_VER = re.compile(r"CPPTRAJ: Version (V.+?) \(AmberTools (V.+?)\)")



# =============== functions =============== #
def check_command(command_name):
	"""
	Function to check for the presence of a command

	Args:
		command_name (str): command name

	Returns:
		str: absolute path for command
	"""
	process = subprocess.Popen(
		"which {0}".format(command_name),
		shell=True,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE
	)
	(stdout, stderr) = process.communicate()

	stdout = stdout.decode("utf-8")
	if stdout == "":
		sys.stderr.write("ERROR: Invalid command ({0}).\n       Check PATH or installed in the system.\n".format(command_name))
		sys.exit(1)
	command_path = stdout.rstrip("\r\n")

	return command_path


def get_cpptraj_ver(cpptraj_path):
	"""
	Function to get cpptraj and AmberTools versions

	Args:
		cpptraj_path (str): cpptraj path

	Returns:
		list: [cpptraj_ver, AmberTools_ver]
	"""
	obj_process = subprocess.run([cpptraj_path, "--version"], stdout=subprocess.PIPE, text=True)
	obj_match = RE_CPPTRAJ_VER.search(obj_process.stdout)
	return [obj_match.group(1), obj_match.group(2)]


def exec_sp(command, operation=False):
	"""
	Function to execute outer program by subprocess module

	Args:
		command (str): command line
		operation (bool, optional): show prompt (Default: False)
	"""
	if operation:
		process = subprocess.Popen(
			command,
			shell=True
		)
	else:
		process = subprocess.Popen(
			command,
			shell=True,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE
		)
	streamdata = process.communicate()

	if process.returncode == 1:
		sys.stderr.write("ERROR: subprocess failed\n    '{0}'.\n".format(command))
		sys.exit(1)



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "Program to convert Gromacs trajectory to AMBER trajectory", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-s", dest="TPR_FILE", metavar="INPUT.tpr", required=True, help="Gromacs run input file")
	parser.add_argument("-x", dest="TRAJECTORY_FILE", metavar="INPUT.<trr|xtc|gro>", required=True, help="Gromacs trajectory file")
	parser.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.<nc|mdcrd|xtc>", required=True, help="output trajectory")
	parser.add_argument("-t", dest="TOP_FILE", metavar="INPUT.top", required=True, help="Gromacs topology file")
	parser.add_argument("-p", dest="PRMTOP_FILE", metavar="OUTPUT.prmtop", required=True, help="Amber topology file")
	parser.add_argument("-sc", dest="TEMP_DIR", metavar="TEMP_DIR", default=".", help="Temporary directory (Default: current dir)")
	parser.add_argument("--separate-mol", dest="SEPARATE_MOL", metavar="MOL_NAME", nargs="+", default=[], help="separate molecules into individual molecules (specify molecule name written in .top file) (periodic boundary condition problem)")

	gmx_option = parser.add_argument_group("gromacs option")
	gmx_option.add_argument("-b", dest="BEGIN", metavar="START_TIME", type=int, help="First frame index to read from trajectory (start from 0)")
	gmx_option.add_argument("-e", dest="END", metavar="END_TIME", type=int, help="Last frame index to read from trajectory (start from 0)")
	gmx_option.add_argument("--offset", dest="OFFSET", metavar="OFFSET", type=int, default=1, help="Output interval (Default: 1)")
	gmx_option.add_argument("--gmx", dest="COMMAND_GMX", metavar="COMMAND_GMX", help="command line path for `gmx` (Default: autodetect)")

	cpptraj_option = parser.add_argument_group("cpptraj option")
	cpptraj_option.add_argument("-mc", dest="CENTER_MASK", metavar="CENTER_MASK", required=True, help="center mask for cpptraj")
	cpptraj_option.add_argument("-ms", dest="STRIP_MASK", metavar="STRIP_MASK", help="strip mask for cpptraj")
	cpptraj_option.add_argument("--cpptraj", dest="COMMAND_CPPTRAJ", metavar="COMMAND_CPPTRAJ", help="command line path for `cpptraj` (Default: autodetect)")

	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")

	args = parser.parse_args()

	# check arguments
	check_exist(args.TPR_FILE, 2)
	check_exist(args.TRAJECTORY_FILE, 2)
	check_exist(args.TOP_FILE, 2)
	command_gmx = args.COMMAND_GMX
	if command_gmx is None:
		command_gmx = check_command(COMMAND_NAME_GMX)

	command_cpptraj = args.COMMAND_CPPTRAJ
	if command_cpptraj is None:
		command_cpptraj = check_command(COMMAND_NAME_CPPTRAJ)

	cpptraj_vers = get_cpptraj_ver(command_cpptraj)
	ambertools_ver = int(cpptraj_vers[1][1:].split(".")[0])

	if os.path.splitext(args.OUTPUT_FILE)[1].lower() == ".xtc":
		if ambertools_ver <= 16:
			sys.stderr.write("ERROR: output of .xtc file is only supported in AmberTools version 17.0 or later.\n")
			sys.exit(1)

	# determine name of temporary file
	tempfile_name = ""
	tempfile_name_full = ""
	with tempfile.NamedTemporaryFile(mode="w", prefix=".trr2nc_", dir=".") as obj_output:
		tempfile_name = os.path.basename(obj_output.name)
	tempfile_name_full = os.path.join(args.TEMP_DIR, tempfile_name)
	delete_files = []

	process_i = 0
	max_process = None
	if os.path.splitext(args.TRAJECTORY_FILE)[1].lower() == ".nc" \
		and os.path.splitext(args.OUTPUT_FILE)[1].lower() == ".nc":
		max_process = 3
	else:
		max_process = 6

	# read topology file
	process_i += 1
	sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Loading topology file"), LOG_COLOR, attrs=["bold"]))
	obj_topol = parmed.gromacs.GromacsTopologyFile(args.TOP_FILE)


	# create .prmtop file
	process_i += 1
	sys.stdout.write(colored("Process ({0}/{1}): {2} => {3}\n".format(process_i, max_process, "Generate prmtop", args.PRMTOP_FILE), LOG_COLOR, attrs=["bold"]))
	check_overwrite(args.PRMTOP_FILE, args.FLAG_OVERWRITE)
	if args.STRIP_MASK is not None:
		obj_topol.strip(args.STRIP_MASK)
	obj_topol.save(args.PRMTOP_FILE)


	# conversion of trajectory file
	trajectory_input = args.TRAJECTORY_FILE
	ndx_file = None
	if os.path.splitext(args.TRAJECTORY_FILE)[1] in [".gro", ".xtc"]:
		# treat PBC
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate trajectory with adjusted molecules across the boundary."), LOG_COLOR, attrs=["bold"]))
		gmx_arg = {
			"-s": args.TPR_FILE,
			"-f": args.TRAJECTORY_FILE,
			"-o": None,
			"-b": args.BEGIN,
			"-e": args.END,
			"-skip": args.OFFSET,
			"-pbc": "whole",
			"-n": None,
		}
		gmx_eof = None

		gmx_arg["-o"] = tempfile_name_full + "_step1_whole.trr"
		if ambertools_ver >= 17:
			gmx_arg["-o"] = tempfile_name_full + "_step1_whole.xtc"
		delete_files.append(gmx_arg["-o"])

		obj_ndx = FileNDX(obj_topol)
		obj_ndx.add_def("Center", args.CENTER_MASK)
		ndx_file = tempfile_name_full + ".ndx"
		gmx_arg["-n"] = ndx_file
		delete_files.append(gmx_arg["-n"])
		if args.STRIP_MASK is None:
			# no strip_mask
			gmx_eof = "<< 'EOF'\nSystem\nEOF"
		else:
			# when strip_mask is specified, create .ndx file
			mask = "!({0})".format(args.STRIP_MASK)
			obj_ndx.add_def("Strip", mask).output_ndx(gmx_arg["-n"])
			gmx_eof = "<< 'EOF'\nStrip\nEOF"

		obj_ndx.output_ndx(gmx_arg["-n"])

		command = " ".join([command_gmx, "trjconv"] + ["{0} {1}".format(o, v) for o, v in gmx_arg.items() if v is not None])
		command += " " + gmx_eof
		exec_sp(command, True)


		# treat cluster in PBC (Molecular collisions occur)
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate trajectory with adjusted molecules pairs."), LOG_COLOR, attrs=["bold"]))
		del(gmx_arg["-b"], gmx_arg["-e"], gmx_arg["-skip"])
		gmx_arg["-pbc"] = "cluster"

		gmx_arg["-f"] = gmx_arg["-o"]
		gmx_arg["-o"] = tempfile_name_full + "_step2_cluster.trr"
		if ambertools_ver >= 17:
			gmx_arg["-o"] = tempfile_name_full + "_step2_cluster.xtc"
		delete_files.append(gmx_arg["-o"])

		if args.STRIP_MASK is None:
			gmx_eof = "<< 'EOF'\nCenter\nSystem\nEOF"
		else:
			gmx_eof = "<< 'EOF'\nCenter\nStrip\nEOF"

		command = " ".join([command_gmx, "trjconv"] + ["{0} {1}".format(o, v) for o, v in gmx_arg.items() if v is not None])
		command += " " + gmx_eof
		exec_sp(command, True)


		# remove collision
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate trajectory with molecular collisions removed."), LOG_COLOR, attrs=["bold"]))
		gmx_arg["-pbc"] = "mol"
		gmx_arg["-ur"] = "compact"

		gmx_arg["-f"] = gmx_arg["-o"]
		gmx_arg["-o"] = tempfile_name_full + "_step3_mol.trr"
		if ambertools_ver >= 17:
			gmx_arg["-o"] = tempfile_name_full + "_step3_mol.xtc"
		delete_files.append(gmx_arg["-o"])
		trajectory_input = gmx_arg["-o"]

		gmx_arg["-center"] = ""
		if args.STRIP_MASK is None:
			gmx_eof = "<< 'EOF'\nCenter\nSystem\nEOF"
		else:
			gmx_eof = "<< 'EOF'\nCenter\nStrip\nEOF"

		command = " ".join([command_gmx, "trjconv"] + ["{0} {1}".format(o, v) for o, v in gmx_arg.items() if v is not None])
		command += " " + gmx_eof
		exec_sp(command, True)


		# output .gro
		if os.path.splitext(args.OUTPUT_FILE)[1].lower() == ".gro":
			process_i += 1
			sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate trajectory with molecular collisions removed."), LOG_COLOR, attrs=["bold"]))

			check_overwrite(args.OUTPUT_FILE, args.FLAG_OVERWRITE)
			del(gmx_arg["-pbc"])
			gmx_arg["-f"] = gmx_arg["-o"]
			gmx_arg["-o"] = args.OUTPUT_FILE
			gmx_arg["-fit"] = "rot+trans"
			if args.STRIP_MASK is None:
				gmx_eof = "<< 'EOF'\nCenter\nCenter\nSystem\nEOF"
			else:
				gmx_eof = "<< 'EOF'\nCenter\nCenter\nStrip\nEOF"

			command = " ".join([command_gmx, "trjconv"] + ["{0} {1}".format(o, v) for o, v in gmx_arg.items() if v is not None])
			command += " " + gmx_eof
			exec_sp(command, True)

			# delete temporary files
			delete_all()
			sys.exit(0)


	# final conversion (rot+trans)
	check_overwrite(args.OUTPUT_FILE, args.FLAG_OVERWRITE)

	temp_in = tempfile_name_full + ".in"
	delete_files.append(temp_in)

	process_i += 1
	sys.stdout.write(colored("Process ({0}/{1}): {2} => {3}\n".format(process_i, max_process, "Generate trajectory with rotated and shifted molecules.", args.OUTPUT_FILE), LOG_COLOR, attrs=["bold"]))
	with open(temp_in, "w") as obj_output:
		obj_output.write("parm {0}\n".format(args.PRMTOP_FILE))
		obj_output.write("trajin {0}\n".format(trajectory_input))
		obj_output.write("unwrap {0}\n".format(args.CENTER_MASK))
		obj_output.write("center {0} mass origin\n".format(args.CENTER_MASK))
		obj_output.write("rms {0} first mass\n".format(args.CENTER_MASK))
		obj_output.write("autoimage\n")
		obj_output.write("trajout {0}\n".format(args.OUTPUT_FILE))
		obj_output.write("go\n")
	exec_sp("{0} -i {1}".format(command_cpptraj, temp_in), True)


	# delete temporary files
	delete_all()
