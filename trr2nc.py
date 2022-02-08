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
MDP_FILE = """
integrator	= steep
emtol		= 1000.0
emstep		= 0.01
nsteps		= 50000
nstlist		= 100
ns_type		= grid
rlist		= 1.0
coulombtype	= PME
rcoulomb	= 1.0
nstlog = 1
pbc = xyz
vdwtype = cut-off
constraints = none
cutoff-scheme = Verlet
"""



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


def output_mdp(output_file):
	"""
	.tpr ファイル作成のためのダミー .mdp ファイルを作成する関数

	Args:
		output_file (str)): 出力先
	"""
	with open(output_file, "w") as obj_output:
		obj_output.write(MDP_FILE)



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "Program to convert Gromacs trajectory to AMBER trajectory", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-s", dest="TPR_FILE", metavar="INPUT.tpr", required=True, help="Gromacs run input file")
	parser.add_argument("-x", dest="TRAJECTORY_FILE", metavar="INPUT.<trr|xtc|gro>", required=True, help="Gromacs trajectory file")
	parser.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.<nc|mdcrd|xtc|pdb>", required=True, help="output trajectory")
	parser.add_argument("-t", dest="TOP_FILE", metavar="INPUT.top", required=True, help="Gromacs topology file")
	parser.add_argument("-p", dest="PRMTOP_FILE", metavar="OUTPUT.prmtop", required=True, help="Amber topology file")
	parser.add_argument("-sc", dest="TEMP_DIR", metavar="TEMP_DIR", default=".", help="Temporary directory (Default: current dir)")
	parser.add_argument("--separate-mol", dest="SEPARATE_MOL", metavar="MOL_NAME", nargs="+", default=[], help="separate molecules into individual molecules (specify molecule name written in .top file) (periodic boundary condition problem)")

	gmx_option = parser.add_argument_group("gromacs option")
	gmx_option.add_argument("-b", dest="BEGIN", metavar="START_TIME", type=int, help="First frame index to read from trajectory (ps) (start from 0)")
	gmx_option.add_argument("-e", dest="END", metavar="END_TIME", type=int, help="Last frame index to read from trajectory (ps) (start from 0)")
	gmx_option.add_argument("-skip", dest="OFFSET", metavar="OFFSET", type=int, default=1, help="Only write every nr-th frame (Default: 1)")
	gmx_option.add_argument("-tu", dest="TIME_UNIT", metavar="TIME_UNIT", default="ps", choices=["fs", "ps", "ns", "us", "ms", "s"], help="Unit for time values: fs, ps, ns, us, ms, s (Default: ps)")
	gmx_option.add_argument("--gmx", dest="COMMAND_GMX", metavar="COMMAND_GMX", help="command line path for `gmx` (Default: autodetect)")

	cpptraj_option = parser.add_argument_group("cpptraj option")
	cpptraj_option.add_argument("-mc", dest="CENTER_MASK", metavar="CENTER_MASK", required=True, help="center mask for cpptraj")
	cpptraj_option.add_argument("-ms", dest="STRIP_MASK", metavar="STRIP_MASK", help="strip mask for cpptraj")
	cpptraj_option.add_argument("--cpptraj", dest="COMMAND_CPPTRAJ", metavar="COMMAND_CPPTRAJ", help="command line path for `cpptraj` (Default: autodetect)")
	cpptraj_option.add_argument("--multi", dest="FLAG_MULTI", action="store_true", default=False, help="Output PDB file for each frame")

	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	parser.add_argument("--keep", dest="FLAG_KEEP", action="store_true", default=False, help="Leave intermediate files")

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
	elif os.path.splitext(args.OUTPUT_FILE)[1].lower() == ".gro":
		max_process = 11
	else:
		max_process = 12


	# read topology file
	process_i += 1
	sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Loading topology file"), LOG_COLOR, attrs=["bold"]))
	obj_topol = parmed.gromacs.GromacsTopologyFile(args.TOP_FILE)


	# conversion of trajectory file
	trajectory_input = args.TRAJECTORY_FILE
	if os.path.splitext(args.TRAJECTORY_FILE)[1].lower() in [".gro", ".xtc"]:

		# create .ndx file
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate initial .ndx file."), LOG_COLOR, attrs=["bold"]))
		obj_ndx1 = FileNDX(obj_topol)
		obj_ndx1.add_def("Center", args.CENTER_MASK)
		ndx_file1 = tempfile_name_full + "1.ndx"
		gmx_eof = None
		if args.STRIP_MASK is None:
			# no strip_mask
			# output
			gmx_eof = "<< 'EOF'\nSystem\nEOF"
		else:
			# when strip_mask is specified, create .ndx file
			mask = "!({0})".format(args.STRIP_MASK)
			obj_ndx1.add_def("Strip", mask).output_ndx(ndx_file1)

			# output
			gmx_eof = "<< 'EOF'\nStrip\nEOF"
		obj_ndx1.output_ndx(ndx_file1)
		if not args.FLAG_KEEP:
			delete_files.append(ndx_file1)


		# create trajectory file with treating PBC
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate trajectory with adjusted molecules across the boundary."), LOG_COLOR, attrs=["bold"]))
		gmx_arg = {
			"-s": args.TPR_FILE,
			"-f": args.TRAJECTORY_FILE,
			"-o": None,
			"-b": args.BEGIN,
			"-e": args.END,
			"-n": ndx_file1,
			"-skip": args.OFFSET,
			"-tu": args.TIME_UNIT,
			"-pbc": "whole",
		}
		step1_whole_trajectory = tempfile_name_full + "_step1_whole.trr"
		if ambertools_ver >= 17:
			step1_whole_trajectory = tempfile_name_full + "_step1_whole.xtc"
		gmx_arg["-o"] = step1_whole_trajectory

		command = " ".join([command_gmx, "trjconv"] + ["{0} {1}".format(o, v) for o, v in gmx_arg.items() if v is not None])
		command += " " + gmx_eof
		exec_sp(command, True)
		if not args.FLAG_KEEP:
			delete_files.append(step1_whole_trajectory)


		# create .gro file for new .tpr file
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate stripped .gro file."), LOG_COLOR, attrs=["bold"]))
		tmp_gro_file = tempfile_name_full + "_tmp.gro"
		gmx_arg = {
			"-s": args.TPR_FILE,
			"-f": args.TRAJECTORY_FILE,
			"-o": tmp_gro_file,
			"-b": 0,
			"-e": 0,
			"-n": ndx_file1,
		}
		command = " ".join([command_gmx, "trjconv"] + ["{0} {1}".format(o, v) for o, v in gmx_arg.items() if v is not None])
		command += " " + gmx_eof
		exec_sp(command, False)
		if not args.FLAG_KEEP:
			delete_files.append(tmp_gro_file)


		# create .top file
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate stripped .top file."), LOG_COLOR, attrs=["bold"]))
		top_file = tempfile_name_full + ".top"
		if args.STRIP_MASK is not None:
			obj_topol.strip(args.STRIP_MASK)
		obj_topol.save(top_file)
		if not args.FLAG_KEEP:
			delete_files.append(top_file)


		# create stripped .ndx file
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate stripped .ndx file."), LOG_COLOR, attrs=["bold"]))
		obj_ndx2 = FileNDX(obj_topol)
		obj_ndx2.add_def("Center", args.CENTER_MASK)
		ndx_file2 = tempfile_name_full + "2.ndx"
		obj_ndx2.output_ndx(ndx_file2)
		if not args.FLAG_KEEP:
			delete_files.append(ndx_file2)


		# create .mdp file
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate stripped .mdp file."), LOG_COLOR, attrs=["bold"]))
		mdp_file = tempfile_name_full + ".mdp"
		output_mdp(mdp_file)
		if not args.FLAG_KEEP:
			delete_files.append(mdp_file)


		# create stripped .tpr file
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate stripped .tpr file."), LOG_COLOR, attrs=["bold"]))
		tpr_file = tempfile_name_full + ".tpr"
		tmp_mdp_file = tempfile_name_full + "_out.mdp"
		gmx_arg = {
			"-f": mdp_file,
			"-c": tmp_gro_file,
			"-o": tpr_file,
			"-p": top_file,
			"-maxwarn": 100,
			"-po": tmp_mdp_file,
		}
		command = " ".join([command_gmx, "grompp"] + ["{0} {1}".format(o, v) for o, v in gmx_arg.items() if v is not None])
		command += " " + gmx_eof
		exec_sp(command, True)
		if not args.FLAG_KEEP:
			delete_files.append(tpr_file)
			delete_files.append(tmp_mdp_file)


		# create trajectory file with treating cluster in PBC (Molecular collisions occur)
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate trajectory with adjusted molecules pairs."), LOG_COLOR, attrs=["bold"]))
		gmx_arg = {
			"-s": tpr_file,
			"-f": step1_whole_trajectory,
			"-o": None,
			"-n": ndx_file2,
			"-pbc": "cluster",
		}
		step2_cluster_trajectory = tempfile_name_full + "_step2_cluster.trr"
		if ambertools_ver >= 17:
			step2_cluster_trajectory = tempfile_name_full + "_step2_cluster.xtc"
		gmx_arg["-o"] = step2_cluster_trajectory
		gmx_eof = "<< 'EOF'\nCenter\nSystem\nEOF"
		command = " ".join([command_gmx, "trjconv"] + ["{0} {1}".format(o, v) for o, v in gmx_arg.items() if v is not None])
		command += " " + gmx_eof
		exec_sp(command, True)
		if not args.FLAG_KEEP:
			delete_files.append(step2_cluster_trajectory)


		# remove collision
		process_i += 1
		sys.stdout.write(colored("Process ({0}/{1}): {2}\n".format(process_i, max_process, "Generate trajectory with molecular collisions removed."), LOG_COLOR, attrs=["bold"]))
		gmx_arg["-pbc"] = "mol"
		gmx_arg["-ur"] = "compact"
		gmx_arg["-center"] = ""

		gmx_arg["-f"] = gmx_arg["-o"]
		step3_mol_trajectory = tempfile_name_full + "_step3_mol.trr"
		if ambertools_ver >= 17:
			step3_mol_trajectory = tempfile_name_full + "_step3_mol.xtc"
		gmx_arg["-o"] = step3_mol_trajectory
		trajectory_input = gmx_arg["-o"]

		gmx_eof = "<< 'EOF'\nCenter\nSystem\nEOF"
		command = " ".join([command_gmx, "trjconv"] + ["{0} {1}".format(o, v) for o, v in gmx_arg.items() if v is not None])
		command += " " + gmx_eof
		exec_sp(command, True)
		if not args.FLAG_KEEP:
			delete_files.append(gmx_arg["-o"])


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


	# create .prmtop
	process_i += 1
	sys.stdout.write(colored("Process ({0}/{1}): {2} => {3}\n".format(process_i, max_process, "Generate prmtop", args.PRMTOP_FILE), LOG_COLOR, attrs=["bold"]))
	check_overwrite(args.PRMTOP_FILE, args.FLAG_OVERWRITE)
	obj_topol.save(args.PRMTOP_FILE)


	# final conversion (rot+trans)
	check_overwrite(args.OUTPUT_FILE, args.FLAG_OVERWRITE)

	process_i += 1
	sys.stdout.write(colored("Process ({0}/{1}): {2} => {3}\n".format(process_i, max_process, "Generate trajectory with rotated and shifted molecules.", args.OUTPUT_FILE), LOG_COLOR, attrs=["bold"]))
	temp_in = tempfile_name_full + ".in"
	with open(temp_in, "w") as obj_output:
		obj_output.write("parm {0}\n".format(args.PRMTOP_FILE))
		obj_output.write("trajin {0}\n".format(trajectory_input))
		obj_output.write("unwrap {0}\n".format(args.CENTER_MASK))
		obj_output.write("center {0} mass origin\n".format(args.CENTER_MASK))
		obj_output.write("rms {0} first mass\n".format(args.CENTER_MASK))
		obj_output.write("autoimage\n")
		if os.path.splitext(args.OUTPUT_FILE)[1].lower() == ".pdb" and args.FLAG_MULTI:
			obj_output.write("trajout {0} multi\n".format(args.OUTPUT_FILE))
		else:
			obj_output.write("trajout {0}\n".format(args.OUTPUT_FILE))
		obj_output.write("go\n")

	if not args.FLAG_KEEP:
		delete_files.append(temp_in)
	exec_sp("{0} -i {1}".format(command_cpptraj, temp_in), True)


	# delete temporary files
	delete_all()
