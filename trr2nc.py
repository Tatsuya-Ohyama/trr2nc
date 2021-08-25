#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
trr2nc
Convert Gromacs trajectory to AMBER trajectory
"""

import sys, os, signal
import argparse
import subprocess
import tempfile
from termcolor import colored
import netCDF4
from tqdm import tqdm

from mods.func_prompt_io import *
from mods.molecule_topology import MoleculeTopology
from mods.file_NDX import FileNDX
from mods.file_GRO import FileGRO

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



# =============== variable =============== #
COMMAND_NAME_GMX = "gmx"
COMMAND_NAME_CPPTRAJ = "cpptraj"
LOG_COLOR = "yellow"



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
	parser = argparse.ArgumentParser(description = "trr2nc.py - Convert trr to nc with treating PBC", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-s", dest="TPR_FILE", metavar="INPUT.tpr", required=True, help="Gromacs run input file")
	parser.add_argument("-x", dest="TRR_FILE", metavar="INPUT.<trr|xtc|gro>", required=True, help="Gromacs trajectory file")
	parser.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.<nc|mdcrd>", required=True, help="output for Amber trajectory (.nc)")
	parser.add_argument("-t", dest="TOP_FILE", metavar="INPUT.top", required=True, help="Gromacs topology file")
	parser.add_argument("-p", dest="PRMTOP_FILE", metavar="OUTPUT.prmtop", required=True, help="Amber topology file")
	parser.add_argument("-sc", dest="TEMP_DIR", metavar="TEMP_DIR", default=".", help="Temporary directory (Default: current dir)")
	parser.add_argument("--separate-mol", dest="SEPARATE_MOL", metavar="MOL_NAME", nargs="+", default=[], help="separate molecules into individual molecules (periodic boundary condition problem)")

	gmx_option = parser.add_argument_group("gromacs option")
	gmx_option.add_argument("-b", dest="BEGIN", metavar="START_TIME", type=int, help="First frame index to read from trajectory (start from 0)")
	gmx_option.add_argument("-e", dest="END", metavar="END_TIME", type=int, help="Last frame index to read from trajectory (start from 0)")
	gmx_option.add_argument("--offset", dest="OFFSET", metavar="OFFSET", type=int, default=1, help="Output interval (Default: 1)")

	cpptraj_option = parser.add_argument_group("cpptraj option")
	cpptraj_option.add_argument("-mc", dest="CENTER_MASK", metavar="CENTER_MASK", help="center mask for cpptraj")
	cpptraj_option.add_argument("-ms", dest="STRIP_MASK", metavar="STRIP_MASK", help="strip mask for cpptraj")

	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")

	args = parser.parse_args()


	# check arguments
	check_exist(args.TPR_FILE, 2)
	check_exist(args.TRR_FILE, 2)
	check_exist(args.TOP_FILE, 2)
	command_gmx = check_command(COMMAND_NAME_GMX)
	command_cpptraj = check_command(COMMAND_NAME_CPPTRAJ)


	# determine name of temporary file
	tempfile_name = ""
	tempfile_name_cd = ""
	with tempfile.NamedTemporaryFile(mode="w", prefix=".trr2nc_", dir=".") as obj_output:
		tempfile_name_cd = os.path.basename(obj_output.name)
		tempfile_name = tempfile_name_cd
	tempfile_name = os.path.join(args.TEMP_DIR, tempfile_name)
	delete_files = []


	# when output is .gro file
	gro_temp = None
	if os.path.splitext(args.OUTPUT_FILE)[1] == ".gro":
		gro_temp = tempfile_name + ".gro"
		delete_files.append(gro_temp)


	# edit topology file
	top = args.TOP_FILE
	tpr = args.TPR_FILE
	if len(args.SEPARATE_MOL) != 0:
		sys.stderr.write(colored("INFO: Modify topology.\n", LOG_COLOR, attrs=["bold"]))
		sys.stderr.flush()
		top = tempfile_name_cd + ".top"
		delete_files.append(top)
		with open(top, "w") as obj_output:
			with open(args.TOP_FILE, "r") as obj_input:
				flag_read = False
				for line_val in obj_input:
					if flag_read:
						if len(line_val.strip()) == 0:
							continue

						data = line_val.strip().split()
						if data[0] in args.SEPARATE_MOL:
							sys.stderr.write(colored("    : molecule {0} is divided into {1} molecules.\n".format(data[0], data[1]), LOG_COLOR))
							for idx in range(int(data[1])):
								obj_output.write("{0:<5} {1} {2:>3}\n".format(data[0], " " * 10, 1))
						else:
							obj_output.write(line_val)

					elif "[ molecules ]" in line_val:
						flag_read = True
						obj_output.write(line_val)
					else:
						obj_output.write(line_val)

		mdp = tempfile_name + ".mdp"
		delete_files.append(mdp)
		with open(mdp, "w") as obj_output:
			obj_output.write("integrator = steep\n")
			obj_output.write("cutoff-scheme = Verlet\n")
			obj_output.write("nsteps = 1\n")
			obj_output.write("constraints = none\n")
			obj_output.write("coulombtype	= PME\n")

		tpr = tempfile_name + ".tpr"
		delete_files.append(tpr)
		mdout = tempfile_name + "_o.mdp"
		delete_files.append(mdout)
		exec_sp("{0} grompp -f {1} -c {2} -p {3} -o {4} -po {5}".format(command_gmx, mdp, args.TPR_FILE, top, tpr, mdout), False)


	# read topology file
	sys.stderr.write(colored("INFO: Loading topology file ... ", LOG_COLOR, attrs=["bold"]), )
	sys.stderr.flush()
	obj_top = MoleculeTopology(top)
	sys.stderr.write(colored("done\n", LOG_COLOR, attrs=["bold"]))


	# conversion of trajectory file (xtc, gro -> trr)
	trr_input = args.TRR_FILE
	nc_output = args.OUTPUT_FILE
	if os.path.splitext(args.TRR_FILE)[1] in [".gro", ".xtc"]:
		gmx_arg_opt = ["-s", "-f", "-o", "-b", "-e", "-skip", "-pbc", "-n", ""]
		gmx_arg_var = [tpr, args.TRR_FILE, None, None, None, args.OFFSET, "mol", None, None]
		trr_output = tempfile_name + ".trr"
		gmx_arg_var[2] = trr_output
		delete_files.append(trr_output)
		sys.stderr.write(colored("INFO: Start convert trajectory file ... \n    {0} -> {1}\n".format(trr_input, trr_output), LOG_COLOR, attrs=["bold"]))

		if args.BEGIN is not None:
			gmx_arg_var[3] = args.BEGIN
		if args.END is not None:
			gmx_arg_var[4] = args.END


		obj_ndx = FileNDX(obj_top)
		ndx = tempfile_name + ".ndx"
		delete_files.append(ndx)
		gmx_arg_var[7] = ndx
		if args.STRIP_MASK is None:
			# no strip_mask
			obj_ndx.output_ndx(ndx)
			gmx_arg_var[8] = "<< 'EOF'\n0\nEOF"
		else:
			# when strip_mask is specified, create .ndx file
			mask = "!({0})".format(args.STRIP_MASK)
			obj_ndx.add_def("Strip", mask).output_ndx(ndx)
			gmx_arg_var[8] = "<< 'EOF'\n1\nEOF"
			obj_top.set_mask(mask)

		command = " ".join([command_gmx, "trjconv"] + ["{0} {1}".format(o, v) for o, v in zip(gmx_arg_opt, gmx_arg_var) if v is not None])
		exec_sp(command, True)
		trr_input = trr_output

		# when output is .gro file
		if gro_temp is not None:
			nc_output = tempfile_name + ".nc"
			delete_files.append(nc_output)
			gmx_arg_var[2] = gro_temp
			gmx_arg_var[3] = 0
			gmx_arg_var[4] = 0
			command = " ".join([command_gmx, "trjconv"] + ["{0} {1}".format(o, v) for o, v in zip(gmx_arg_opt, gmx_arg_var) if v is not None])
			exec_sp(command, False)

		sys.stderr.write(colored("INFO: Finished creating trajectory file\n", LOG_COLOR, attrs=["bold"]))


	# create .prmtop file
	sys.stderr.write(colored("INFO: Creating prmtop ({file}) ... ".format(file=args.PRMTOP_FILE), LOG_COLOR, attrs=["bold"]))
	sys.stderr.flush()
	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.PRMTOP_FILE)
	obj_top.save_file(args.PRMTOP_FILE)
	sys.stderr.write(colored("done\n".format(file=args.PRMTOP_FILE), LOG_COLOR, attrs=["bold"]))


	# convert to .nc file
	if args.FLAG_OVERWRITE == False:
		check_overwrite(nc_output)
	temp_in = tempfile_name + ".in"
	delete_files.append(temp_in)
	sys.stderr.write(colored("INFO: Converting AMBER trajectory ({file})\n".format(file=args.OUTPUT_FILE), LOG_COLOR, attrs=["bold"]))
	with open(temp_in, "w") as obj_output:
		obj_output.write("parm {0}\n".format(args.PRMTOP_FILE))
		obj_output.write("trajin {0}\n".format(trr_input))
		if args.CENTER_MASK:
			obj_output.write("autoimage\n")
			obj_output.write("unwrap {0}\n".format(args.CENTER_MASK))
			obj_output.write("center {0} mass origin\n".format(args.CENTER_MASK))
			obj_output.write("rms {0} first mass\n".format(args.CENTER_MASK))
		obj_output.write("trajout {0}\n".format(nc_output))
		obj_output.write("go\n")
	exec_sp("{0} -i {1}".format(command_cpptraj, temp_in), True)


	# convert to .gro file
	if gro_temp is not None:
		sys.stderr.write(colored("INFO: Converting .nc to .gro. ... ", LOG_COLOR, attrs=["bold"]))
		sys.stderr.flush()
		obj_gro = FileGRO(gro_temp)
		trr = netCDF4.Dataset(nc_output)
		for frame_idx in tqdm(range(len(trr.dimensions["frame"])), ascii=True, leave=False):
			coord = [(x / 10).round(3).tolist() for x in trr.variables["coordinates"][frame_idx]]
			if frame_idx == 0:
				obj_gro.set_data(frame_idx, "coord", coord)
			else:
				obj_gro.add_frame(coord)

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)
		obj_gro.output_gro(args.OUTPUT_FILE)
		sys.stderr.write(colored("done.\n", LOG_COLOR, attrs=["bold"]))


	# delete temporary files
	delete_all()
