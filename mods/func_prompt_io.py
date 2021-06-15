#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
prompt I/O function module
"""

import sys
import os



# =============== function =============== #
def check_overwrite(file):
	"""
	Function to check for overwrite

	Args:
		file(str): target file path

	Returns:
		None
	"""
	if os.path.exists(file):
		# Warn if file exists
		sys.stderr.write("WARN: %s exists. Overwrite it? (y/N): " % file)
		sys.stderr.flush()
		user = sys.stdin.readline().strip()

		if user != "y":
			# When there is no permission to overwrite
			sys.exit(0)
		else:
			# If there is permission to overwrite, delete it in advance considering conflicts
			os.remove(file)


def check_exist(path, mode, flag_exit=True):
	"""
	Function to check for file existence

	Args:
		path(str): target file path
		mode(int): 1(existence) / 2(existence for file) / 3(existence for dir)
		flag_exit(bool): Exit if not present (Default: True)

	Returns:
		(bool) or exit(None)
	"""
	if path is None:
		sys.stderr.write("ERROR: Path is not specified (None value).\n")
		if flag_exit:
			sys.exit(1)
		else:
			return False
	if mode == 1:
		if not os.path.exists(path):
			sys.stderr.write("ERROR: No such path (%s)\n" % path)
			if flag_exit:
				sys.exit(1)
			else:
				return False
	elif mode == 2:
		if not os.path.isfile(path):
			sys.stderr.write("ERROR: No such file (%s)\n" % path)
			if flag_exit:
				sys.exit(1)
			else:
				return False
	elif mode == 3:
		if not os.path.isdir(path):
			sys.stderr.write("ERROR: No such directory (%s)\n" % path)
			if flag_exit:
				sys.exit(1)
			else:
				return False
	else:
		sys.stderr.write("ERROR: Subroutine error: Not specified mode\n")
		sys.exit(1)
	return True
