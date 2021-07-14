# trr2nc

## 概要
Gromacs のトラジェクトリファイルを Amber のトラジェクトリファイルに変換するプログラム


## 使用方法
```sh
$ trr2nc.py [-h] -s INPUT.tpr -x INPUT.<trr|xtc|gro> -o OUTPUT.<nc|mdcrd> -t INPUT.top -p OUTPUT.prmtop [-sc TEMP_DIR] [--separate-mol MOL_NAME [MOL_NAME ...]] [-b START_TIME] [-e END_TIME] [--offset OFFSET] [-mc CENTER_MASK] [-ms STRIP_MASK] [-O]
```

* optional arguments:
	* `-h`, `--help`
		: ヘルプメッセージを表示して終了する。
	* `-s INPUT.tpr`
		: Gromacs の .tpr ファイル (Input)
	* `-x INPUT.<trr|xtc|gro>`
		: Gromacs のトラジェクトリファイル (Input)
	* `-o OUTPUT.<nc|mdcrd>`
		: Amber のトラジェクトリファイル (Output)
	* `-t INPUT.top`
		: Gromacs のトポロジーファイル (Input)
	* `-p OUTPUT.prmtop`
		: Amber のトポロジーファイル (Output)
	* `-sc TEMP_DIR`
		: 一時ファイルを置くためのディレクトリ (Default: カレントディレクトリ)
	* `--separate-mol MOL_NAME [MOL_NAME ...]`
		: .top ファイル内の `[ molecules ]` の分子を 1 つずつ分割するための分子名 (周期境界条件対策)
	* `-O`
		: プロンプトを出さずに上書きする。

* gromacs option:
	* `-b START_TIME`
		: 読み込み開始フレームインデックス (start from 0)
	* `-e END_TIME`
		: 読み込み終了フレームインデックス (start from 0)
	* `--offset OFFSET`
		: 出力するフレームの間隔 (Default: 1)

* cpptraj option:
	* `-mc CENTER_MASK`
		: トラジェクトリの中心に配置する原子群の Ambermask
	* `-ms STRIP_MASK`
		: トラジェクトリから削除する原子群の Ambermask


## 動作要件
* Python3
	* netCDF4
	* numpy
	* parmed
	* pathlib
	* pickle
	* termcolor
	* tqdm


## License
This software is released under the MIT License, see LICENSE.


## Authors
* Tatsuya Ohyama


## ChangeLog
* Ver. 17.9 (2021-07-14)
	* モジュール内のスタイルも PEP8 スタイルに変更した。
* Ver. 17.8 (2021-07-14)
	* PEP8 スタイルに変更した。
	* mercurial から git に変更した。
	* `README.md` を追加した。
