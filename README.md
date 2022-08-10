# trr2nc

## 概要
Gromacs のトラジェクトリファイルを変換するプログラム


## 使用方法
```sh
$ trr2nc.py [-h] -s INPUT.tpr -x INPUT.<trr|xtc|gro> -o OUTPUT.<nc|mdcrd|xtc|pdb> -t INPUT.top -p OUTPUT.prmtop [-sc TEMP_DIR] [--gmx COMMAND_GMX] [-b START_TIME] [-e END_TIME] [-skip OFFSET] [-tu TIME_UNIT] [--separate-mol MOL_NAME [MOL_NAME ...]] [--cpptraj COMMAND_CPPTRAJ] -mc CENTER_MASK [-ms STRIP_MASK] [--multi] [--leave-atom LEAVE_ATOM_MASK] [-O] [--keep]
```

* Basic options:
	* `-h`, `--help`
		: ヘルプメッセージを表示して終了する。
	* `-s INPUT.tpr`
		: Gromacs の .tpr ファイル (Input)
	* `-x INPUT.<trr|xtc|gro>`
		: Gromacs のトラジェクトリファイル (Input)
	* `-o OUTPUT.<nc|mdcrd|xtc>`
		: Amber のトラジェクトリファイル (Output)
	* `-t INPUT.top`
		: Gromacs のトポロジーファイル (Input)
	* `-p OUTPUT.prmtop`
		: Amber のトポロジーファイル (Output)
	* `-sc TEMP_DIR`
		: 一時ファイルを置くためのディレクトリ (Default: カレントディレクトリ)
	* `-O`
		: プロンプトを出さずに上書きする。
	* `--keep`
		: 中間生成ファイルを残す。

* Gromacs option:
	* `--gmx COMMAND_GMX`
		: `gmx` コマンドパス (Default: 自動検出)
	* `-b START_TIME`
		: 読み込み開始時間 (ps) (start from 0)
	* `-e END_TIME`
		: 読み込み終了時間 (ps) (start from 0)
	* `-skip OFFSET`
		: 出力するフレームの間隔 (Default: 1)
	* `-tu`
		: 時間の単位 (Default: ps)
	* `--separate-mol MOL_NAME [MOL_NAME ...]`
		: .top ファイル内の `[ molecules ]` の分子を 1 つずつ分割するための分子名 (周期境界条件対策)

* cpptraj option:
	* `--cpptraj COMMAND_CPPTRAJ`
		: `cpptraj` コマンドパス (Default: 自動検出)
	* `-mc CENTER_MASK`
		: トラジェクトリの中心に配置する原子群の Ambermask
	* `-ms STRIP_MASK`
		: トラジェクトリから削除する原子群の Ambermask
	* `--multi`
		: 各フレーム毎に .pdb ファイルに出力する。
	* `--leave-atom`
		: 残す原子の Amber mask (生体分子から一定距離の水分子の切り出し等で使用する。出力は .pdb ファイルのみ使用可。例: `:1-20<:5.0`)


## pdb_separator.py
トラジェクトリを .pdb ファイルに変換する際に誤って一つのファイルにまとめてしまった (`--multi` オプションを付け忘れた) 場合の救済プログラム

### 使用方法
```sh
$ pdb_separator.py [-h] -i INPUT.pdb -o OUTPUT.pdb [-O]
```

* `-h`, `--help`
	: ヘルプメッセージを表示して終了する。
* `-i INPUT.pdb`
	: 入力ファイル
* `-o OUTPUT.pdb`
	: 出力ファイルの接頭辞
* `-O`
	: プロンプトを出さずに上書きする。


## 動作要件
* Python3
	* numpy
	* parmed
	* termcolor


## License
This software is released under [the MIT License](https://opensource.org/licenses/mit-license.php).


## Authors
* Tatsuya Ohyama


## ChangeLog
### Ver. 19.6 (2022-08-10)
* ボックス内の溶媒分子の座標がおかしくなるバグを修正した。

### Ver. 19.5 (2022-02-08)
* `--keep` オプションを追加した。
* `--leave-atom` オプションを追加した。

### Ver. 19.4 (2022-01-25)
* `-tu` オプションを追加した。
* コマンドライン引数の `-b` と `-e` の説明を変更した。

### Ver. 19.3 (2022-01-25)
* `pdb_separator.py` を追加した。

### Ver. 19.2 (2022-01-24)
* `--multi` オプションを追加した。

### Ver. 19.1 (2022-01-24)
* 中間ファイル `mdout.mdp` が残るバグを修正した。

### Ver. 19.0 (2021-12-06)
* `-ms` オプションを指定した際に構造が崩れる致命的なバグを修正した。

### Ver. 18.0 (2021-10-27)
* `cpptraj` と `gmx` コマンドを指定できるようにした。
* .xtc 出力時に構造がおかしくなる致命的なバグを修正した。

### Ver. 17.10 (2021-10-12)
* AmberTools 17 以降のインストール環境で、.xtc ファイルの出力および中間ファイルとして .xtc ファイルを出力するようにした。
* 使用モジュールを整理した。

### Ver. 17.9 (2021-07-14)
* モジュール内のスタイルも PEP8 スタイルに変更した。

### Ver. 17.8 (2021-07-14)
* PEP8 スタイルに変更した。
* mercurial から git に変更した。
* `README.md` を追加した。
