#!/usr/bin/env python3

import os
import gzip
import glob
import re
import argparse
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
import concurrent.futures
import itertools
from datetime import datetime

# ------------------------- 默认配置 -------------------------
DEFAULT_BARCODE_PATHS = {
    'T5': Path('/share/home/sxzhang25/01.scMulit/02.e_ca_QTL/script/index/T5'),
    'T7': Path('/share/home/sxzhang25/01.scMulit/02.e_ca_QTL/script/index/T7'),
    'I5': Path('/share/home/sxzhang25/01.scMulit/02.e_ca_QTL/script/index/I5'),
    'I7': Path('/share/home/sxzhang25/01.scMulit/02.e_ca_QTL/script/index/I7')
}

# ------------------------- 核心函数 -------------------------
@lru_cache(maxsize=1024)
def parse_range(range_str: str, range_type: str) -> list:
    """
    解析范围表达式并缓存结果
    """
    expanded = []
    try:
        sub_ranges = range_str.split(";")

        for sub in sub_ranges:
            if ":" in sub:
                start_str, end_str = sub.split(":")
                prefix_start, num_start = start_str.split("-")
                prefix_end, num_end = end_str.split("-")

                if prefix_start != prefix_end:
                    raise ValueError(f"区间前后缀不匹配: {sub}")

                for i in range(int(num_start), int(num_end) + 1):
                    expanded.append(f"{prefix_start}-{i}")
            else:
                expanded.append(sub)

        return expanded
    except Exception as e:
        raise ValueError(f"解析范围失败: {range_str}，错误信息: {str(e)}") from e

def load_barcode(file_path: Path) -> dict:
    """加载条码文件到字典"""
    if not file_path.exists():
        raise FileNotFoundError(f"条码文件不存在: {file_path}")

    barcode = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                barcode[parts[0]] = parts[1]
    return barcode

def generate_barcode_map(sample_file: Path, barcode_files: dict) -> dict:
    """生成barcode到样本的映射关系"""
    barcode_db = {
        'T5': load_barcode(barcode_files['T5']),
        'T7': load_barcode(barcode_files['T7']),
        'I5': load_barcode(barcode_files['I5']),
        'I7': load_barcode(barcode_files['I7'])
    }

    barcode_map = {}
    conflicts = []

    with open(sample_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 5:
                print(f"跳过第{line_num}行：列数不足")
                continue

            sample = parts[0]
            try:
                t5 = parse_range(parts[1], 'T5')
                t7 = parse_range(parts[2], 'T7')
                i5 = parse_range(parts[3], 'I5')
                i7 = parse_range(parts[4], 'I7')
            except ValueError as e:
                print(f"第{line_num}行错误: {str(e)}")
                continue

            for combo in itertools.product(t7, t5, i7, i5):
                try:
                    full_bc = (
                        barcode_db['T7'][combo[0]] +
                        barcode_db['T5'][combo[1]] +
                        barcode_db['I7'][combo[2]] +
                        barcode_db['I5'][combo[3]]
                    )
                except KeyError as e:
                    raise ValueError(f"缺失条码索引: {e} (样本: {sample})")

                if full_bc in barcode_map:
                    if barcode_map[full_bc] != sample:
                        conflicts.append( (full_bc, barcode_map[full_bc], sample) )
                    continue

                barcode_map[full_bc] = sample

    if conflicts:
        print("\n发现条码冲突：")
        for bc, s1, s2 in conflicts:
            print(f"条码 {bc} 被 {s1} 和 {s2} 共享")
        raise ValueError("存在冲突条码，请检查样本配置")

    return barcode_map

def process_single_file(args: tuple):
    """处理单个fastq文件对（并行执行）"""
    r1_path, r2_path, barcode_map, output_dir, r1_flag, ext = args
    handles = defaultdict(lambda: {'R1': None, 'R2': None, 'count': 0})

    total_input_reads = 0

    # 利用传入的 ext 和 r1_flag 精准提取前缀
    # 例如：原名 "sample_R1.fastq.gz"
    # 第一步去后缀：变成 "sample_R1"
    # 第二步去标识：变成 "sample"
    file_name = Path(r1_path).name
    input_prefix = file_name.replace(ext, '').replace(r1_flag, '')
    # 注意：如果你的标识是 "_R1_"，替换后可能会留下前后连字符的问题，

    try:
        with gzip.open(r1_path, 'rt') as r1_in, \
             gzip.open(r2_path, 'rt') as r2_in:

            while True:
                r1_header = r1_in.readline()
                if not r1_header:
                    break
                total_input_reads += 1
                r1_seq = r1_in.readline().strip()
                _ = r1_in.readline()
                r1_qual = r1_in.readline().strip()

                r2_header = r2_in.readline()
                r2_seq = r2_in.readline().strip()
                _ = r2_in.readline()
                r2_qual = r2_in.readline().strip()

                # 提取并组合barcode
                header_16bp = r1_header.split()[-1].split(':')[-1].replace('+', '')
                bc_part = r2_seq[:8] + r1_seq[:8]
                full_bc = bc_part + header_16bp

                # 匹配样本
                sample = barcode_map.get(full_bc)
                if not sample:
                    continue

                # 截取序列（去除前30bp）
                new_r1_seq = r1_seq[30:]
                new_r2_seq = r2_seq[30:]
                new_r1_qual = r1_qual[30:]
                new_r2_qual = r2_qual[30:]

                # 获取文件句柄，直接输出为压缩格式 (compresslevel=4 提升写入速度)
                if not handles[sample]['R1']:
                    out_r1_name = f"{sample}.{input_prefix}_R1.fastq.gz"
                    handles[sample]['R1'] = gzip.open(
                         Path(output_dir) / out_r1_name, 'at', compresslevel=4
                    )
                if not handles[sample]['R2']:
                    out_r2_name = f"{sample}.{input_prefix}_R2.fastq.gz"
                    handles[sample]['R2'] = gzip.open(
                        Path(output_dir) / out_r2_name, 'at', compresslevel=4
                    )

                # 写入数据（保持原版格式，不插入 UMI）
                handles[sample]['R1'].write(
                    f"@{full_bc}:{r1_header[1:]}{new_r1_seq}\n+\n{new_r1_qual}\n"
                )
                handles[sample]['R2'].write(
                    f"@{full_bc}:{r2_header[1:]}{new_r2_seq}\n+\n{new_r2_qual}\n"
                )
                handles[sample]['count'] += 1

    finally:
        # 关闭并保存所有压缩文件
        for h in handles.values():
            if h['R1']: h['R1'].close()
            if h['R2']: h['R2'].close()

    sample_counts = {sample: h['count'] for sample, h in handles.items() if h['count'] > 0}
    return {"input_reads": total_input_reads, "sample_counts": sample_counts}

# ------------------------- 命令行接口 -------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="单细胞测序数据分拣工具",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-i', '--input-dir', required=True, help="原始fastq文件目录")
    parser.add_argument('-o', '--output-dir', required=True, help="输出目录")
    parser.add_argument('-s', '--sample-file', required=True, help="样本配置文件路径")
    parser.add_argument('-t', '--threads', type=int, default=1, help="并行进程数")

    parser.add_argument('--r1-flag', default='_R1_', help="R1文件的特征标识，例如 '_R1_', '_1'")
    parser.add_argument('--r2-flag', default='_R2_', help="R2文件的特征标识，例如 '_R2_', '_2'")
    parser.add_argument('--ext', default='.fastq.gz', help="输入文件的扩展名，例如 '.fastq.gz', '.fq.gz'")

    for bc_type in ['T5', 'T7', 'I5', 'I7']:
        parser.add_argument(
            f'--{bc_type.lower()}',
            type=Path,
            default=DEFAULT_BARCODE_PATHS[bc_type],
            help=f"{bc_type}条码文件路径"
        )

    return parser.parse_args()

def validate_arguments(args):
    if not Path(args.input_dir).is_dir():
        raise NotADirectoryError(f"输入目录不存在: {args.input_dir}")
    var_file = Path(args.sample_file)
    if not var_file.is_file():
        raise FileNotFoundError(f"样本配置文件不存在: {var_file}")
    for bc_type in ['T5', 'T7', 'I5', 'I7']:
        bc_path = getattr(args, bc_type.lower())
        if not bc_path.exists():
            raise FileNotFoundError(f"{bc_type}条码文件不存在: {bc_path}")

def main():
    args = parse_arguments()

    try:
        validate_arguments(args)
    except Exception as e:
        print(f"参数错误: {str(e)}")
        return

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # 记录并打印开始时间
    start_time = datetime.now()
    print(f"\n[{start_time.strftime('%Y-%m-%d %H:%M:%S')}] 脚本开始运行...")

    try:
        print("正在生成barcode映射表...")
        barcode_files = {'T5': args.t5, 'T7': args.t7, 'I5': args.i5, 'I7': args.i7}
        barcode_map = generate_barcode_map(Path(args.sample_file), barcode_files)
        print(f"已加载 {len(barcode_map)} 个唯一barcode")
    except Exception as e:
        print(f"生成barcode映射失败: {str(e)}")
        return

    file_pairs = []
    search_pattern = f'*{args.r1_flag}*{args.ext}'
    print(f"正在输入目录中搜索匹配模式: {search_pattern}")

    for r1_path in glob.glob(str(Path(args.input_dir) / search_pattern)):
        p = Path(r1_path)
        if args.r1_flag not in p.name:
            print(f"警告：文件名中未找到R1标识 {p.name}")
            continue
        r2_path = str(p.with_name(p.name.replace(args.r1_flag, args.r2_flag, 1)))
        if Path(r2_path).exists():
            file_pairs.append( (r1_path, r2_path) )
        else:
            print(f"警告：缺失配对文件 {r2_path}")

    if not file_pairs:
        print("错误：未找到有效的fastq文件对")
        return

    print(f"开始处理 {len(file_pairs)} 个文件对（使用 {args.threads} 进程）...")

    global_total_input_reads = 0
    global_sample_counts = defaultdict(int)

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for r1, r2 in file_pairs:
            futures.append(
                executor.submit(
                    process_single_file,
                    (r1, r2, barcode_map, args.output_dir, args.r1_flag, args.ext)
                )
            )

        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                global_total_input_reads += result["input_reads"]

                for sample, count in result["sample_counts"].items():
                    global_sample_counts[sample] += count

                print(f"已处理 {result['input_reads']} 条原始reads...")
            except Exception as e:
                print(f"处理失败: {str(e)}")

    print("\n开始生成统计报告...")
    global_total_valid_reads = sum(global_sample_counts.values())

    # 获取脚本执行时的当前运行目录 (Current Working Directory)
    run_dir = Path.cwd()

    # 获取当前输出目录的名称（保留用来做文件名后缀，方便区分不同任务）
    dir_name = Path(args.output_dir).resolve().name

    # 组合新的路径和文件名：运行目录 / FQstat.输出目录名称
    stat_file_path = run_dir / f"FQstat.{dir_name}"

    with open(stat_file_path, "w") as stat_f:
        stat_f.write(f"总reads数\t{global_total_input_reads}\n")
        stat_f.write(f"有效reads数\t{global_total_valid_reads}\n")

        split_rate = (global_total_valid_reads / global_total_input_reads * 100) if global_total_input_reads > 0 else 0
        stat_f.write(f"拆分率\t{split_rate:.2f}%\n")

        for sample, count in sorted(global_sample_counts.items()):
            sample_pct = (count / global_total_valid_reads * 100) if global_total_valid_reads > 0 else 0
            stat_f.write(f"{sample}\t{count}\t{sample_pct:.2f}%\n")

    print(f"统计报告已保存至: {stat_file_path}")

    # 记录结束时间并计算总耗时
    end_time = datetime.now()
    elapsed_time = end_time - start_time

    print(f"\n[{end_time.strftime('%Y-%m-%d %H:%M:%S')}] 所有拆分与压缩任务已完成！")
    # 把 elapsed_time 转换为不带微秒的格式，看起来更清爽
    print(f"任务总耗时: {str(elapsed_time).split('.')[0]}")

if __name__ == "__main__":
    main()